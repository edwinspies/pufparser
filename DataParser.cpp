#include "DataParser.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <list>
#include <algorithm>
#include <sys/stat.h>
#include <tuple>

using namespace std;

/// @brief Delivers the index of the std csv file; can be modified to read other csv structured puf data
enum CSV_STRUCTURE {
  SAMPLE_NUMBER,
  DATE,
  BOARD_TYPE,
  BOARD_ID,
  UNDEFINED_COLUMN_4,
  ADDRESS,
  NO_OF_BITS,
  UNDEFINED_COLUMN_7,
  UNDEFINED_COLUMN_8,
  RAWDATA
};

/// @brief Defines the index for the sample tuple
enum TUPLE_STRUCTURE {
  TUPLE_SAMPLE_NUMBER,
  TUPLE_BOARD_ID,
  TUPLE_RAWDATA
};

DataParser::DataParser(const string &fileName) {
  getDataFromCSV(fileName);
}

tuple<int, string, string> DataParser::getNextLineAndSplitIntoTokens(istream &str) {
  vector<string> result;
  tuple<int, string, string> returnTriple;
  size_t pos{};
  int sampleNumber;

  string line;
  getline(str, line);

  stringstream lineStream(line);
  string cell;

  while (getline(lineStream, cell, ';')) {
	result.push_back(cell);
  }
  // This checks for a trailing comma with no data after it.
  if (!lineStream && cell.empty()) {
	// If there was a trailing comma then add an empty element.
	result.emplace_back("");
  }

  //if an input line (for example the last) is empty or shorter than expected,
  // accessing the index at 9 would result in an error
  //TODO replace 9 with a dynamic value so other files can be processed safely
  if (result.size() >= 9) {
	try {
	  sampleNumber = stoi(result.at(CSV_STRUCTURE::SAMPLE_NUMBER), &pos);
	}
	catch (std::invalid_argument const &ex) {
	  cout << "std::invalid_argument::what(): " << ex.what() << '\n';
	  sampleNumber = -1;
	}

	get<TUPLE_SAMPLE_NUMBER>(returnTriple) = sampleNumber;
	get<TUPLE_BOARD_ID>(returnTriple) = result.at(CSV_STRUCTURE::BOARD_ID);
	get<TUPLE_RAWDATA>(returnTriple) = result.at(CSV_STRUCTURE::RAWDATA);
  }

  return returnTriple;
}

/// @brief Reads in a CSV that has
void DataParser::getDataFromCSV(const string &fileName) {
  filebuf fb;

  //if this function is called twice, clear all items before to remove old data
  p_listOfSamples.clear();

  if (fb.open(fileName, ios::in)) {
	istream is(&fb);
	while (is) {
	  p_listOfSamples.emplace_back(getNextLineAndSplitIntoTokens(is));
	}
	fb.close();
  }
  p_listOfSamples.pop_back();

}

/// @brief takes the list of all samples and proceeds to output them into the Python ND function format
void DataParser::processAndOutputDataToNDFormat() {
  cout << "This will take a while depending on the no. of samples" << endl;

  int check;
  check = mkdir("data", 0777);
  //is no dir was created, return
  if (check == -1 && errno != EEXIST) {
	cout << "data directory could not be created";
	return;
  }

  for (const tuple<int, string, string> &singleSample : p_listOfSamples) {
	writeDeviceDataIntoFile(singleSample);
  }
}

/// @brief receives a single device sample and processes this data to output a file for the Python ND Function
void DataParser::writeDeviceDataIntoFile(const tuple<int, string, string> &data) {
  //create dir
  int check;
  const string &dirname = "data/" + get<TUPLE_BOARD_ID>(data);

  check = mkdir(dirname.c_str(), 0777);
  //is no dir was created, return
  if (check == -1 && errno != EEXIST) {
	cout << "directory could not be created";
	return;
  }

  // Create and open a text file
  const string &filename = to_string(get<TUPLE_SAMPLE_NUMBER>(data));
  ofstream MyFile(dirname + "/" + filename);

  //process data bits into comma separated bits
  const string &deviceData = get<TUPLE_RAWDATA>(data);
  string commaSeparatedDeviceData = commaSeparateData(deviceData);

  // Write to the file
  MyFile << commaSeparatedDeviceData;

  // Close the file
  MyFile.close();
}

/// @brief Reads in the raw output of the PUF and converts it into comma separated data
string DataParser::commaSeparateData(const string &deviceData) {
  string returnString;

  int i = 0;
  while (deviceData.compare(i, 1, "\0")) {
	returnString += deviceData.at(i);
	if (deviceData.compare(i + 1, 1, "\0")) {
	  returnString += ",";
	}
	i++;
  }
  return returnString;
}

/// @param [in] data List that contains samples (ideally only from the same board ID)
/// @brief values close to 1 have a high tendency to be 1, analog for 0
double *DataParser::getProbabilityOfIndex(const list<tuple<int, string, string>> &samplesOfUniqueDevice) {
  int arraySize = 4096;

  auto *returnArray = (double *) (malloc(arraySize * sizeof(double)));

  for (tuple<int, string, string> sample : samplesOfUniqueDevice) {
	for (int i = 0; i < arraySize && get<TUPLE_RAWDATA>(sample).size(); ++i) {
	  if (get<TUPLE_RAWDATA>(sample).at(i) == '1') {
		returnArray[i]++;
	  }
	}
  }

  unsigned long sizeOfData = samplesOfUniqueDevice.size();
  for (int i = 0; i < arraySize; ++i) {
	returnArray[i] /= sizeOfData;
  }

  //return nullptr;
  return returnArray;
}
/// @brief Will write individual pixels in the file between 0-255,
/// depending on the corresponding values in the double array, ranging vom 0-1
void DataParser::outputGraph(double *array, const string &addFileName) {
  const string &filename = "picture";
  ofstream pictureFile(filename + "_" + addFileName + ".pgm");
  int possibleValues = 255;

  pictureFile << "P2" //Image format
			  << "\n"
			  << "64 64" //Image size; //TODO make this dynamic
			  << "\n"
			  << possibleValues //No. of possible values
			  << "\n";

  for (int i = 0; i < 4096; i++) {
	pictureFile << int(array[i] * possibleValues) << " ";
  }

  free(array);
}

/// @brief Returns all Samples that have a given board ID
list<tuple<int, string, string>> DataParser::extractSamplesByBoardID(const string &boardID) {
  list<tuple<int, string, string>> returnList;

  for (tuple<int, string, string> sample : p_listOfSamples) {
	if (get<TUPLE_BOARD_ID>(sample) == boardID) {
	  returnList.emplace_back(sample);
	}
  }

  return returnList;
}

///@brief Returns a set of (unique) board IDs that are located in the list of samples
set<string> DataParser::extractAllBoardIDs() {
  set<string> returnSet;
  for (tuple<int, string, string> sample : p_listOfSamples) {
	returnSet.insert(get<TUPLE_BOARD_ID>(sample));
  }

  return returnSet;
}






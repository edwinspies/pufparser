//
// Created by edwin on 22.11.2022.
//

#include "DataParser.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
#include <list>
#include <algorithm>
#include <sys/stat.h>
#include <tuple>

using namespace std;

enum structure {
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

//triple
//0 = sampleNumber
//1 = Board ID
//2 = Raw Data

DataParser::DataParser(const string &fileName) {
  getDataFromCSV(fileName);
}

tuple<int, string, string> DataParser::getNextLineAndSplitIntoTokens(istream &str) {
  vector<string> result;

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

  tuple<int, string, string> returnTriple;
  size_t pos{};
  int sampleNumber;

  //if an element (for example the last) is empty or shorter, accessing the index at 9 would result in an error
  if (result.size() >= 9) {
	try {
	  sampleNumber = stoi(result.at(0), &pos);
	}
	catch (std::invalid_argument const &ex) {
	  cout << "std::invalid_argument::what(): " << ex.what() << '\n';
	  sampleNumber = -1;
	}

	get<0>(returnTriple) = sampleNumber;
	get<1>(returnTriple) = result.at(3);
	get<2>(returnTriple) = result.at(9);
  }

  return returnTriple;
}

void DataParser::getDataFromCSV(const string &fileName) {
  filebuf fb;
  if (fb.open(fileName, ios::in)) {
	istream is(&fb);
	while (is) {
	  p_listOfSamples.emplace_back(getNextLineAndSplitIntoTokens(is));
	}
	fb.close();
  }
  p_listOfSamples.pop_back();

}

//takes all samples and processes them
void DataParser::processAndOutputDataToNDFormat() {
  for (const tuple<int, string, string>& singleSample : p_listOfSamples) {
	writeDeviceDataIntoFile(singleSample);
  }
}
//receives a single device sample and processes this data to output a files for the Python ND Function
void DataParser::writeDeviceDataIntoFile(const tuple<int, string, string>& data) {
  //create dir
  int check;
  const string& dirname = get<1>(data);

  check = mkdir(dirname.c_str(),0777);
  //is no dir was created, return
  if(check == -1 && errno!=EEXIST)  {
	cout << "directory could not be created";
	return;
  }

  // Create and open a text file
  const string& filename = to_string(get<0>(data));
  ofstream MyFile(dirname + "/" + filename);

  //process data bits into comma separated bits
  const string& deviceData = get<2>(data);
  string commaSeparatedDeviceData = commaSeparateData(deviceData);

  // Write to the file
  MyFile << commaSeparatedDeviceData;

  // Close the file
  MyFile.close();
}

string DataParser::commaSeparateData(const string& deviceData) {
  string returnString;

  int i = 0;
  while(deviceData.compare(i, 1, "\0")) {
	returnString+=deviceData.at(i);
	if(deviceData.compare(i+1, 1, "\0")) {
	  returnString+=",";
	}
	i++;
  }
  return returnString;
}

/// @param [in] data List that contains samples (ideally only from the same board ID)
/// @brief values close to 1 have a high tendency to be 1, analog for 0
double* DataParser::getProbabilityOfIndex(const list<tuple<int, string, string>>& samplesOfUniqueDevice) {
  int arraySize = 4096;

  auto *returnArray = (double *) (malloc(arraySize * sizeof(double)));

  for (tuple<int, string, string> sample : samplesOfUniqueDevice) {
	for (int i = 0; i < arraySize && get<2>(sample).size(); ++i) {
	  if (get<2>(sample).at(i) == '1') {
		returnArray[i]++;
	  }
	}
  }

  unsigned long sizeOfData = samplesOfUniqueDevice.size();
  for (int i = 0; i < arraySize; ++i) {
	returnArray[i]/=sizeOfData;
  }

  //return nullptr;
return returnArray;
}
/// @brief 1 will output white, 0 will output black; rest is in between
void DataParser::outputGraph(double *array, const string& addFileName) {
  const string& filename = "picture";
  ofstream pictureFile(filename + "_" + addFileName + ".pgm");

  pictureFile   << "P2"
  		   << "\n"
		   << "64 64"
		   << "\n"
		   << "255"
		   << "\n";

  for(int i = 0; i < 4096; i++) {
	pictureFile << int (array[i]*255) << " ";
  }

  free(array);
}
list<tuple<int, string, string>> DataParser::extractSamplesByBoardID(const string& boardID) {
  list<tuple<int, string, string>> returnList;

  for (tuple<int, string, string> sample : p_listOfSamples) {
	if (get<1>(sample) == boardID) {
	  returnList.emplace_back(sample);
	}
  }

  return returnList;
}
set<string> DataParser::extractAllBoardIDs() {
  set<string> returnSet;
  for (tuple<int, string, string> sample : p_listOfSamples) {
	returnSet.insert(get<1>(sample));
  }

  return returnSet;
}






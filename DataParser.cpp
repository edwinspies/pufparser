#include "DataParser.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <list>
#include <algorithm>
#include <sys/stat.h>
#include <tuple>
#include <bitset>

using namespace std;

DataParser::DataParser(const string &fileName, bool fileFormat) {
  altFileFormat = fileFormat;
  getDataFromCSV(fileName);
  currentSample = 0;
}

DataParser::~DataParser() {
  p_listOfSamples.clear();
}

/// @brief Reads in a CSV that has been saved in the format that is available at the Grenoble University; check enum CSV_STRUCTURE
void DataParser::getDataFromCSV(const string &fileName) {
  filebuf fb;

  //if this function is called twice, clear all items before to remove old data
  p_listOfSamples.clear();
  currentSample = 0;

  if (fb.open(fileName, ios::in)) {
	istream is(&fb);
	while (is) {
	  p_listOfSamples.emplace_back(getNextLineAndSplitIntoTokens(is));
	}
	fb.close();
  } else {
	cout << "File not found " + fileName;
  }

  if (altFileFormat) p_listOfSamples.pop_front();
  p_listOfSamples.pop_back();
}

/// @brief Reads a single line of the csv and splits it into a tuple
tuple<int, string, string, string> DataParser::getNextLineAndSplitIntoTokens(istream &str) {
  vector<string> result;
  tuple<int, string, string, string> returnTuple;

  string line;
  getline(str, line);

  stringstream lineStream(line);
  string cell;

  if (altFileFormat) {
	while (getline(lineStream, cell, '"')) {
	  result.push_back(cell);
	}
  } else {
	while (getline(lineStream, cell, ';')) {
	  result.push_back(cell);
	}
  }

  // This checks for a trailing comma with no data after it.
  if (!lineStream && cell.empty()) {
	// If there was a trailing comma then add an empty element.
	result.emplace_back("");
  }

  if (altFileFormat) {
	returnTuple = getTupleFromSergioCSVToken(result);
  } else {
	returnTuple = getTupleFromStdCSVToken(result);
  }

  return returnTuple;
}

tuple<int, string, string, string> DataParser::getTupleFromStdCSVToken(const vector<string> &result) {
  tuple<int, string, string, string> returnTuple;
  int sampleNumber;
  size_t pos{};

  //if an input line (for example the last) is empty or shorter than expected,
  // accessing the index at 9 would result in an error
  if (result.size() >= CSV_STRUCTURE::STD_CSV_LAST) {
	try {
	  sampleNumber = stoi(result.at(CSV_STRUCTURE::SAMPLE_NUMBER), &pos);
	}
	catch (std::invalid_argument const &ex) {
	  cout << "std::invalid_argument::what(): " << ex.what() << '\n';
	  sampleNumber = -1;
	}

	get<TUPLE_SAMPLE_NUMBER>(returnTuple) = sampleNumber;
	get<TUPLE_BOARD_ID>(returnTuple) = result.at(CSV_STRUCTURE::BOARD_ID);
	get<TUPLE_ADDRESS>(returnTuple) = result.at(CSV_STRUCTURE::ADDRESS);
	get<TUPLE_RAWDATA>(returnTuple) = result.at(CSV_STRUCTURE::RAWDATA);
  }

  return returnTuple;
}

tuple<int, string, string, string> DataParser::getTupleFromSergioCSVToken(const vector<string> &unsanitizedResult) {
  tuple<int, string, string, string> returnTuple;
  vector<string> result;
  int sampleNumber;
  size_t pos{};
  string rawData;

  for (const auto &i : unsanitizedResult) {
	if (i != "," && !i.empty()) {
	  result.emplace_back(i);
	}
  }

  //if an input line (for example the last) is empty or shorter than expected,
  // accessing the index at 9 would result in an error
  if (result.size() >= SERGIO_CSV_STRUCTURE::SERGIO_CSV_LAST) {
	sampleNumber = currentSample++;
	get<TUPLE_SAMPLE_NUMBER>(returnTuple) = sampleNumber;
	get<TUPLE_BOARD_ID>(returnTuple) = result.at(SERGIO_CSV_STRUCTURE::BOARD_ID_SERGIO);
	get<TUPLE_ADDRESS>(returnTuple) = result.at(SERGIO_CSV_STRUCTURE::ADDRESS_SERGIO);
	rawData = convertRawDataSergio(result.at(SERGIO_CSV_STRUCTURE::RAWDATA_SERGIO));
	get<TUPLE_RAWDATA>(returnTuple) = rawData;
  }

  return returnTuple;
}

string DataParser::convertRawDataSergio(const string &data) {
  string localData = data;
  std::string delimiter = ",";
  vector<string> tokens;
  string returnString;

  size_t pos = 0;
  std::string token;
  while ((pos = localData.find(delimiter)) != std::string::npos) {
	token = localData.substr(0, pos);
	tokens.emplace_back(token);
	localData.erase(0, pos + delimiter.length());
  }

  for (const string &s : tokens) {
	string binary = bitset<8>(stoi(s)).to_string(); //to binary
	returnString += binary;
  }

  return returnString;
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

  for (const tuple<int, string, string, string> &singleSample : p_listOfSamples) {
	writeDeviceDataIntoFile(singleSample);
  }
}

/// @param data Requires a single sample Tuple that will be written into a file
/// @brief receives a single device sample and processes this data to output a file for the Python ND Function
void DataParser::writeDeviceDataIntoFile(const tuple<int, string, string, string> &data) {
  //create dir
  int check;
  const string &dirname = "data/" + get<TUPLE_BOARD_ID>(data) + "_" + get<TUPLE_ADDRESS>(data);

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

/// @param deviceRawData Requires a raw bit string of a single device
/// @brief Reads in the raw output of the PUF and converts it into comma separated data
string DataParser::commaSeparateData(const string &deviceRawData) {
  string returnString;

  int i = 0;
  while (deviceRawData.compare(i, 1, "\0")) {
	returnString += deviceRawData.at(i);
	if (deviceRawData.compare(i + 1, 1, "\0")) {
	  returnString += ",";
	}
	i++;
  }
  return returnString;
}

/// @param data List that contains samples (ideally only from the same board ID & starting address)
/// @brief values close to 1 have a high tendency to be 1, analog for 0
double *DataParser::getProbabilityOfIndex(const list<tuple<int, string, string, string>> &samplesOfUniqueDevice) {
  int arraySize = 4096;

  auto *returnArray = (double *) (calloc(sizeof(double), arraySize));

  for (tuple<int, string, string, string> sample : samplesOfUniqueDevice) {
	for (int i = 0; (i < arraySize) && (i < get<TUPLE_RAWDATA>(sample).size()); ++i) {
	  if (get<TUPLE_RAWDATA>(sample).at(i) == '1') {
		returnArray[i]++;
	  }
	}
  }

  unsigned long sizeOfData = samplesOfUniqueDevice.size();
  for (int i = 0; i < arraySize; ++i) {
	returnArray[i] /= (double) sizeOfData;
  }

  return returnArray;
}

/// @param samplesOfUniqueDevice List of samples that have the same unique board ID
/// @brief This function processes samples of a unique device and splits those samples by their starting address,
/// resulting in a separate pixel picture for each device & starting address
void DataParser::outputGraph(const list<tuple<int, string, string, string>> &samplesOfUniqueDevice) {
  //The List contains only samples for a single device, but the device is probed at different starting addresses,
  //so we need to output a graph for samples having the same board id && starting address
  list<tuple<int, string, string, string>> localSamplesOfDevice(samplesOfUniqueDevice);

  while (!localSamplesOfDevice.empty()) {
	list<tuple<int, string, string, string>> samplesOfDeviceWithEqualAddress;
	string startingAddress = get<TUPLE_ADDRESS>(localSamplesOfDevice.front());
	auto i = localSamplesOfDevice.begin();
	while (i != localSamplesOfDevice.end()) {
	  if (startingAddress == get<TUPLE_ADDRESS>(i.operator*())) {
		samplesOfDeviceWithEqualAddress.emplace_back(i.operator*());
		localSamplesOfDevice.erase(i++);  // alternatively, i = items.erase(i);
	  } else {
		++i;
	  }
	}
	outputSingleImage(samplesOfDeviceWithEqualAddress);
  }
}

/// @brief Will write individual pixels in the file between 0 to [possiblePixelValues],
/// depending on the corresponding values in the double array, ranging vom 0-1
void DataParser::outputSingleImage(const list<tuple<int, string, string, string>> &samplesOfDeviceWithEqualAddress) {
  int check;
  int possiblePixelValues = 255;
  const string &subfolder = "pictures/";
  const string &boardID = get<TUPLE_BOARD_ID>(samplesOfDeviceWithEqualAddress.front());
  const string &address = get<TUPLE_ADDRESS>(samplesOfDeviceWithEqualAddress.front());

  const string &filename = subfolder + "picture_" + boardID + "_" + address;

  check = mkdir(subfolder.c_str(), 0777);
  //is no dir was created, return
  if (check == -1 && errno != EEXIST) {
	cout << "directory could not be created";
	return;
  }

  double *array;
  array = getProbabilityOfIndex(samplesOfDeviceWithEqualAddress);

  ofstream pictureFile(filename + ".pgm");

  pictureFile << "P2" //Image format
			  << "\n"
			  << "64 64" //Image size; //TODO make this dynamic
			  << "\n"
			  << possiblePixelValues //No. of possible values
			  << "\n";

  for (int i = 0; i < 4096; i++) {
	pictureFile << int(array[i] * possiblePixelValues) << " ";
  }

  free(array);
}

/// @brief Returns all Samples that have a given board ID
list<tuple<int, string, string, string>> DataParser::extractSamplesByBoardID(const string &boardID) {
  list<tuple<int, string, string, string>> returnList;

  for (tuple<int, string, string, string> sample : p_listOfSamples) {
	if (get<TUPLE_BOARD_ID>(sample) == boardID) {
	  returnList.emplace_back(sample);
	}
  }

  return returnList;
}

///@brief Returns a set of (unique) board IDs that are located in the list of samples
set<string> DataParser::extractAllBoardIDs() {
  set<string> returnSet;
  for (tuple<int, string, string, string> sample : p_listOfSamples) {
	returnSet.insert(get<TUPLE_BOARD_ID>(sample));
  }

  return returnSet;
}
bool DataParser::testRawData() {
  int sampleID;
  string boardID;
  string rawData;

  //////Case 1
  sampleID = 3207;
  boardID = "0x30314710303537322F80380";
  rawData =
	  "0001101101101000000100110010101100011100110110000111101101101000000110110110100000011001011010010011101101101000000110100110100000010011010001100101101100000000000100110100010000011110001110110000011100100010000000101111101000000011111100111101101101000011000110010100000000111011011010001001100001101000001110110110100000011010011010000001001101000110010110110000000000010011010001000001111000111011000000001111101000000011111100100111101101101000000110110110100000001010010000110001101001100001001000011110000000111011011010000001101101101000000110100010101100011100110110000111101101101000000110110110100011011001011010000011101101101000000110100110100000010011010001100101101100000000000100110100010000111100001110110000011100100010000000101111101000000011111100111101101101000011000110010100000000111011011010001001100001101000001110110110100000011010011010000001001101000110010110110000000000010011010001000011110000111011000000001111101000000011111100100111101101101000000110110110100000001010010000111101101001100000000000001110000000000000101111110011101101101000000110110110100000010000001010110000001111010000001110110110100000011011011010000001000100101011001000011101000100010101010010110101101101101000000000111111010000000000000000110000000000101011000110111101000100010010010010110101101101101000000100010100101001000011111101000000000000000011010100110110000000111011011010000001101101101000000100000010101100010001110100010000111001001011000110110110100000001110010010101010001011111011000000110010001110011010000011000001001101000110100110110000000000010011010001000101101100000000101110110110000000000010111000001011101101101000000000010011101110111011011000001011101101101000000000000010101111111001110100010111101101101000000000000010001010000011111110000100100000100000111110110111101100011000010001100001010000110111101111010100011010000000101111000111000001000111000000000010011100000001010000000000100000000000000000000010000010000011110111100001101101000011100000001011010110000100101100000000000010101111011110000110000000000000001000111111101101100000000000000010001110111011011000000111101101101000000110110110100000011011011010000000001111110000010000000000001101000000001010110011110011010000011110110110100000011011011010001001101001101000011110110110100000011011011010000100001011110000000000010000001010011010011000000001110101001011000110110110100000011101010010101010001011111011000000110010001110011010000011000001001101000110010110110000000000010011010001001011101101100000000000101110000010111011011010000000000100111011101110110110000010111011011010000000000000101011111110011101000111111111111101111000010011111011111110000110000000011000111000001111111111110111100000001111101100000010010001101111101101101000110100110001101000000010001010110001000111011001011110110110100011011011011011000100001111110000000100000000001001111011011010001101101001100100011110110110100000011011011011010100001111110000000000010000001001111011011010000001101001100101011110110110100000000000001000101000001111111000010010000010000000000001001000110000011111100000011110110110100000011011011010000001101101101000000000111111000001000000000000110100000000101011110111111101000100000000001000110001100001000110000100000011011110111101010001101000000010111101000000001011111100001000000000000000000000100000100000111101111000011011010000111000000010110101100001001011000000000000101011110111100001100000011110110110100001011011011010101111101101100000111110110110100011011011011011000000001111110000010100000000001100000000001010110011110111010001111110110110100011011011011011000100001111110100000000000111001011111011011010001101101001100100111110110110100000011011011010001001101101101000000000111111000001000000010100110000000000101011001011001101000111111011011010001001001111111000001001000011000000000000001010110010011111010001111110110110100000011011011010000001101101101011000000111111000011111000011100110000000000101011000001101101000011111011011010000001101101101000100110110110100000000011111101001000000001100011000000000010101100011001110100011111101101101000";

  for (tuple<int, string, string, string> sample : p_listOfSamples) {
	if (get<TUPLE_SAMPLE_NUMBER>(sample) == sampleID &&
		get<TUPLE_BOARD_ID>(sample) == boardID) {
	  cout << get<TUPLE_RAWDATA>(sample) << endl << endl << endl;
	  cout << rawData << endl;
	  if (get<TUPLE_RAWDATA>(sample) == rawData) {
		return true;
	  }

	}
  }
  return false;
}








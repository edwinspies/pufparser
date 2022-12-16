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
#include <valarray>

using namespace std;

DataParser::DataParser(const string &fileName, bool fileFormat) {
  altFileFormat = fileFormat;
  currentSample = 0;
  getDataFromCSV(fileName);
}

DataParser::~DataParser() {
  p_listOfSamples.clear();
}

/// @brief Reads in a CSV that has been saved in the format that is available at the Grenoble University; check enum STD_CSV_STRUCTURE
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
	returnTuple = getTupleFromAlternateCSVToken(result);
  } else {
	returnTuple = getTupleFromStdCSVToken(result);
  }

  //Transform the board id and address string toUpper starting from the third char the leave the leading '0x' intact
  if (get<TUPLE_BOARD_ID>(returnTuple).size() > 2 && get<TUPLE_ADDRESS>(returnTuple).size() > 2) {
	std::transform(get<TUPLE_BOARD_ID>(returnTuple).begin() + 2,
				   get<TUPLE_BOARD_ID>(returnTuple).end(),
				   get<TUPLE_BOARD_ID>(returnTuple).begin() + 2,
				   ::toupper);
	std::transform(get<TUPLE_ADDRESS>(returnTuple).begin() + 2,
				   get<TUPLE_ADDRESS>(returnTuple).end(),
				   get<TUPLE_ADDRESS>(returnTuple).begin() + 2,
				   ::toupper);
  }

  return returnTuple;
}

tuple<int, string, string, string> DataParser::getTupleFromStdCSVToken(const vector<string> &result) {
  tuple<int, string, string, string> returnTuple;
  int sampleNumber;
  size_t pos{};

  //if an input line (for example the last) is empty or shorter than expected,
  // accessing the index at 9 would result in an error
  if (result.size() >= STD_CSV_STRUCTURE::STD_CSV_LAST) {
	try {
	  sampleNumber = stoi(result.at(SAMPLE_NUMBER_STD), &pos);
	}
	catch (std::invalid_argument const &ex) {
	  cout << "std::invalid_argument::what(): " << ex.what() << '\n';
	  sampleNumber = -1;
	}

	get<TUPLE_SAMPLE_NUMBER>(returnTuple) = sampleNumber;
	get<TUPLE_BOARD_ID>(returnTuple) = result.at(BOARD_ID_STD);
	get<TUPLE_ADDRESS>(returnTuple) = result.at(ADDRESS_STD);
	get<TUPLE_RAW_DATA>(returnTuple) = result.at(RAW_DATA_STD);
  }

  return returnTuple;
}

tuple<int, string, string, string> DataParser::getTupleFromAlternateCSVToken(const vector<string> &unsanitizedResult) {
  tuple<int, string, string, string> returnTuple;
  vector<string> result;
  int sampleNumber = currentSample++;
  string rawData;

  for (const auto &i : unsanitizedResult) {
	if (i != "," && !i.empty()) {
	  result.emplace_back(i);
	}
  }

  //if an input line (for example the last) is empty or shorter than expected,
  // accessing the index at 9 would result in an error
  if (result.size() >= ALTERNATE_CSV_LAST) {
	//since there is no sample number in the csv, the dataparser object has its own counter to count samples
	get<TUPLE_SAMPLE_NUMBER>(returnTuple) = sampleNumber;
	get<TUPLE_BOARD_ID>(returnTuple) = result.at(BOARD_ID_ALTERNATE);
	get<TUPLE_ADDRESS>(returnTuple) = result.at(ADDRESS_ALTERNATE);
	rawData = convertRawDataAlternate(result.at(RAW_DATA_ALTERNATE));
	get<TUPLE_RAW_DATA>(returnTuple) = rawData;
  }

  return returnTuple;
}

/// @brief Converts the comma separated 8bit decimal numbers into single bit strings
/// @details Needed because raw data bytes in this CSV format are also comma separated
/// (without this function every raw data byte would be split into its own token)
string DataParser::convertRawDataAlternate(const string &data) {
  string returnString, token, binary;
  string localData = data;
  const string delimiter = ",";
  vector<string> tokens;
  size_t pos;

  while ((pos = localData.find(delimiter)) != string::npos) {
	token = localData.substr(0, pos);
	tokens.emplace_back(token);
	localData.erase(0, pos + delimiter.length());
  }

  for (const string &s : tokens) {
	binary = bitset<8>(stoi(s)).to_string(); //to binary
	returnString += binary;
  }

  return returnString;
}

/// @brief Takes the list of all samples and proceeds to output them into the Python ND function format
/// @details Depending on the no. of samples, this will take a few minutes since many folders and file are created
void DataParser::processAndOutputDataToNDFormat() {
  cout << "This will take a while depending on the no. of samples" << endl;

  int check;
  check = mkdir("data", 0777);
  //if no dir was created, return
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
  int check;
  ofstream outFile;
  const string &dirname = "data/" + get<TUPLE_BOARD_ID>(data) + "_" + get<TUPLE_ADDRESS>(data);

  check = mkdir(dirname.c_str(), 0777);
  //is no dir was created, return
  if (check == -1 && errno != EEXIST) {
	cout << "directory could not be created";
	return;
  }

  //instead of writing a file per sample id, call the file sample and
  const string &filename = "sample";
  outFile.open(dirname + "/" + filename, std::ios_base::app);

  //process data bits into comma separated bits
  const string &deviceData = get<TUPLE_RAW_DATA>(data);
  string commaSeparatedDeviceData = commaSeparateData(deviceData);

  // Write to the file
  outFile << commaSeparatedDeviceData << endl;

  // Close the file
  outFile.close();
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
  //remove last ','
  returnString.pop_back();
  return returnString;
}

/// @param data List that contains samples (ideally only from the same board ID & starting address)
/// @brief values close to 1 have a high tendency to be 1, analog for 0
double *DataParser::getProbabilityOfIndex(const list<tuple<int, string, string, string>> &samplesOfUniqueDevice) {
  const int arraySize = 4096;

  auto *returnArray = (double *) (calloc(sizeof(double), arraySize));
  if (returnArray == nullptr) {
	cout << "Calloc failed; aborting...";
	exit(-1);
  }

  for (tuple<int, string, string, string> sample : samplesOfUniqueDevice) {
	for (int i = 0; (i < arraySize) && (i < get<TUPLE_RAW_DATA>(sample).size()); ++i) {
	  if (get<TUPLE_RAW_DATA>(sample).at(i) == '1') returnArray[i]++;
	}
  }

  unsigned long sizeOfData = samplesOfUniqueDevice.size();
  for (int i = 0; i < arraySize; ++i) {
	returnArray[i] /= (double) sizeOfData;
  }

  return returnArray;
}

void DataParser::prepareBinEntrop() {
  set<string> allBoardIDs = extractAllBoardIDs();
  list<tuple<int, string, string, string>> samples;
  list<tuple<int, string, string, string>> samples2;

  list<list<tuple<int, string, string, string>>> groupedSamples;
  list<list<tuple<int, string, string, string>>> groupedSamples2;

  samples = extractSamplesByBoardID("0x30314710303537322F80380");
  samples2 = extractSamplesByBoardID("0x3430471836733632294300");

  groupedSamples = groupSamplesByAddress(samples);
  groupedSamples2 = groupSamplesByAddress(samples2);

  for (const auto &boardID : allBoardIDs) {
	samples = extractSamplesByBoardID(boardID);

	groupedSamples = groupSamplesByAddress(samples);
  }

  calcBinaryEntropy(groupedSamples.front(), groupedSamples2.front());
  calcBinaryEntropy(groupedSamples.back(), groupedSamples2.back());
  calcBinaryEntropy(groupedSamples.front(), groupedSamples2.back());
}

void DataParser::calcBinaryEntropy(const list<tuple<int, string, string, string>> &firstBoard,
								   const list<tuple<int, string, string, string>> &secondBoard) {
  set<string> allBoardIDs = extractAllBoardIDs();
  list<tuple<int, string, string, string>> samples;
  double *arr, *arr2;
  list<list<tuple<int, string, string, string>>> groupedSamples;
  const string &subfolder = "entropy/";
  const string &filename =
	  subfolder + "entropy_" + get<TUPLE_BOARD_ID>(firstBoard.front()) + "_" + get<TUPLE_BOARD_ID>(secondBoard.front())
		  + "_" + get<TUPLE_ADDRESS>(firstBoard.front()) + "_" + get<TUPLE_ADDRESS>(secondBoard.front());

  int check = mkdir(subfolder.c_str(), 0777);
  //is no dir was created, return
  if (check == -1 && errno != EEXIST) {
	cout << "directory could not be created";
	return;
  }

  arr = getProbabilityOfIndex(firstBoard);
  for (int i = 0; i < 4096; ++i) {
	if (arr[i] != 0.0 && arr[i] != 1.0) {
	  //binary entropy
	  arr[i] = (-arr[i]) * std::log2(arr[i]) - (1 - arr[i]) * std::log2(1 - arr[i]);
	}
  }

  vector<double> first;
  for (int i = 0; i < 4096; i++) {
	first.emplace_back(arr[i]);
  }
  std::sort(first.begin(), first.end());

  arr2 = getProbabilityOfIndex(secondBoard);
  for (int i = 0; i < 4096; ++i) {
	if (arr2[i] != 0.0 && arr2[i] != 1.0) {
	  //binary entropy
	  arr2[i] = (-arr2[i]) * std::log2(arr2[i]) - (1 - arr2[i]) * std::log2(1 - arr2[i]);
	}
  }

  vector<double> second;
  for (int i = 0; i < 4096; i++) {
	second.emplace_back(arr2[i]);
  }
  std::sort(second.begin(), second.end());

  ofstream entropyFile(filename);

  for (int i = 0; i < 4096; i++) {
	entropyFile << "bit " << i << " like, unlike, difference\t"
				<< first[i]
				<< ((first[i] == 0.0 || first[i] == 1.0) ? "\t\t\t" : "\t\t")
				<< second[i]
				<< ((second[i] == 0.0 || second[i] == 1.0) ? "\t\t\t" : "\t\t")
				<< std::abs(second[i] - first[i]) << std::endl;
  }

  entropyFile.close();

  free(arr);
  free(arr2);
}

list<list<tuple<int, string, string, string>>> DataParser::groupSamplesByAddress(const list<tuple<int,
																								  string,
																								  string,
																								  string>> &samplesOfUniqueDevice) {
  list<tuple<int, string, string, string>> localSamplesOfDevice(samplesOfUniqueDevice);
  list<list<tuple<int, string, string, string>>> returnList;

  //loop until the last sample is taken out of the list and into samplesOfDeviceWithEqualAddress
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
	returnList.emplace_back(samplesOfDeviceWithEqualAddress);
  }

  return returnList;
}

/// @param samplesOfUniqueDevice List of samples that have the same unique board ID
/// @brief This function processes samples of a unique device and splits those samples by their starting address,
/// resulting in a separate pixel picture for each device & starting address
void DataParser::outputGraph(const list<tuple<int, string, string, string>> &samplesOfUniqueDevice) {
  //The List contains only samples for a single device, but the device is probed at different starting addresses,
  //so we need to output a graph for samples having the same board id && starting address
  list<tuple<int, string, string, string>> localSamplesOfDevice(samplesOfUniqueDevice);

  //loop until the last sample is taken out of the list and into samplesOfDeviceWithEqualAddress
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
  double *array;
  const int possiblePixelValues = 63;
  const string &subfolder = "pictures/";
  const string noOfSamples = to_string(samplesOfDeviceWithEqualAddress.size());
  const string &boardID = get<TUPLE_BOARD_ID>(samplesOfDeviceWithEqualAddress.front());
  const string &address = get<TUPLE_ADDRESS>(samplesOfDeviceWithEqualAddress.front());

  const string &filename = subfolder + "picture_" + boardID + "_" + address + "_" + noOfSamples;

  check = mkdir(subfolder.c_str(), 0777);
  //is no dir was created, return
  if (check == -1 && errno != EEXIST) {
	cout << "directory could not be created";
	return;
  }

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

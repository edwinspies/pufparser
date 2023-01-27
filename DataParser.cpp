#include "DataParser.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
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
	throw std::invalid_argument("File not found " + fileName);
  }

  if (altFileFormat) p_listOfSamples.pop_front();
  p_listOfSamples.pop_back();
}

/// @brief Reads a single line of the csv and splits it into a tuple
bitBlock DataParser::getNextLineAndSplitIntoTokens(istream &str) {
  vector<string> result;
  bitBlock returnBlock;

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
	returnBlock = getTupleFromAlternateCSVToken(result);
  } else {
	returnBlock = getTupleFromStdCSVToken(result);
  }

  //Transform the board id and address string toUpper starting from the third char the leave the leading '0x' intact
  if (returnBlock.boardID.size() > 2 && returnBlock.address.size() > 2) {
	std::transform(returnBlock.boardID.begin() + 2,
				   returnBlock.boardID.end(),
				   returnBlock.boardID.begin() + 2,
				   ::toupper);
	std::transform(returnBlock.address.begin() + 2,
				   returnBlock.address.end(),
				   returnBlock.address.begin() + 2,
				   ::toupper);
  }

  return returnBlock;
}

bitBlock DataParser::getTupleFromStdCSVToken(const vector<string> &result) {
  bitBlock returnBlock;
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

	returnBlock.sampleNumber = sampleNumber;
	returnBlock.boardID = result.at(BOARD_ID_STD);
	returnBlock.address = result.at(ADDRESS_STD);
	returnBlock.rawData = result.at(RAW_DATA_STD);
  }

  return returnBlock;
}

bitBlock DataParser::getTupleFromAlternateCSVToken(const vector<string> &unsanitizedResult) {
  bitBlock returnBlock;
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
	returnBlock.sampleNumber = sampleNumber;
	returnBlock.boardID = result.at(BOARD_ID_ALTERNATE);
	returnBlock.address = result.at(ADDRESS_ALTERNATE);
	rawData = convertRawDataAlternate(result.at(RAW_DATA_ALTERNATE));
	returnBlock.rawData = rawData;
  }

  return returnBlock;
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

  createFolder("data");

  for (const bitBlock &singleSample : p_listOfSamples) {
	writeDeviceDataIntoFile(singleSample);
  }
}

/// @param data Requires a single sample Tuple that will be written into a file
/// @brief receives a single device sample and processes this data to output a file for the Python ND Function
void DataParser::writeDeviceDataIntoFile(const bitBlock &data) {
  ofstream outFile;
  const string &dirname = "data/" + data.boardID + "_" + data.address;

  createFolder(dirname);

  //instead of writing a file per sample id, call the file sample and
  const string &filename = "sample";
  outFile.open(dirname + "/" + filename, std::ios_base::app);

  //process data bits into comma separated bits
  const string &deviceData = data.rawData;
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
void DataParser::getProbabilityOfIndex(double *array, const int arraySize, const list<bitBlock> &samplesOfUniqueDevice) {
  for (bitBlock sample : samplesOfUniqueDevice) {
	for (int i = 0; (i < arraySize) && (i < sample.rawData.size()); ++i) {
	  if (sample.rawData.at(i) == '1') array[i]++;
	}
  }

  unsigned long sizeOfData = samplesOfUniqueDevice.size();
  for (int i = 0; i < arraySize; ++i) {
	array[i] /= (double) sizeOfData;
  }

}

void DataParser::createFolder(const string &folderName) {
  const int check = mkdir(folderName.c_str(), 0777);
  if (check == -1 && errno != EEXIST) {
	throw std::invalid_argument("folder could not be created with error code: " + to_string(errno));
  }
}

void DataParser::prepareBinaryEntropyOutput() {
  set<string> allBoardIDs = extractAllBoardIDs();
  list<bitBlock> samples;
  list<list<bitBlock>> groupedSamples;

  for (const auto &boardID : allBoardIDs) {
	//if (i++ > 10) break;
	samples = extractSamplesByBoardID(boardID);

	groupedSamples = groupSamplesByAddress(samples);

	for (const auto &listSameAddressSamples : groupedSamples) {
	  calcBinaryEntropy(listSameAddressSamples);
	}
  }
}

void DataParser::calcBinaryEntropy(const list<bitBlock> &firstBoard) {
  set<string> allBoardIDs = extractAllBoardIDs();
  const string &generalFolder = "entropy/";
  //const string &addressSubfolder = firstBoard.front().address + "/";
  const string &filename =
	  generalFolder +
	  //addressSubfolder +
	  firstBoard.front().boardID + "_" +
	  firstBoard.front().address + "_" +
	  to_string(firstBoard.size());

  createFolder(generalFolder);
  //createFolder(generalFolder + addressSubfolder);

  auto *arr = (double *) (calloc(sizeof(double), arraySize));
  if (arr == nullptr) {
	throw std::invalid_argument("Calloc failed; aborting...");
  }

  getProbabilityOfIndex(arr, arraySize, firstBoard);
  for (int i = 0; i < arraySize; ++i) {
	if (arr[i] != 0.0 && arr[i] != 1.0) {
	  //TODO check what happens with 1
	  //binary entropy
	  arr[i] = (-arr[i]) * std::log2(arr[i]) - (1 - arr[i]) * std::log2(1 - arr[i]);
	} else {
	  arr[i] = 0;
	}
  }

  vector<double> first;
  for (int i = 0; i < arraySize; i++) {
	first.emplace_back(arr[i]);
  }
  free(arr);

  ofstream entropyFile(filename);

  for (int i = 0; i < arraySize; i++) {
	entropyFile << first[i] << std::endl;
  }

  entropyFile.close();
}

list<list<bitBlock>> DataParser::groupSamplesByAddress(const list<bitBlock> &samplesOfUniqueDevice) {
  list<bitBlock> localSamplesOfDevice(samplesOfUniqueDevice);
  list<list<bitBlock>> returnList;

  //loop until the last sample is taken out of the list and into samplesOfDeviceWithEqualAddress
  while (!localSamplesOfDevice.empty()) {
	list<bitBlock> samplesOfDeviceWithEqualAddress;
	string startingAddress = localSamplesOfDevice.front().address;
	auto i = localSamplesOfDevice.begin();
	while (i != localSamplesOfDevice.end()) {
	  if (startingAddress == i.operator*().address) {
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

void DataParser::outputProbability(const list<bitBlock> &samplesOfUniqueDevice) {
  list<bitBlock> localSamplesOfDevice(samplesOfUniqueDevice);

  //loop until the last sample is taken out of the list and into samplesOfDeviceWithEqualAddress
  while (!localSamplesOfDevice.empty()) {
	list<bitBlock> samplesOfDeviceWithEqualAddress;
	string startingAddress = localSamplesOfDevice.front().address;
	auto i = localSamplesOfDevice.begin();
	while (i != localSamplesOfDevice.end()) {
	  if (startingAddress == i.operator*().address) {
		samplesOfDeviceWithEqualAddress.emplace_back(i.operator*());
		localSamplesOfDevice.erase(i++);  // alternatively, i = items.erase(i);
	  } else {
		++i;
	  }
	}
	outputSingleProbability(samplesOfDeviceWithEqualAddress);
  }
}

void DataParser::outputSingleProbability(const list<bitBlock> &samplesOfDeviceWithEqualAddress) {
  const string subfolder = "probability/";
  const string noOfSamples = to_string(samplesOfDeviceWithEqualAddress.size());
  const string boardID = samplesOfDeviceWithEqualAddress.front().boardID;
  const string address = samplesOfDeviceWithEqualAddress.front().address;
  const string &filename = subfolder + "prob_" + boardID + "_" + address + "_" + noOfSamples;

  createFolder(subfolder);

  auto *array = callocDoubleArray(arraySize);

  getProbabilityOfIndex(array, arraySize, samplesOfDeviceWithEqualAddress);
  ofstream probabilityFile(filename);

  for (int i = 0; i < arraySize; i++) {
	probabilityFile << array[i] << std::endl;
  }
}

/// @param samplesOfUniqueDevice List of samples that have the same unique board ID
/// @brief This function processes samples of a unique device and splits those samples by their starting address,
/// resulting in a separate pixel picture for each device & starting address
void DataParser::outputGraph(const list<bitBlock> &samplesOfUniqueDevice, const int markBit) {
  //The List contains only samples for a single device, but the device is probed at different starting addresses,
  //so we need to output a graph for samples having the same board id && starting address
  list<bitBlock> localSamplesOfDevice(samplesOfUniqueDevice);

  //loop until the last sample is taken out of the list and into samplesOfDeviceWithEqualAddress
  while (!localSamplesOfDevice.empty()) {
	list<bitBlock> samplesOfDeviceWithEqualAddress;
	string startingAddress = localSamplesOfDevice.front().address;
	auto i = localSamplesOfDevice.begin();
	while (i != localSamplesOfDevice.end()) {
	  if (startingAddress == i.operator*().address) {
		samplesOfDeviceWithEqualAddress.emplace_back(i.operator*());
		localSamplesOfDevice.erase(i++);  // alternatively, i = items.erase(i);
	  } else {
		++i;
	  }
	}
	outputSingleImage(samplesOfDeviceWithEqualAddress, markBit);
  }
}

/// @brief Will write individual pixels in the file between 0 to [possiblePixelValues],
/// depending on the corresponding values in the double array, ranging vom 0-1
void DataParser::outputSingleImage(const list<bitBlock> &samplesOfDeviceWithEqualAddress, const int markBit) {
  const int possiblePixelValues = 255;
  const string subfolder = "pictures/";
  const string noOfSamples = to_string(samplesOfDeviceWithEqualAddress.size());
  const string boardID = samplesOfDeviceWithEqualAddress.front().boardID;
  const string address = samplesOfDeviceWithEqualAddress.front().address;

  const string &filename = subfolder + "picture_" + boardID + "_" + address + "_" + noOfSamples;

  createFolder(subfolder);

  auto *array = callocDoubleArray(arraySize);

  getProbabilityOfIndex(array, arraySize, samplesOfDeviceWithEqualAddress);
  ofstream pictureFile(filename + ".pgm");

  pictureFile << "P3" //Image format
			  << "\n"
			  << "64 64" //Image dimensions; //TODO make this dynamic
			  << "\n"
			  << possiblePixelValues //No. of possible values
			  << "\n";

  for (int i = 0; i < arraySize; i++) {
	if (markBit == i-1 || markBit == i+1) {
	  pictureFile
	  << 255 << " "
	  << "0" << " "
	  << "0" << " ";
	} else {
	  pictureFile
	  <<  int(array[i] * possiblePixelValues) << " "
	  << int(array[i] * possiblePixelValues) << " "
	  << int(array[i] * possiblePixelValues) << " ";
	}

  }

  free(array);
}

double *DataParser::callocDoubleArray(int size) {
  auto *array = (double *) (calloc(sizeof(double), size));
  if (array == nullptr) {
	throw std::invalid_argument("Calloc failed; aborting...");
  }

  return array;
}

/// @brief Returns all Samples that have a given board ID
list<bitBlock> DataParser::extractSamplesByBoardID(const string &boardID) {
  list<bitBlock> returnList;

  for (const bitBlock &sample : p_listOfSamples) {
	if (sample.boardID == boardID) {
	  returnList.emplace_back(sample);
	}
  }

  return returnList;
}

///@brief Returns a set of (unique) board IDs that are located in the list of samples
set<string> DataParser::extractAllBoardIDs() {
  set<string> returnSet;
  for (const bitBlock &sample : p_listOfSamples) {
	returnSet.insert(sample.boardID);
  }

  return returnSet;
}

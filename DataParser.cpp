#include "DataParser.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <bitset>
#include <valarray>
#include <array>
#include <cstdlib>
#include <cstdio>

using namespace std;

DataParser::DataParser(const string &fileName, bool fileFormat) {
  altFileFormat = fileFormat;
  currentSample = 0;
  getDataFromCSV(fileName);
}

DataParser::~DataParser() {
  p_vectorOfSamples.clear();
}

/// @brief Reads in a CSV that has been saved in the format that is available at the Grenoble University; check enum STD_CSV_STRUCTURE
void DataParser::getDataFromCSV(const string &fileName) {
  filebuf fb;

  //if this function is called twice, clear all items before to remove old data
  p_vectorOfSamples.clear();

  if (fb.open(fileName, ios::in)) {
	istream is(&fb);
	while (is) {
	  p_vectorOfSamples.emplace_back(getNextLineAndSplitIntoTokens(is));
	}
	fb.close();
  } else {
	throw std::invalid_argument("File not found " + fileName);
  }

  p_vectorOfSamples.erase(p_vectorOfSamples.begin());
  p_vectorOfSamples.erase(p_vectorOfSamples.end());
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

double binary_entropy(double value) {
  if (value != 0.0 && value != 1.0) {
	//TODO check what happens with 1
	//binary entropy
	value = (-value) * std::log2(value) - (1 - value) * std::log2(1 - value);
  } else {
	value = 0;
  }

  return value;
}

bool sortbitblock(const bitBlock &a, const bitBlock &b) { return a.sampleNumber < b.sampleNumber; }

//bool sortbitranking()

int DataParser::tryTo3DData() {
  //10 boards á 4096 bits and 7 samples
  //double bitMatrix[10][4096*64][7]; //this will lead to a segmentation fault
  const int maxSamples = MAX_SAMPLES;
  const int maxBoards = MAX_BOARDS;
  //bitMatrix[maxBoards][4096*64][maxSamples]
  std::vector<std::array<std::array<bool, maxSamples>, 4096 * 64 >> bitMatrix(maxBoards);
  //this vector will help to find out the board id because the 3D matrix has no additional space to store this info
  std::vector<std::string> helperData;


  //temporary variable to convert a char to an int
  char c[2] = {0};

  set<string> s = extractAllBoardIDs();
  vector<bitBlock> b;

  int currentBoardCount = 0;

  //this will iterate through all boards and get the individual samples from a specific board
  //continues to place the bits of this board into the 3D matrix
  for (const string &board : s) {
	b = extractSamplesByBoardID(board);
	std::sort(b.begin(), b.end(), sortbitblock);

	//check that a board has more samples than "maxSamples"; skips boards with less
	if (b.size() < 64 * maxSamples) continue;
	helperData.emplace_back(b.front().boardID);
	for (int currentSampleCount = 0; currentSampleCount < 7; currentSampleCount++) {
	  for (int i = 0; i < 64; i++) { //each board is represented by 64 memory regions
		for (int j = 0; j < b[i].rawData.size(); j++) { //each memory region consists of 4096 bits
		  if (b[i].rawData.size() != 4096) throw std::invalid_argument("size should be 4096");
		  c[0] = b[i].rawData.at(j);
		  bitMatrix[currentBoardCount][i * 4096 + j][currentSampleCount] = atoi(c);
		}
	  }
	  //remove already used samples (this will remove 64 bitblocks which equals to a readout of the whole board)
	  // -> the next 64 items in the vector will be those from e.g. another day
	  b.erase(b.begin(), b.begin() + 64);
	}
	if (currentBoardCount++ >= maxBoards) break;
  }

  calcMetrics32Increments(bitMatrix, helperData);


//  int count = 0;
//  for (int d = 5; d < 50; d) {
//	for (int i = 0; i < (4096*64)-d-1; i++) {
//	  count += bitMatrix[0][i][0] ^ bitMatrix[0][i+d][0];
//	}
//	cout << "This is the result of autocorrelation("<< d << "):" << count << std::endl;
//
//	double x5 = (2*(count-((4096*64 - (double) d)/2)))/(std::sqrt(4096*64-d));
//	cout << "x5 of " << d << " is: " << x5 << std::endl;
//	d=(d+3)*2;
//  }



////////////////// Uniqueness
//	double totalCount = 0;
//	for (int u = 0; u < maxBoards-1; u++) {
//	  for (int v = u+1; v < maxBoards; v++) {
//		double hdCount = 0;
//		for(int i = 0; i < 4096*64; i++) {
//		  hdCount += bitMatrix[u][i][0] ^ bitMatrix[v][i][0];
//		}
//		hdCount = hdCount / (4096*64);
//		hdCount*=100;
//		totalCount += hdCount;
//	  }
//	}
//
//	double result;
//	result = ((double) 2*totalCount) / (maxBoards*(maxBoards-1));
//
//	cout << "Uniqueness: " << result << std::endl;





  //double reliability[maxBoards][4096 * 64] = {};
  std::vector<std::array<double, 4096 * 64 >> reliability(maxBoards);
  double uniformity[maxBoards] = {};

  //FOR EVERY BOARD: do Bit evaluation
  for (int board = 0; board < maxBoards; board++) {
	//TODO maybe make the uniformity calculation analog to the probability difference calculation where
	//RELIABILITY / UNIFORMITY
	for (int i = 0; i < bitMatrix[board].size(); i++) {
	  uniformity[board] += bitMatrix[board][i][0];

	  for (int sample = 0; sample < maxSamples; sample++) {
		reliability[board][i] += bitMatrix[board][i][sample];
	  }
	  reliability[board][i] /= maxSamples;
	}
	uniformity[board] /= (double) bitMatrix[board].size();
  }



  //double bitaliasing[maxBoards][4096 * 64] = {0};
  std::vector<std::array<double, 4096 * 64 >> bitaliasing(maxBoards);
  //do a pairwise calculation of bit aliasing
  for (int board1 = 0; board1 < maxBoards; ++board1) {
	for (int board2 = 0; board2 < maxBoards; ++board2) {
	  if (board1 == board2) continue;

	  //do the calculation for every cell of the PUF
	  for (int i = 0; i < bitMatrix[i].size(); ++i) {
		bitaliasing[board1][i] += std::abs(bitMatrix[board1][i][0] - bitMatrix[board2][i][0]);
	  }
	}

	//now take the average of this pairwise comparison
	for (int i = 0; i < bitMatrix[i].size(); ++i) {
	  bitaliasing[board1][i] /= (maxBoards-1);
	}
  }



//  int countSmaller03 = 0;
//  int countBetween0305 = 0;
//  int countBetween0507 = 0;
//  int countGreater07 = 0;
//  for(int i = 0; i < 4096*64; i++) {
//	if (bitaliasing[0][i] < 0.3) {
//	  countSmaller03++;
//	} else if (bitaliasing[0][i] > 0.3 && bitaliasing[0][i] < 0.5) {
//	  countBetween0305++;
//	} else if (bitaliasing[0][i] >= 0.5 && bitaliasing[0][i] < 0.7) {
//	  countBetween0507++;
//	} else if (bitaliasing[0][i] >= 0.7) {
//	  countGreater07++;
//	}
//  }
//
//  cout << "Values smaller than 0.3: " << countSmaller03 << std::endl;
//  cout << "Values between 0.3 and 0.5: " << countBetween0305 << std::endl;
//  cout << "Values between 0.5 and 0.7: " << countBetween0507 << std::endl;
//  cout << "Values greater than 0.7: " << countGreater07 << std::endl;


/*  double bitAliasingAverage[4096*64];
  double reliabilityAverage[4096*64];
  double uniformityAverage = 0;
  //average over all boards
  for (int i = 0; i < 4096*64; i++) {
	for(int j = 0; j < maxBoards; j++) {
	  bitAliasingAverage[i] += bitaliasing[j][i];
	  reliabilityAverage[i] += reliability[j][i];
	}
	bitAliasingAverage[i] /= maxBoards;
	reliabilityAverage[i] /= maxBoards;
  }

  for (int i = 0; i < maxBoards; i++) {
	uniformityAverage += uniformity[i];
  }
  uniformityAverage /= maxBoards;

  ofstream bitRankingFileAverage("bitranks/bitranking_for_board_ALL");
  for (int i = 0; i < 4096 * 64; i++) {
	bitRankingFileAverage << (std::min(bitAliasingAverage[i], uniformityAverage) - binary_entropy(reliabilityAverage[i]))
				   << std::endl;
  }
  bitRankingFileAverage.close();*/


  //outputBitRanksAllBoards(maxBoards, bitaliasing, uniformity, reliability, helperData);

  //outputBitRanksGlobalAverage(maxBoards, bitaliasing, uniformity, reliability, helperData);

  return 0;
}

void DataParser::calcMetrics32Increments(const std::vector<std::array<std::array<bool, MAX_SAMPLES>, 4096 * 64 >> &bitMatrix,
										 const std::vector<std::string> &helperData) {
  const int maxBoards = MAX_BOARDS;
  const int maxSamples = MAX_SAMPLES;

  std::vector<std::array<double, 4096 * 32 >> reliability(maxBoards);
  int reliabilityIndexCounter = 0;
  double uniformity[maxBoards] = {};

  //FOR EVERY BOARD: do Bit evaluation
  for (int board = 0; board < maxBoards; board++) {
	//TODO maybe make the uniformity calculation analog to the probability difference calculation where
	//RELIABILITY / UNIFORMITY
	bool useBits = true;
	for (int i = 0; i < bitMatrix[board].size(); i++) {
	  if (i % 32 == 0) useBits = !useBits;
	  if (useBits) {
	  uniformity[board] += bitMatrix[board][i][0];

	  for (int sample = 0; sample < maxSamples; sample++) {
		reliability[board][reliabilityIndexCounter] += bitMatrix[board][i][sample];
	  }
	  reliability[board][reliabilityIndexCounter] /= maxSamples;
	  reliabilityIndexCounter++;
	  }
	}
	reliabilityIndexCounter = 0;
	uniformity[board] /= 4096*32;
  }

  std::vector<std::array<double, 4096 * 32 >> bitaliasing(maxBoards);
  int bitaliasingIndexCounter = 0;
  //do a pairwise calculation of bit aliasing
  for (int board1 = 0; board1 < maxBoards; ++board1) {
	for (int board2 = 0; board2 < maxBoards; ++board2) {
	  if (board1 == board2) continue;

	  bool useBits = true;
	  //do the calculation for every cell of the PUF
	  for (int i = 0; i < bitMatrix[i].size(); ++i) {
		if (i % 32 == 0) useBits = !useBits;
		if (useBits) {
		  bitaliasing[board1][bitaliasingIndexCounter] += std::abs(bitMatrix[board1][i][0] - bitMatrix[board2][i][0]);
		  bitaliasingIndexCounter++;
		}
	  }
	  bitaliasingIndexCounter = 0;
	}

	//now take the average of this pairwise comparison
	for (int i = 0; i < 4096 * 32; ++i) {
	  bitaliasing[board1][i] /= (maxBoards-1);
	}
  }


  createFolder("bitranks");
  //Output the BIT RANKING
  double currentBitrank;
  for (int board = 0; board < maxBoards; board++) {
	ofstream bitRankingFile("bitranks/32_increment_bitranking_for_board_" + helperData[board]);
	for (int i = 0; i < 4096 * 32; i++) {
	  //this can provide an actual rank by sorting; instead write all values directly into the output file
	  //bitRanking.emplace_back(i, std::min(bitaliasing[i], uniformity) - binary_entropy(reliability[i]));
	  currentBitrank = (std::min(bitaliasing[board][i], uniformity[board]) - binary_entropy(reliability[board][i]));
	  bitRankingFile << currentBitrank << std::endl;
	}
	bitRankingFile.close();
  }


}

void DataParser::sortCellsByBitRanks() {
  
}

void DataParser::outputBitRanksAllBoards(int maxBoards, const std::vector<std::array<double, 4096 * 64 >> &bitaliasing,
							 double uniformity[],
							 const std::vector<std::array<double, 4096 * 64 >> &reliability,
							 const std::vector<std::string> &helperData) {
  createFolder("bitranks");
  //Output the BIT RANKING
  double currentBitrank;
  for (int board = 0; board < maxBoards; board++) {
	ofstream bitRankingFile("bitranks/bitranking_for_board_" + helperData[board]);
	for (int i = 0; i < 4096 * 64; i++) {
	  //this can provide an actual rank by sorting; instead write all values directly into the output file
	  //bitRanking.emplace_back(i, std::min(bitaliasing[i], uniformity) - binary_entropy(reliability[i]));
	  currentBitrank = (std::min(bitaliasing[board][i], uniformity[board]) - binary_entropy(reliability[board][i]));
	  bitRankingFile << currentBitrank << std::endl;
	}
	bitRankingFile.close();
  }
}

void DataParser::outputBitRanksGlobalAverage(int maxBoards, const std::vector<std::array<double, 4096 * 64 >> &bitaliasing,
										 double uniformity[],
										 const std::vector<std::array<double, 4096 * 64 >> &reliability,
										 const std::vector<std::string> &helperData) {
  //////////////////////// Takes the end result of the bit ranking metric and takes the average over all boards //////////
  createFolder("bitranks");
  //Output the BIT RANKING
  double bitRankingMatrix[4096*64] = {};
  double currentBitrank;
  for (int board = 0; board < maxBoards; board++) {
	for (int i = 0; i < 4096 * 64; i++) {
	  //this can provide an actual rank by sorting; instead write all values directly into the output file
	  currentBitrank = (std::min(bitaliasing[board][i], uniformity[board]) - binary_entropy(reliability[board][i]));
	  bitRankingMatrix[i] += currentBitrank;
	}
  }


  //takes average over generated bitranking values
  ofstream bitRankingFile("bitranks/bitranking_for_board_endresultAverage");
  for(int i = 0; i < 4096*64; i++) {
	bitRankingMatrix[i] /= maxBoards;
	bitRankingFile << bitRankingMatrix[i] << std::endl;
  }
  bitRankingFile.close();

  ///////////////////////////////////////////////////

  ////////////////////////// Takes the average of the individual metrics before applying the bit ranking metric /////////

  double bitAliasingAverage[4096*64];
  double reliabilityAverage[4096*64];
  double uniformityAverage = 0;
  //average over all boards
  for (int i = 0; i < 4096*64; i++) {
	for(int j = 0; j < maxBoards; j++) {
	  bitAliasingAverage[i] += bitaliasing[j][i];
	  reliabilityAverage[i] += reliability[j][i];
	}
	bitAliasingAverage[i] /= maxBoards;
	reliabilityAverage[i] /= maxBoards;
  }

  for (int i = 0; i < maxBoards; i++) {
	uniformityAverage += uniformity[i];
  }
  uniformityAverage /= maxBoards;

  ofstream bitRankingFileAverage("bitranks/bitranking_for_board_individualMetricsAverage");
  for (int i = 0; i < 4096 * 64; i++) {
	bitRankingFileAverage << (std::min(bitAliasingAverage[i], uniformityAverage) - binary_entropy(reliabilityAverage[i]))
				   << std::endl;
  }
  bitRankingFileAverage.close();


}

void calculateReliability() {

}

/// @brief Takes the vector of all samples and proceeds to output them into the Python ND function format
/// @details Depending on the no. of samples, this will take a few minutes since many folders and file are created
void DataParser::processAndOutputDataToNDFormat() {
  cout << "This will take a while depending on the no. of samples" << endl;

  createFolder("data");

  for (const bitBlock &singleSample : p_vectorOfSamples) {
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

/// @param data Vector that contains samples (ideally only from the same board ID & starting address)
/// @brief values close to 1 have a high tendency to be 1, analog for 0
void DataParser::getProbabilityOfIndex(double *array,
									   const int arraySize,
									   const vector<bitBlock> &samplesOfUniqueDevice) {
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
  vector<bitBlock> samples;
  list<vector<bitBlock>> groupedSamples;

  for (const auto &boardID : allBoardIDs) {
	//if (i++ > 10) break;
	samples = extractSamplesByBoardID(boardID);

	groupedSamples = groupSamplesByAddress(samples);

	for (const auto &vectorSameAddressSamples : groupedSamples) {
	  calcBinaryEntropy(vectorSameAddressSamples);
	}
  }
}

void DataParser::calcBinaryEntropy(const vector<bitBlock> &firstBoard) {
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

list<vector<bitBlock>> DataParser::groupSamplesByAddress(const vector<bitBlock> &samplesOfUniqueDevice) {
  vector<bitBlock> localSamplesOfDevice(samplesOfUniqueDevice);
  list<vector<bitBlock>> returnList;

  //loop until the last sample is taken out of the list and into samplesOfDeviceWithEqualAddress
  while (!localSamplesOfDevice.empty()) {
	vector<bitBlock> samplesOfDeviceWithEqualAddress;
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

void DataParser::outputProbability(const vector<bitBlock> &samplesOfUniqueDevice) {
  vector<bitBlock> localSamplesOfDevice(samplesOfUniqueDevice);

  //loop until the last sample is taken out of the vector and into samplesOfDeviceWithEqualAddress
  while (!localSamplesOfDevice.empty()) {
	vector<bitBlock> samplesOfDeviceWithEqualAddress;
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

void DataParser::outputSingleProbability(const vector<bitBlock> &samplesOfDeviceWithEqualAddress) {
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

/// @param samplesOfUniqueDevice Vector of samples that have the same unique board ID
/// @brief This function processes samples of a unique device and splits those samples by their starting address,
/// resulting in a separate pixel picture for each device & starting address
void DataParser::outputGraph(const vector<bitBlock> &samplesOfUniqueDevice, const int markBit) {
  //The Vector contains only samples for a single device, but the device is probed at different starting addresses,
  //so we need to output a graph for samples having the same board id && starting address
  vector<bitBlock> localSamplesOfDevice(samplesOfUniqueDevice);

  //loop until the last sample is taken out of the vector and into samplesOfDeviceWithEqualAddress
  while (!localSamplesOfDevice.empty()) {
	vector<bitBlock> samplesOfDeviceWithEqualAddress;
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
void DataParser::outputSingleImage(const vector<bitBlock> &samplesOfDeviceWithEqualAddress, const int markBit) {
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
	if (markBit == i - 1 || markBit == i + 1) {
	  pictureFile
		  << 255 << " "
		  << "0" << " "
		  << "0" << " ";
	} else {
	  pictureFile
		  << int(array[i] * possiblePixelValues) << " "
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
vector<bitBlock> DataParser::extractSamplesByBoardID(const string &boardID) {
  vector<bitBlock> returnVector;

  for (const bitBlock &sample : p_vectorOfSamples) {
	if (sample.boardID == boardID) {
	  returnVector.emplace_back(sample);
	}
  }

  return returnVector;
}

///@brief Returns a set of (unique) board IDs that are located in the vector of samples
set<string> DataParser::extractAllBoardIDs() {
  set<string> returnSet;
  for (const bitBlock &sample : p_vectorOfSamples) {
	returnSet.insert(sample.boardID);
  }

  return returnSet;
}

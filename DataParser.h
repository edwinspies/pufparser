#ifndef PUFPARSER_DATAPARSER_H_
#define PUFPARSER_DATAPARSER_H_

#include <string>
#include <vector>
#include <list>
#include <set>
#include <array>

#define MAX_SAMPLES 10
#define MAX_BOARDS 48

using namespace std;

/// @brief Delivers the index of the std csv file; can be modified to read other csv structured puf data
enum STD_CSV_STRUCTURE {
  SAMPLE_NUMBER_STD,
  DATE_STD,
  BOARD_TYPE_STD,
  BOARD_ID_STD,
  UNDEFINED_COLUMN_4,
  ADDRESS_STD,
  NO_OF_BITS_STD,
  UNDEFINED_COLUMN_7,
  UNDEFINED_COLUMN_8,
  RAW_DATA_STD,
  STD_CSV_LAST = RAW_DATA_STD
  //STD_CSV_LAST needs to be updated if the structure changes,
  //this is required to not hardcode any index values in other functions
};

/// @brief Delivers the index of the csv file that ALTERNATE provided; can be modified to read other csv structured puf data
enum ALTERNATE_CSV_STRUCTURE {
  BOARD_TYPE_ALTERNATE,
  BOARD_ID_ALTERNATE,
  POSITION_IN_CHAIN_ALTERNATE,
  ADDRESS_ALTERNATE,
  RAW_DATA_ALTERNATE,
  DATE_ALTERNATE,
  ALTERNATE_CSV_LAST = DATE_ALTERNATE
  //ALTERNATE_CSV_LAST needs to be updated if the structure changes,
  //this is required to not hardcode any index values in other functions
};

struct bitBlock{
  int sampleNumber = 0;
  string boardID;
  string address;
  string rawData;
};

class DataParser {
 public:
  explicit DataParser(const string &fileName, bool altFileFormat);
  ~DataParser();

  void processAndOutputDataToNDFormat();
  vector<bitBlock> extractSamplesByBoardID(const string &boardID);
  static void outputGraph(const vector<bitBlock> &samplesOfUniqueDevice, int markBit);
  set<string> extractAllBoardIDs();
  void prepareBinaryEntropyOutput();
  static void outputProbability(const vector<bitBlock> &samplesOfUniqueDevice);

  int tryTo3DData();

 private:
  void getDataFromCSV(const string &fileName);
  bitBlock getNextLineAndSplitIntoTokens(istream &str);
  static bitBlock getTupleFromStdCSVToken(const vector<string> &result);
  bitBlock getTupleFromAlternateCSVToken(const vector<string> &result);
  static string convertRawDataAlternate(const string &data);
  static void writeDeviceDataIntoFile(const bitBlock &data);
  static string commaSeparateData(const string &deviceData);
  void calcBinaryEntropy(const vector<bitBlock> &firstBoard);
  static void getProbabilityOfIndex(double *array, int arraySize, const vector<bitBlock> &samplesOfUniqueDevice);
  static list<vector<bitBlock>> groupSamplesByAddress(const vector<bitBlock> &samplesOfUniqueDevice);
  void outputBitRanksAllBoards(int maxBoards, const std::vector<std::array<double, 4096 * 64 >> &bitaliasing,
							   double uniformity[],
							   const std::vector<std::array<double, 4096 * 64 >> &reliability, const std::vector<std::string> &helperData);
  void outputBitRanksGlobalAverage(int maxBoards, const std::vector<std::array<double, 4096 * 64 >> &bitaliasing,
											   double uniformity[],
											   const std::vector<std::array<double, 4096 * 64 >> &reliability,
											   const std::vector<std::string> &helperData);
  void outputBitRanks32IncrementGlobalAverage(int maxBoards, const std::vector<std::array<double, 4096 * 32 >> &bitaliasing,
														  double uniformity[],
														  const std::vector<std::array<double, 4096 * 32 >> &reliability,
														  const std::vector<std::string> &helperData);
  void calcMetrics32Increments(const std::vector<std::array<std::array<bool, MAX_SAMPLES>, 4096 * 64 >> &bitMatrix,
							   const std::vector<std::string> &helperData);

  static void outputSingleProbability(const vector<bitBlock> &samplesOfDeviceWithEqualAddress);
  static void outputSingleImage(const vector<bitBlock> &samplesOfDeviceWithEqualAddress, int markBit);
  static void createFolder(const string& folderName);
  static double *callocDoubleArray(int size);
  void sortCellsByBitRanks();

  vector<bitBlock> p_vectorOfSamples;
  bool altFileFormat;
  int currentSample;
  static const int arraySize = 4096;
};

#endif //PUFPARSER_DATAPARSER_H_
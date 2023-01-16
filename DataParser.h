#ifndef PUFPARSER_DATAPARSER_H_
#define PUFPARSER_DATAPARSER_H_

#include <string>
#include <vector>
#include <list>
#include <set>

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
  int sampleNumber;
  string boardID;
  string address;
  string rawData;
};

class DataParser {
 public:
  explicit DataParser(const string &fileName, bool altFileFormat);
  ~DataParser();

  void processAndOutputDataToNDFormat();
  list<bitBlock> extractSamplesByBoardID(const string &boardID);
  static void outputGraph(const list<bitBlock> &samplesOfUniqueDevice, int markBit);
  set<string> extractAllBoardIDs();
  void prepareBinaryEntropyOutput();

 private:
  void getDataFromCSV(const string &fileName);
  bitBlock getNextLineAndSplitIntoTokens(istream &str);
  static bitBlock getTupleFromStdCSVToken(const vector<string> &result);
  bitBlock getTupleFromAlternateCSVToken(const vector<string> &result);
  static string convertRawDataAlternate(const string &data);
  static void writeDeviceDataIntoFile(const bitBlock &data);
  static string commaSeparateData(const string &deviceData);
  void calcBinaryEntropy(const list<bitBlock> &firstBoard);
  static void getProbabilityOfIndex(double *array, int arraySize, const list<bitBlock> &samplesOfUniqueDevice);
  static list<list<bitBlock>> groupSamplesByAddress(const list<bitBlock> &samplesOfUniqueDevice);
  static void outputSingleImage(const list<bitBlock> &samplesOfDeviceWithEqualAddress, int markBit);
  static void createFolder(const string& folderName);

  list<bitBlock> p_listOfSamples;
  bool altFileFormat;
  int currentSample;
  static const int arraySize = 4096;
};

#endif //PUFPARSER_DATAPARSER_H_
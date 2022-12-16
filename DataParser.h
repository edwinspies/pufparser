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

/// @brief Defines the index for the sample tuple which is standardized for all read csv formats
enum TUPLE_STRUCTURE {
  TUPLE_SAMPLE_NUMBER,
  TUPLE_BOARD_ID,
  TUPLE_ADDRESS,
  TUPLE_RAW_DATA
};

class DataParser {
 public:
  explicit DataParser(const string &fileName, bool altFileFormat);
  ~DataParser();

  void processAndOutputDataToNDFormat();
  list<tuple<int, string, string, string>> extractSamplesByBoardID(const string &boardID);
  static void outputGraph(const list<tuple<int, string, string, string>> &samplesOfUniqueDevice);
  set<string> extractAllBoardIDs();
  void prepareBinEntrop();

 private:
  void getDataFromCSV(const string &fileName);
  tuple<int, string, string, string> getNextLineAndSplitIntoTokens(istream &str);
  static tuple<int, string, string, string> getTupleFromStdCSVToken(const vector<string> &result);
  tuple<int, string, string, string> getTupleFromAlternateCSVToken(const vector<string> &result);
  static string convertRawDataAlternate(const string &data);
  static void writeDeviceDataIntoFile(const tuple<int, string, string, string> &data);
  static string commaSeparateData(const string &deviceData);
  void calcBinaryEntropy(const list<tuple<int, string, string, string>> &firstBoard,
						 const list<tuple<int, string, string, string>> &secondBoard);
  static double *getProbabilityOfIndex(const list<tuple<int, string, string, string>> &samplesOfUniqueDevice);
  static list<list<tuple<int, string, string, string>>> groupSamplesByAddress(const list<tuple<int,
																							   string,
																							   string,
																							   string>> &samplesOfUniqueDevice);
  static void outputSingleImage(const list<tuple<int, string, string, string>> &samplesOfDeviceWithEqualAddress);

  list<tuple<int, string, string, string>> p_listOfSamples;
  bool altFileFormat;
  int currentSample;
};

#endif //PUFPARSER_DATAPARSER_H_
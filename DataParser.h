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

/// @brief Delivers the index of the csv file that sergio provided; can be modified to read other csv structured puf data
enum SERGIO_CSV_STRUCTURE {
  BOARD_TYPE_SERGIO,
  BOARD_ID_SERGIO,
  POSITION_IN_CHAIN_SERGIO,
  ADDRESS_SERGIO,
  RAW_DATA_SERGIO,
  DATE_SERGIO,
  SERGIO_CSV_LAST = DATE_SERGIO
  //SERGIO_CSV_LAST needs to be updated if the structure changes,
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

 private:
  void getDataFromCSV(const string &fileName);
  tuple<int, string, string, string> getNextLineAndSplitIntoTokens(istream &str);
  static tuple<int, string, string, string> getTupleFromStdCSVToken(const vector<string> &result);
  tuple<int, string, string, string> getTupleFromSergioCSVToken(const vector<string> &result);
  static void writeDeviceDataIntoFile(const tuple<int, string, string, string> &data);
  static string commaSeparateData(const string &deviceData);
  static double *getProbabilityOfIndex(const list<tuple<int, string, string, string>> &samplesOfUniqueDevice);
  static void outputSingleImage(const list<tuple<int, string, string, string>> &samplesOfDeviceWithEqualAddress);
  static string convertRawDataSergio(const string &data);

  list<tuple<int, string, string, string>> p_listOfSamples;
  bool altFileFormat;
  int currentSample;
};

#endif //PUFPARSER_DATAPARSER_H_
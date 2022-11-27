#ifndef PUFPARSER_DATAPARSER_H_
#define PUFPARSER_DATAPARSER_H_

#include <string>
#include <vector>
#include <list>
#include <set>

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
  RAWDATA,
  STD_CSV_LAST = RAWDATA //Needs to be updated if the structure changes, this is required to not hardcode any index values in other functions
};

/// @brief Delivers the index of the csv file that sergio provided; can be modified to read other csv structured puf data
enum SERGIO_CSV_STRUCTURE {
  BOARD_TYPE_SERGIO,
  BOARD_ID_SERGIO,
  POSITION_IN_CHAIN_SERGIO,
  ADDRESS_SERGIO,
  RAWDATA_SERGIO,
  DATE_SERGIO,
  SERGIO_CSV_LAST = DATE_SERGIO //Needs to be updated if the structure changes, this is required to not hardcode any index values in other functions
};

/// @brief Defines the index for the sample tuple which is standardized for all read csv formats
enum TUPLE_STRUCTURE {
  TUPLE_SAMPLE_NUMBER,
  TUPLE_BOARD_ID,
  TUPLE_ADDRESS,
  TUPLE_RAWDATA
};

class DataParser {
 public:
  explicit DataParser(const string &fileName, bool altFileFormat);
  ~DataParser();

  void processAndOutputDataToNDFormat();

  list<tuple<int, string, string, string>> extractSamplesByBoardID(const string &boardID);
  static void outputGraph(const list<tuple<int, string, string, string>> &samplesOfUniqueDevice);
  set<string> extractAllBoardIDs();

  ///////////////////TEST///////////////////
  bool testRawData();

 private:
  void getDataFromCSV(const string &fileName);
  tuple<int, string, string, string> getNextLineAndSplitIntoTokens(istream &str);
  static tuple<int, string, string, string> getTupleFromStdCSVToken(const vector<string>& result);
  tuple<int, string, string, string> getTupleFromSergioCSVToken(const vector<string>& result);
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
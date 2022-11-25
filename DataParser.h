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
  RAWDATA
};

/// @brief Defines the index for the sample tuple
enum TUPLE_STRUCTURE {
  TUPLE_SAMPLE_NUMBER,
  TUPLE_BOARD_ID,
  TUPLE_RAWDATA
};

class DataParser {
 public:
  explicit DataParser(const string &fileName);
  ~DataParser();

  void processAndOutputDataToNDFormat();

  list<tuple<int, string, string>> extractSamplesByBoardID(const string &boardID);
  static void outputGraph(const list<tuple<int, string, string>> &samplesOfUniqueDevice);
  set<string> extractAllBoardIDs();

 private:
  void getDataFromCSV(const string &fileName);
  static tuple<int, string, string> getNextLineAndSplitIntoTokens(istream &str);
  static void writeDeviceDataIntoFile(const tuple<int, string, string> &data);
  static string commaSeparateData(const string &deviceData);
  static double *getProbabilityOfIndex(const list<tuple<int, string, string>> &samplesOfUniqueDevice);

  list<tuple<int, string, string>> p_listOfSamples;
};

#endif //PUFPARSER_DATAPARSER_H_
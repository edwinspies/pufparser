#ifndef PUFPARSER_DATAPARSER_H_
#define PUFPARSER_DATAPARSER_H_

#include <string>
#include <vector>
#include <list>
#include <set>
using namespace std;

class DataParser {
 public:
  explicit DataParser(const string &fileName);

  void processAndOutputDataToNDFormat();
  double *getProbabilityOfIndex(const list<tuple<int, string, string>> &samplesOfUniqueDevice);
  list<tuple<int, string, string>> extractSamplesByBoardID(const string &boardID);
  static void outputGraph(double *array, const string &addFileName);
  set<string> extractAllBoardIDs();

 private:
  void getDataFromCSV(const string &fileName);
  static tuple<int, string, string> getNextLineAndSplitIntoTokens(istream &str);
  static void writeDeviceDataIntoFile(const tuple<int, string, string> &data);
  static string commaSeparateData(const string &deviceData);

  list<tuple<int, string, string>> p_listOfSamples;
};

#endif //PUFPARSER_DATAPARSER_H_
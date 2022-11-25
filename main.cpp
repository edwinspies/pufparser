#include "DataParser.h"

#include <iostream>
#include <list>
#include <algorithm>
#include <tuple>
#include <chrono>

using namespace std;

int main() {
  set<string> allBoardIDs;
  list<tuple<int, string, string>> samplesOfUniqueBoardID;

  auto start = std::chrono::steady_clock::now();

  auto *dp = new DataParser("sram-result.csv");
  //output a file for each sample in the respective folder
  //Use this function to generate the data for the python script
  //dp->processAndOutputDataToNDFormat();

  //These lines of code will generate a picture for the first board ID that is in the list
  allBoardIDs = dp->extractAllBoardIDs();

  for (const string& s: allBoardIDs) {
	samplesOfUniqueBoardID = dp->extractSamplesByBoardID(s);
	DataParser::outputGraph(samplesOfUniqueBoardID);
  }

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

  delete(dp);
  return 0;
}



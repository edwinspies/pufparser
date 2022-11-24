#include "DataParser.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <list>
#include <algorithm>
#include <tuple>

using namespace std;

int main() {
  double *probArray;
  set<string> allBoardIDs;
  list<tuple<int, string, string>> samplesOfUniqueBoardID;

  auto *dp = new DataParser("sram-result.csv");

  //output a file for each sample in the respective folder
  //Use this function to generate the data for the python script
  //dp->processAndOutputDataToNDFormat();


  //These lines of code will generate a picture for the first board ID that is in the list
  allBoardIDs = dp->extractAllBoardIDs();
  samplesOfUniqueBoardID = dp->extractSamplesByBoardID(*allBoardIDs.begin());
  probArray = dp->getProbabilityOfIndex(samplesOfUniqueBoardID);

  DataParser::outputGraph(probArray, get<1>(samplesOfUniqueBoardID.front()));

  std::cout << "Hello, World!" << std::endl;
  return 0;
}



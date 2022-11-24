#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
#include <list>
#include <algorithm>

#include "DataParser.h"

using namespace std;

//vector structure
/*
 * 0. sample number
 * 1. date
 * 2. board type
 * 3. board id
 * 
 */

enum structure {
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

int main() {
  double *probArray;
  set<string> allBoardIDs;
  list<tuple<int, string, string>> samplesOfUniqueBoardID;

  auto *dp = new DataParser("sram-result.csv");

  //output a file for each sample in the respective folder
  //dp->processAndOutputDataToNDFormat();

  allBoardIDs = dp->extractAllBoardIDs();
  samplesOfUniqueBoardID = dp->extractSamplesByBoardID(*allBoardIDs.begin());
  probArray = dp->getProbabilityOfIndex(samplesOfUniqueBoardID);

  DataParser::outputGraph(probArray, get<1>(samplesOfUniqueBoardID.front()));

  std::cout << "Hello, World!" << std::endl;
  return 0;
}



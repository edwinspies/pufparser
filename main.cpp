#include "DataParser.h"

#include <iostream>
#include <list>
#include <algorithm>
#include <tuple>
#include <chrono>
#include <getopt.h>
#include <cstring>

using namespace std;

void operations(const string &filename, bool fileformat) {
  DataParser *dp;
  set<string> allBoardIDs;
  list<tuple<int, string, string, string>> samplesOfUniqueBoardID;

  auto start = std::chrono::steady_clock::now();

  dp = new DataParser(filename, fileformat);


  //output a file for each sample in the respective folder
  //Use this function to generate the data for the python script
  //dp->processAndOutputDataToNDFormat();

  //These lines of code will generate a picture for the first board ID that is in the list
  allBoardIDs = dp->extractAllBoardIDs();

  int i = 0;
  for (const string &s : allBoardIDs) {
	if (i > 0) break;
	samplesOfUniqueBoardID = dp->extractSamplesByBoardID(s);
	DataParser::outputGraph(samplesOfUniqueBoardID);
	i++;
  }

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

  delete (dp);
}

int main(int argc, char **argv) {
  int opt;
  string filename;
  bool altFileFormat = false;

  if (argc < 2) {
	cout << "Wrong input parameter, try again" << endl;
	cout << "Usage: " << argv[0] << "-f {CSV filename} -i {sergio or request}" << endl;
  }

  while ((opt = getopt(argc, argv, "f:i:")) != -1) {
	switch (opt) {
	  case 'f': filename = optarg;
		break;
	  case 'i':
		if (strcmp(optarg, "alternative") == 0) {
		  altFileFormat = true;
		  break;
		}
		if (strcmp(optarg, "request") == 0) {
		  altFileFormat = false;
		  break;
		}
	  default: /* '?' */
		fprintf(stderr, "Usage: %s -f {CSV filename} -i {request or alternative}\n", argv[0]);
		exit(-1);
	}
  }

  operations(filename, altFileFormat);

  return 0;
}



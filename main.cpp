#include "DataParser.h"

#include <iostream>
#include <list>
#include <algorithm>
#include <tuple>
#include <chrono>
#include <getopt.h>
#include <cstring>

using namespace std;
DataParser *dp;

string readInputFromConsole(const string& outputMessage) {
  std::cout << outputMessage << endl << flush;
  std::string input;
  std::getline(std::cin, input);

  return input;
}

void processUserinputForPictureGen(const string& inputBoardID, const string& inputAddress, const set<string>& boardIDs) {
  set<string> localBoardIDs;
  list<tuple<int, string, string, string>> samplesOfUniqueBoardID;

  //TODO convert input strings toUpper so they match the internal structure
  string inputBoardIDUpper = inputBoardID;
  string inputAddressUpper = inputAddress;


  std::transform(inputBoardIDUpper.begin()+2, inputBoardIDUpper.end(),inputBoardIDUpper.begin()+2, ::toupper);
  std::transform(inputAddressUpper.begin()+2, inputAddressUpper.end(),inputAddressUpper.begin()+2, ::toupper);

  //load localBoardIDs either with the given input or all board IDs
  if (inputBoardIDUpper != "all") {
	localBoardIDs.insert(inputBoardIDUpper);
  } else {
	localBoardIDs = boardIDs;
  }

  if (inputAddressUpper != "all") {
	for (const string &s : localBoardIDs) {
	  samplesOfUniqueBoardID = dp->extractSamplesByBoardID(s);
	  list<tuple<int, string, string, string>> l;
	  for(auto sample : samplesOfUniqueBoardID) {
		if (get<TUPLE_ADDRESS>(sample) == inputAddressUpper) {
		   l.emplace_back(sample);
		}
	  }
	  DataParser::outputGraph(l);
	}
  } else {
	for (const string &s : localBoardIDs) {
	  samplesOfUniqueBoardID = dp->extractSamplesByBoardID(s);
	  DataParser::outputGraph(samplesOfUniqueBoardID);
	}
  }


}

void operations(const string &filename, bool fileformat) {
  set<string> allBoardIDs;
  int actionNumber;
  string inputBoardID, inputAddress;
  size_t pos{};
  list<tuple<int, string, string, string>> samplesOfUniqueBoardID;

  dp = new DataParser(filename, fileformat);

//just let the user input commands and exit when he wants to
  while(true) {
	std::cout << "Enter number to execute action" << endl << flush;
	std::cout << "1 - Generate folder structure for specific board id & starting address to use it in the Python Histogram Script (not yet implemented)" << endl << flush;
	std::cout << "2 - Generate folder structure to use it in the Python Histogram Script" << endl << flush;
	std::cout << "3 - Generate averaged images for specific board id & starting address" << endl << flush;
	std::cout << "4 - Generate all averaged images per board id & starting address" << endl << flush;
	std::cout << "5 - exit" << endl << flush;
	std::string input;
	std::getline(std::cin, input);

	try {
	  actionNumber = stoi(input, &pos);
	}
	catch (std::invalid_argument const &ex) {
	  cout << "std::invalid_argument::what(): " << ex.what() << '\n';
	  actionNumber = -1;
	}

	auto start = std::chrono::steady_clock::now();
	allBoardIDs = dp->extractAllBoardIDs();

	cout << "no of board ids: " << allBoardIDs.size() << endl;

	switch(actionNumber) {
	  case 1: //specific folder structure
		//enter board id (maybe all)
		inputBoardID = readInputFromConsole("Please enter Board ID (or \"all\"");
		inputAddress = readInputFromConsole("Please enter Memory Address (or \"all\"");
		//enter address (maybe all)
		break;
	  case 2: //full folder structure
		dp->processAndOutputDataToNDFormat();
		break;
	  case 3: //specific image
		//enter board id (maybe all)
		inputBoardID = readInputFromConsole("Please enter Board ID (or \"all\")");
		//enter address (maybe all)
		inputAddress = readInputFromConsole("Please enter Memory Address (or \"all\")");
		processUserinputForPictureGen(inputBoardID, inputAddress, allBoardIDs);
		break;
	  case 4: //all board ids and all addresses
		for (const string &s : allBoardIDs) {
		  samplesOfUniqueBoardID = dp->extractSamplesByBoardID(s);
		  DataParser::outputGraph(samplesOfUniqueBoardID);
		}
		break;
	  case 5: //exit
		delete dp;
		return;
	  default:
		std::cout << "No right number was entered, please try again" << endl << flush;
	}

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
  }

}

int main(int argc, char **argv) {
  int opt;
  string filename;
  bool altFileFormat = false;

  if (argc < 2) {
	cout << "Wrong input parameter, try again" << endl;
	cout << "Usage: " << argv[0] << "-f {CSV filename} -i {request or alternative}" << endl;
	cout << "[request] ist equal to a file that was requested on the github page" << endl;
	cout << "[alternative] ist equal to a file that was provided via email" << endl;
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
  delete dp;

  return 0;
}



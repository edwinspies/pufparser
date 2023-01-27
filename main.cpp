#include "DataParser.h"

#include <iostream>
#include <list>
#include <algorithm>
#include <chrono>
#include <getopt.h>
#include <cstring>

using namespace std;
DataParser *dp;

string readInputFromConsole(const string &outputMessage) {
  std::string input;

  std::cout << outputMessage << std::endl << std::flush;
  std::getline(std::cin, input);

  return input;
}

void processUserInputForPictureGen(const string &inputBoardID,
								   const string &inputAddress,
								   const set<string> &boardIDs) {
  set<string> localBoardIDs;
  list<bitBlock> samplesOfUniqueBoardID;

  //load localBoardIDs either with the given input or all board IDs
  if (inputBoardID != "all") {
	localBoardIDs.insert(inputBoardID);
  } else {
	localBoardIDs = boardIDs;
  }

  if (inputAddress != "all") {
	for (const string &s : localBoardIDs) {
	  samplesOfUniqueBoardID = dp->extractSamplesByBoardID(s);
	  list<bitBlock> l;
	  for (const bitBlock &sample : samplesOfUniqueBoardID) {
		if (sample.address == inputAddress) {
		  l.emplace_back(sample);
		}
	  }
	  DataParser::outputGraph(l, -100);
	}
  } else {
	for (const string &s : localBoardIDs) {
	  samplesOfUniqueBoardID = dp->extractSamplesByBoardID(s);
	  DataParser::outputGraph(samplesOfUniqueBoardID, -100);
	}
  }

}

void operations(const string &filename, bool fileformat) {
  set<string> allBoardIDs;
  int actionNumber;
  string inputBoardID, inputAddress;
  size_t pos{};
  list<bitBlock> samplesOfUniqueBoardID;

  dp = new DataParser(filename, fileformat);

//just let the user input commands and exit when he wants to
  while (true) {
	std::cout << "Enter number to execute action" << endl << flush;
	std::cout
		<< "1 - Generate folder structure for specific board id & starting address to use it in the Python Histogram Script (not yet implemented)"
		<< std::endl << std::flush;
	std::cout << "2 - Generate folder structure to use it in the Python Histogram Script" << endl << flush;
	std::cout << "3 - Generate averaged images for specific board id & starting address" << endl << flush;
	std::cout << "4 - Generate all averaged images per board id & starting address" << endl << flush;
	std::cout << "5 - Generate all entropy files" << endl << flush;
	std::cout << "6 - Generate all probabilities per board id & starting address" << endl << flush;
	std::cout << "7 - exit" << endl << flush;
	std::string input;
	std::getline(std::cin, input);

	try {
	  actionNumber = stoi(input, &pos);
	}
	catch (std::invalid_argument const &ex) {
	  std::cout << "std::invalid_argument::what(): " << ex.what() << '\n';
	  actionNumber = -1;
	}

	auto start = std::chrono::steady_clock::now();
	allBoardIDs = dp->extractAllBoardIDs();

	std::cout << "no of board ids: " << allBoardIDs.size() << endl;

	switch (actionNumber) {
	  case 1: //specific folder structure
		inputBoardID = readInputFromConsole("Please enter Board ID (or \"all\"");
		inputAddress = readInputFromConsole("Please enter Memory Address (or \"all\"");
		//TODO
		break;
	  case 2: //full folder structure
		dp->processAndOutputDataToNDFormat();
		break;
	  case 3: //specific image
		inputBoardID = readInputFromConsole("Please enter Board ID (or \"all\")");
		inputAddress = readInputFromConsole("Please enter Memory Address (or \"all\")");
		processUserInputForPictureGen(inputBoardID, inputAddress, allBoardIDs);
		break;
	  case 4: //all images
		for (const string &s : allBoardIDs) {
		  samplesOfUniqueBoardID = dp->extractSamplesByBoardID(s);
		  DataParser::outputGraph(samplesOfUniqueBoardID, -100);
		}
		break;
	  case 5: dp->prepareBinaryEntropyOutput();
	  	break;
	  case 6: //output Probability files
		for (const string &s : allBoardIDs) {
		  samplesOfUniqueBoardID = dp->extractSamplesByBoardID(s);
		  DataParser::outputProbability(samplesOfUniqueBoardID);
		}
		break;
	  case 7: //exit
		delete dp;
		return;
	  default: std::cout << "No right number was entered, please try again" << endl << flush;
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
	std::cout << "Wrong input parameter, try again" << std::endl;
	std::cout << "Usage: " << argv[0] << "-f {CSV filename} -i {request or alternative}" << std::endl;
	std::cout << "[request] ist equal to a file that was requested on the github page" << std::endl;
	std::cout << "[alternative] ist equal to a file that was provided via email" << std::endl;
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



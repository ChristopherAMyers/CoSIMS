#pragma once
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <random>
#include "myStructs.h"

#ifdef WIN32
	#define osParam 0
#else
	#define osParam 1
#endif
using namespace std;


class fileParams
{
	int length_of_int(long long int); //returns number of characters in an integer
	bool isDecimal(const char*); //returns true if character is a decimal
	bool isInt(const char*); //returns true if character is an integer
	bool getBool(const char*); //converts a character array to boolean
	void setDefaultParams(string); //sets the default program parameters
	void importParams(string, char*, string); //import changes to parameters from input file
	void setFileInfo(string, string, char*); //sets file locations
	bool checkMolFileType(); //checks if a valid molecule file was presented

	string osSpacer; //uses a '/' for unix and a '\\' for windows
	string molFileExt; //molecule file extension

	char ifBuffer[400]; //reading file line buffer
	char ifPrintOut[10000]; //holds strings to later print to file
	userOptions changed = {}; //lost of booleans set of user decides to change default parapeters
	bool debugging = false; //if debugging features were enables

public:
	fileParams();
	~fileParams();

	bool interpretCMD(int, char *[], char*); //interprets command line options
	void printInfo(); //prints back the program parameters

	string getMolFileExt(); //returns molecule file extension
	bool getNameChanged(); //returns true if project name has changed from input file
	bool getChargeChanged(); //returns true if charge was changed from input file
};

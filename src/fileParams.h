/* 
* Copyright (C) 2019  Christopher A. Myers
* 
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
* If you have found this code usefull, please cite the research paper
* associated with this package.
*/

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

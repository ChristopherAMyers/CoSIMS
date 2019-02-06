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

#include "fileParams.h"


fileParams::fileParams()
{
	//osSpzcer is used to tell the difference between windows
	//and linux file systems
	if (osParam) osSpacer = "/";
	else osSpacer = "\\";
}

fileParams::~fileParams()
{}

//Returns the number of characters in an integer
int fileParams::length_of_int(long long int num)
{
  int lengthCount = 0;
  for(; num != 0; num /= 10, lengthCount++);
  return lengthCount;
}

//This function determines if a char array is a decimal or not
bool fileParams::isDecimal(const char *inChar)
{
	for (int i = 0; i < (signed)strlen(inChar); i++)
		if (!isdigit(inChar[i]) && inChar[i] != '.')
			return false;
	return true;
}

//This function determines if a char array is an integer or not
bool fileParams::isInt(const char *inChar)
{
    for (int i = 0; i < (signed)strlen(inChar); i++)
            if(!isdigit(inChar[i]))
                    return false;
    return true;
}

//This function converts a character array to a boolean
bool fileParams::getBool(const char *inChar)
{
	bool rtnValue;

	string inStr = inChar;
	::transform(inStr.begin(), inStr.end(), inStr.begin(), ::tolower);
	stringstream ss;
	ss << inStr;
	if (isInt(inChar))
		ss >> noboolalpha >> rtnValue;
	else
		ss >> boolalpha >> rtnValue;

	return rtnValue;
}

bool fileParams::interpretCMD(int argc, char *argv[], char* currentPath)
{
	string molFile = "", inputFile = "", atomFile = "";
        time_t now = time(0); // current date/time based on current system
        char* dt = ctime(&now); // convert now to string form
        char buffer[100];
        char printOut[1000] = "";
        sprintf(buffer, " Starting Program \n");
        strcat(printOut, buffer);
        sprintf(buffer, " Current Date/Time: %s\n", dt);
        strcat(printOut, buffer);
        sprintf(buffer, " Importing and assigning molecule information...\n");
        strcat(printOut, buffer);

	//determins if user changed settings are provided
	//set default to false
	changed = {};

	//loop through each command line arguments
	for (int i = 1; i < argc; i+= 2)
	{

		//checks if the option is in a form "-a", where a is a single character
		if (argv[i][0] == '-' && argv[i+1][0] != '-' && strlen(argv[i]) ==2)
		{
			//checks number of threads
			///CHANGE ME
			if (argv[i][0] == 'n')
			{
				//if the next value is infact an integral number of threads
				if (isInt(argv[i + 1]) && atoi(argv[i + 1]) != 0)
					numThreads_ext = atoi(argv[i + 1]);
			}

			//checks the input-file name
			else if (argv[i][1] == 'i')
				inputFile = argv[i + 1];

			//checks the molecule file name
			else if (argv[i][1] == 'm')
				molFile = argv[i + 1];

			//checks the forcefield atom types file name
			else if (argv[i][1] == 'a')
				atomFile = argv[i + 1];

			//checks the number of openMP threads
			else if (argv[i][1] == 'n')
			{
				if (isInt(argv[i + 1]))
				{
					numThreads_ext = atoi(argv[i + 1]);
					changed.threads = true;
				}
				else
					return false;
			}

			else if (argv[i][1] == 's')
			{
				if (isInt(argv[i + 1]))
				{
					seed_ext = atoll(argv[i + 1]);
					changed.seed = true;
				}
				else
					return false;
			}

			//return error if incorrect argument
			else return false;
		}

		//return error is not in the proper format
		else return false;
	}


	//if an input file was provided at program call
	if (inputFile.empty())
		setDefaultParams(molFile);
	//else, apply default values
	else
		importParams(inputFile, currentPath, molFile);



	//set the remaining file settings
	setFileInfo(molFile, atomFile, currentPath);

        cout << printOut;
        fileStream_ext.open(omegaFilePath_ext.c_str(), std::ios_base::out);
        fileStream_ext << endl << printOut;
        strcpy(printOut, "");

        cout << ifPrintOut;
        fileStream_ext << ifPrintOut;

        sprintf(printOut, " Importing files complete.\n");
        cout << printOut;
        fileStream_ext << printOut;
        fileStream_ext.close();



	//if no molecule file was detected, return with error
	//this must be done last so that the error can be
	//printed to the output file as well as the console
	if (molFile.empty())
	{
		char buffer[100];
		sprintf(buffer, "ERROR: No molecule file was given. Check command-line aguments\n");
		cout << buffer;
		fileStream_ext.open(omegaFilePath_ext, ios_base::app);
		fileStream_ext << buffer;
		fileStream_ext.close();
		return false;
	}
	if(checkMolFileType())
	  return true;
	else
	  return false;

}

void fileParams::setDefaultParams(string molFileName)
{
  string tempName = "";
          bool foundDot = false;

		  //if no name is specified in input file, then choose the molecule file name
		  //the name will be the text before the first period of the molecule file name
		  //loop through each character of the file name until it reaches a period
		  //and then break if at the end of the file name or reaches a '\' or '/'
          for (int i = molFileName.size() - 1; i >= 0; i--)
          {
                  if ((molFileName.compare(i, 1, "/") == 0) || molFileName.compare(i, 1, "\\") == 0)
                          break;
                  if (molFileName.compare(i, 1, ".") == 0)
                  {
                          foundDot = true;
                          tempName = "";
                  }
                  else
                          tempName += molFileName.at(i);
          }
          reverse(tempName.begin(), tempName.end());

          //apply default values
          name_ext = tempName; //project name
          mass_ext = 4; //reduced mass of system 
          subDir_ext = ""; //output directory
          dt_ext = 0.01; //trajectory integration time step in pico-sed
          temp_ext = 298; //temperature of system in Kelvin
          trajNum_ext = 50000; //number of trajectories per integration
          iterations_ext = 10; //number of CCS integrals to average over
		  if(!changed.threads)
			numThreads_ext = 1; //number of OpenMP threads to run on
		  if (!changed.seed)
		  {
			  random_device rd;
			  seed_ext = rd(); //this will be applied later in monteTraj.cpp
		  }
          printMol_ext = false; //Debugging: print rotated/centered molecule
          printTraj_ext = false; //Debugging: 
          printTrajOnFail_ext = false;
          printData_ext = false;
          printDataOnFail_ext = false;
          dispersionCutoff_ext = true;
          dispersion_radius_ext = 0; //this will be applied in molecule.cpp
          multipole_ext = true;
          multipole_radius_ext = 0; //this will be applied in molecule2.cpp
          multipole_order_ext = 1;
          printRate_ext = 20000;
		  displayProg_ext = true;
          maxClusterSize_ext = 4;
          totalCharge_ext = 0;
          charge_ext = false;
          gMax_ext = 35; //this will be applied later
          gMin_ext = 5; //this will be applied later
          projection_ext = true;
}

void fileParams::importParams(string fileName, char *currentPath, string molFileName)
{
        //holds the parameter name and value
        string option;
        string value;
        string name;
        string inputFile;
        ifstream inFile;

        strcpy(ifPrintOut, "");

        string dblParam;
        bool warningFlag = false;

        //holds the line number
        int lineNum = 1;

	if (!fileName.compare(0, 1, osSpacer))
		inputFile = fileName;
	else
		inputFile = currentPath + osSpacer + fileName;

	sprintf(ifBuffer, " Reading input file: %s\n", inputFile.c_str());
	strcat(ifPrintOut, ifBuffer);
	inFile.open(inputFile.c_str());

	//apply default values
	setDefaultParams(molFileName);


	string line;
	
	while(getline(inFile, line))
	{

		istringstream stream (line);
		
		if (!line.empty())
		{
			option.clear();
			value.clear();
			stream >> option >> value;
			warningFlag = false;
			if (value.empty())
			{
				sprintf(ifBuffer, "\t ERROR: Line number %i in %s is incomplete;\n"
					"\t This line wil be ignored.\n\n"
					, lineNum, fileName.c_str());
				strcat(ifPrintOut, ifBuffer);
			}

			else if (option[0] != '#')
			{
				if (option.compare("name") == 0)
				{
					if (changed.name) { warningFlag = true; dblParam = "NAME"; }
					name_ext = value.c_str();
					changed.name = true;
				}
				else if (option.compare("dir") == 0)
				{
					if (changed.dir) { warningFlag = true; dblParam = "DIR"; }
					subDir_ext = value.c_str();
					changed.dir = true;
				}
				else if (option.compare("dt") == 0)
				{
					if (changed.dt) { warningFlag = true; dblParam = "DT"; }
					dt_ext = atof(value.c_str());
					changed.dt = true;
				}
				else if (option.compare("temp") == 0)
				{
					if (changed.temp) { warningFlag = true; dblParam = "TEMP"; }
					temp_ext = atof(value.c_str());
					changed.temp = true;
				}
				else if (option.compare("traj") == 0)
				{
					if (changed.traj) { warningFlag = true; dblParam = "TRAJ"; }
					trajNum_ext = atof(value.c_str()) * 1000;
					changed.traj = true;
				}
				else if (option.compare("iter") == 0)
				{
					if (changed.iter) { warningFlag = true; dblParam = "ITER"; }
					iterations_ext = atof(value.c_str());
					changed.iter = true;
				}
				else if (option.compare("threads") == 0)
				{
					if (changed.threads) { warningFlag = true; dblParam = "THREADS"; }
					numThreads_ext = (int)atof(value.c_str());
					changed.threads = true;
				}
				else if (option.compare("seed") == 0)
				{
					if (changed.seed) { warningFlag = true; dblParam = "SEED"; }
					seed_ext = atoll(value.c_str());
					changed.seed = true;
				}
				else if (option.compare("gmax") == 0)
				{
					if (changed.gMax) { warningFlag = true; dblParam = "GMAX"; }
					gMax_ext = atof(value.c_str());
					changed.gMax = true;
					debugging = true;
				}
				else if (option.compare("gmin") == 0)
				{
					if (changed.gMin) { warningFlag = true; dblParam = "GMIN"; }
					gMin_ext = atof(value.c_str());
					changed.gMin = true;
					debugging = true;
				}
				else if (option.compare("print_mol") == 0)
				{
					if (changed.print_mol) { warningFlag = true; dblParam = "PRINT_MOL"; }
					printMol_ext = stoll(value.c_str());
					changed.print_mol = true;
					debugging = true;
				}
				else if (option.compare("proj") == 0)
				{
					if (changed.proj) { warningFlag = true; dblParam = "PROJ"; }
					projection_ext = (int)atof(value.c_str());
					changed.proj = true;
					debugging = true;
				}
				else if (option.compare("print_traj") == 0)
				{
					if (changed.print_traj) { warningFlag = true; dblParam = "PRINT_TRAJ"; }
					printTraj_ext = getBool(value.c_str());
					changed.print_traj = true;
					debugging = true;
				}
				else if (option.compare("print_traj_fail") == 0)
				{
					if (changed.print_traj_fail) { warningFlag = true; dblParam = "PRINT_TRAJ_FAIL"; }
					printTrajOnFail_ext = getBool(value.c_str());
					changed.print_traj_fail = true;
					debugging = true;
				}
				else if (option.compare("print_data") == 0)
				{
					if (changed.print_data) { warningFlag = true; dblParam = "PRINT_DATA"; }
					printData_ext = getBool(value.c_str());
					changed.print_data = true;
					debugging = true;
				}
				else if (option.compare("print_data_fail") == 0)
				{
					if (changed.print_data_fail) { warningFlag = true; dblParam = "PRINT_DATA_FAIL"; }
					printDataOnFail_ext = getBool(value.c_str());
					changed.print_data_fail = true;
					debugging = true;
				}
				else if (option.compare("dispersion_cutoff") == 0)
				{
					if (changed.dispersion_cutoff) { warningFlag = true; dblParam = "DISPERSION_CUTOFF"; }
					dispersionCutoff_ext = getBool(value.c_str());
					changed.dispersion_cutoff = true;
				}
				else if (option.compare("dispersion_radius") == 0)
				{
					if (changed.dispersion_radius) { warningFlag = true; dblParam = "DISPERSION_RADIUS"; }
					dispersion_radius_ext = atof(value.c_str());
					changed.dispersion_radius = true;
				}
				else if (option.compare("multipole_radius") == 0)
				{
					if (changed.multipole_radius) { warningFlag = true; dblParam = "MULTIPOLE_RADIUS"; }
					multipole_radius_ext = atof(value.c_str());
					changed.multipole_radius = true;
				}
				else if (option.compare("multipole_order") == 0)
				{
					if (changed.multipole_order) { warningFlag = true; dblParam = "MULTIPOLE_ORDER"; }
					multipole_order_ext = atoi(value.c_str());
					changed.multipole_order = true;
				}
				else if (option.compare("print_rate") == 0)
				{
					if (changed.print_rate) { warningFlag = true; dblParam = "PRINT_RATE"; }
					printRate_ext = atof(value.c_str());
					displayProg_ext = true;
					changed.print_rate = true;
				}
				else if (option.compare("multipole") == 0)
				{
					if (changed.multipole) { warningFlag = true; dblParam = "MULTIPOLE"; }
					multipole_ext = getBool(value.c_str());
					changed.multipole = true;
				}
				else if (option.compare("max_cluster_size") == 0)
				{
					if (changed.max_cluster_size) { warningFlag = true; dblParam = "MAX_CLUSTER_SIZE"; }
					maxClusterSize_ext = (int)atoi(value.c_str());
					changed.max_cluster_size = true;
				}
				else if (option.compare("charge") == 0)
				{
					if (changed.charge) { warningFlag = true; dblParam = "CHARGE"; }
					totalCharge_ext = atof(value.c_str());
					changed.charge = true;
					if (totalCharge_ext != 0) charge_ext = true;
					else charge_ext = false;
				}
				else
				{

					sprintf(ifBuffer, "\tERROR: Can not interpret option '%s' on \n\tline number %i in %s;\n"
						"\tThis line wil be ignored.\n\n"
						, option.c_str(), lineNum, fileName.c_str());
					strcat(ifPrintOut, ifBuffer);

				}

				if (warningFlag)
				{
					sprintf(ifBuffer, "\tWARNING: Parameter %s already specified at line number %i;\n"
						"\tMost recent value will be applied.\n\n", dblParam.c_str(), lineNum);
					strcat(ifPrintOut, ifBuffer);
				}
			}
		}
		lineNum++;

	}

	if (!changed.print_rate)
		printRate_ext = (int)(trajNum_ext / 5.0);
	else
		if (printRate_ext > trajNum_ext || printRate_ext == 0)
		{
			printRate_ext = trajNum_ext;
			displayProg_ext = false;
		}


	if(!subDir_ext.empty())
		if (subDir_ext.compare(subDir_ext.size() - 1, 1, osSpacer))
		{
			subDir_ext += osSpacer;
		}

	stringstream ss;
	ss << (int)totalCharge_ext;
	name_ext += "_" + ss.str();
}


void fileParams::printInfo()
{


  char buffer[200];
  char printOut[10000] = {};

  long long int width = max((long long)trajNum_ext, seed_ext);
  width = max((long long)iterations_ext, width);
  int length = length_of_int(width);
  length = max((int)name_ext.length(), length);
  length = max(length, 8);

  sprintf(buffer, "\n OmegaCal will use the following settings:\n");
  strcat(printOut, buffer);

  if(width > 20)
    sprintf(buffer, "\t Project name: %s\n", name_ext.c_str());
  else
    sprintf(buffer, "\t Project name: %*s\n", length, name_ext.c_str());
  sprintf(buffer, "\t Project name:                                              %*s\n", length, name_ext.c_str());
  strcat(printOut, buffer);

  sprintf(buffer, "\t Number of CCS integrals:                                   %*i\n", length, iterations_ext);
  strcat(printOut, buffer);

  sprintf(buffer, "\t Number of trajectories per CCS integral:                   %*d\n", length , trajNum_ext);
  strcat(printOut, buffer);

  sprintf(buffer, "\t Trajectory time step (pico sec):                           %*.5f\n", length, dt_ext);
  strcat(printOut, buffer);

  sprintf(buffer, "\t Temperature of the molecule system (Kelvin):               %*.2f\n", length, temp_ext);
  strcat(printOut, buffer);

  sprintf(buffer, "\t Random number seed for 64-bit Mersenne-Twister:            %*lli\n", length, seed_ext);
  strcat(printOut, buffer);

  sprintf(buffer, "\t Total charge of molecule (e):                              %*.3g\n", length, totalCharge_ext);
  strcat(printOut, buffer);


  if (charge_ext)
          if (multipole_ext)
                  if (dispersionCutoff_ext)
                    sprintf(buffer, "\t Multipole expansion and dispersion cut-off will be used.\n");
                  else
                    sprintf(buffer, "\t Multipole expansion and full dispersion interactions will be used.\n");
          else
                  if (dispersionCutoff_ext)
                    sprintf(buffer, "\t Exact electrostatics and dispersion cut-off will be used.\n");
                  else
                    sprintf(buffer, "\t Exact electrostatics and full dispersion interactions will be used.\n");
  else
          if (dispersionCutoff_ext)
                    sprintf(buffer, "\t No-charges present: molecule: only dispersion cutoff interactions will be used.\n");
          else
                    sprintf(buffer, "\t No-charges present: only (full) dispersion interactions will be used.\n");

  strcat(printOut, buffer);




  if(debugging)
  {
    sprintf(buffer, "\n Debugging features have been initiated:\n"
        "\t These features may not be complete or properly working.\n");
    strcat(printOut, buffer);

    if (printTraj_ext)
    {
            if (printTrajOnFail_ext)
              sprintf(buffer, "\t Program will print failed trajectories to\n");
            else
              sprintf(buffer, "\t Program will print successful trajectories to\n");
            strcat(printOut, buffer);
    }
    if (printData_ext)
    {

            if (printDataOnFail_ext)
              sprintf(buffer, "\t Program will print failed data to \n");
            else
              sprintf(buffer, "\t Program will print successful data to\n");
            strcat(printOut, buffer);

            sprintf(buffer, "\t %s\n", dataFilePath_ext.c_str());
            strcat(printOut, buffer);
    }
    if(printMol_ext)
    {
      sprintf(buffer, "\t Program will print rotated molecule coordinates to \n");
      strcat(printOut, buffer);

      sprintf(buffer, "\t %s\n", xyzFilePath_ext.c_str());
      strcat(printOut, buffer);
    }

  }

  sprintf(buffer, "\n Program input and output:\n");
  strcat(printOut, buffer);

  sprintf(buffer, "\t Program will print all outputs to:\n");
  strcat(printOut, buffer);

  sprintf(buffer, "\t %s\n", omegaFilePath_ext.c_str());
  strcat(printOut, buffer);

  sprintf(buffer, "\t Molecule will be imported from:\n");
  strcat(printOut, buffer);

  sprintf(buffer, "\t %s\n", molFilePath_ext.c_str());
  strcat(printOut, buffer);


  fileStream_ext.open(omegaFilePath_ext.c_str(), ios_base::app);
  cout << printOut << endl;
  fileStream_ext << printOut << endl;
  fileStream_ext.close();


}

void fileParams::setFileInfo(string molFilePath, string atomFilePath, char *currentPath)
{
	//load path to the molecule file by checking if given file
	//is a full directory, sub directory, not specified, or working directory
	if (!molFilePath.compare(0, 1, osSpacer))
		molFilePath_ext = molFilePath;
	else
		molFilePath_ext = currentPath + osSpacer + molFilePath;

	if(!atomFilePath.empty())
		if (!atomFilePath.compare(0, 1, osSpacer))
				atomFilePath_ext = atomFilePath;
		else
				atomFilePath_ext = currentPath + osSpacer + atomFilePath;

	//set the path to store the output files by checking if input files
	//is a full directory, sub directory, or not specified (default to working directory)
	string path;
	if (subDir_ext.empty())
		path = currentPath + osSpacer;
	else if (!subDir_ext.compare(0, 1, osSpacer))
	{
		path = subDir_ext;
	}
	else
		path = currentPath + osSpacer + subDir_ext;


	fileLocaiton_ext = currentPath + osSpacer;
	dataFilePath_ext = path + "data_" + name_ext + ".txt";
	omegaFilePath_ext = path + "omega_" + name_ext + ".txt";
	xyzFilePath_ext = path + "xyz_" + name_ext + ".xyz";
	printMolName_ext = path + "outMol_" + name_ext + ".xyz";


	fileStream_ext.open(omegaFilePath_ext.c_str(), ios::out);
	fileStream_ext.close();

}



bool fileParams::checkMolFileType()
{
  int length = molFilePath_ext.length();

  if(molFilePath_ext.at(length - 4) != '.')
    return false;

  molFileExt = molFilePath_ext.substr(length - 3, 3);
  transform(molFileExt.begin(), molFileExt.end(), molFileExt.begin(), ::tolower);
  if(molFileExt.compare("pdb") == 0 || molFileExt.compare("mfj") == 0)
    return true;
  else
    return false;

}

string fileParams::getMolFileExt()
{
  return molFileExt;
}

bool fileParams::getNameChanged()
{ return changed.name; }

bool fileParams::getChargeChanged()
{ return changed.charge; }

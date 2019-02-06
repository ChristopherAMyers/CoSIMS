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

#include "molecule2.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
using namespace std;


molecule2::molecule2()
{


		//importPDB();

}


molecule2::~molecule2(void)
{
  delete[] data;
}


void molecule2::assignDefaultFF()
{
  //all masses are taken from Handbook of Chemistry and Physics, 84th edition

    forceField tmpFF;

    tmpFF.name = "Hydrogen";
    tmpFF.symbol = "H";
    tmpFF.mass = 1.0079;
    tmpFF.sigma = 2.38;
    tmpFF.epsilon = 0.01498936;
    tmpFF.atomCount = 0;
    tmpFF.intMass = 1;
    fieldData.push_back(tmpFF);

    tmpFF.name = "Carbon";
    tmpFF.symbol = "C";
    tmpFF.mass = 12.010;
    tmpFF.sigma = 3.043;
    tmpFF.epsilon = 0.03090114;
    tmpFF.atomCount = 0;
    tmpFF.intMass = 12;
    fieldData.push_back(tmpFF);

    tmpFF.name = "Nitrogen";
    tmpFF.symbol = "N";
    tmpFF.mass = 14.006;
    tmpFF.sigma = 3.043;
    tmpFF.epsilon = 0.03090114;
    tmpFF.atomCount = 0;
    tmpFF.intMass = 14;
    fieldData.push_back(tmpFF);

    tmpFF.name = "Oxygen";
    tmpFF.symbol = "O";
    tmpFF.mass = 15.999;
    tmpFF.sigma = 3.043;
    tmpFF.epsilon = 0.03090114;
    tmpFF.atomCount = 0;
    tmpFF.intMass = 16;
    fieldData.push_back(tmpFF);

    tmpFF.name = "Sodium";///
    tmpFF.symbol = "Na";
    tmpFF.mass = 22.98977;
    tmpFF.sigma = 3.043;
    tmpFF.epsilon = 0.03090114;
    tmpFF.atomCount = 0;
    tmpFF.intMass = 23;
    fieldData.push_back(tmpFF);

    tmpFF.name = "Silicon";
    tmpFF.symbol = "Si";
    tmpFF.mass = 28.086;
    tmpFF.sigma = 3.5;
    tmpFF.epsilon = 0.03113175;
    tmpFF.atomCount = 0;
    tmpFF.intMass = 28;
    fieldData.push_back(tmpFF);

    tmpFF.name = "Phosphorus";
    tmpFF.symbol = "P";
    tmpFF.mass = 30.97396;
    tmpFF.sigma = 3.5;
    tmpFF.epsilon = 0.03113175;
    tmpFF.atomCount = 0;
    tmpFF.intMass = 31;
    fieldData.push_back(tmpFF);

    tmpFF.name = "Sulfur";
    tmpFF.symbol = "S";
    tmpFF.mass = 32.07;
    tmpFF.sigma = 3.5;
    tmpFF.epsilon = 0.03113175;
    tmpFF.atomCount = 0;
    tmpFF.intMass = 32;
    fieldData.push_back(tmpFF);

    tmpFF.name = "Iron";
    tmpFF.symbol = "Fe";
    tmpFF.mass = 55.84;
    tmpFF.sigma = 3.5;
    tmpFF.epsilon = 0.03113175;
    tmpFF.atomCount = 0;
    tmpFF.intMass = 56;
    fieldData.push_back(tmpFF);
}


void molecule2::importFile()
{

  

  //open forceField data file
  ifstream fieldFile;
  int size = (signed)fieldData.size();
  fieldFile.open(atomFilePath_ext.c_str());
  fileStream_ext.open(omegaFilePath_ext.c_str(), ios_base::app);
  cout << endl;
  cout << "Atom parameters are located at: " << endl;
  cout << atomFilePath_ext.c_str() << endl;

  fileStream_ext << endl;
  fileStream_ext << "Atom parameters are located at: " << endl;
  fileStream_ext << atomFilePath_ext.c_str() << endl;
  fileStream_ext.close();

  int index;

  //temp variables to hold data from line of file
  string name, symbol, mass, sigma, epsilon, intMass;


  //read in data file line by line and assign values
  while(fieldFile >> name >> symbol >> mass >> sigma >> epsilon >> intMass)
  {
	//convert atom name and type to lowercase for easy comparison
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    //transform(symbol.begin(), symbol.end(), symbol.begin(), ::tolower);

	//check if atom types already exist
	//if so, overwrite that atom's properties
	//else, add a new atom type
    index = fieldData.size();
    size = (signed)fieldData.size();

    if(name.at(0) != '#')
    {
      for (int i = 0; i < size; i ++)
      {

        if(symbol.compare(fieldData[i].symbol) == 0 || (fieldData[i].intMass == atoi(intMass.c_str())))
        {
          index = i;
          break;
        }
      }
      if(index == size)
        {
          fieldData.push_back(forceField());
        }

      fieldData[index].name = name;
      fieldData[index].symbol = symbol;
      fieldData[index].mass = atof(mass.c_str());
      fieldData[index].sigma = atof(sigma.c_str());
      fieldData[index].epsilon = atof(epsilon.c_str());
      fieldData[index].atomCount = 0;
      fieldData[index].intMass = atoi(intMass.c_str());
    }
  }


}

void molecule2::setCenter()
{
	//sets the center of mass of molecule to the origin
	for (int i = 0; i < numAtoms; i++)
	{
		//moves the original coords to the origin
		data[i].coords.X -= (centerX);
		data[i].coords.Y -= (centerY);
		data[i].coords.Z -= (centerZ);

		//moves the coords used by other classes to the origin
		data[i].origCoords.X = data[i].coords.X;
		data[i].origCoords.Y = data[i].coords.Y;
		data[i].origCoords.Z = data[i].coords.Z;
	}
}

void molecule2::findCenter()
{
	//finds the center of mass of molecule
	centerX = 0; centerY = 0; centerZ = 0;
	for (int i = 0; i < numAtoms; i++)
	{
		centerX += data[i].mass*data[i].coords.X;
		centerY += data[i].mass*data[i].coords.Y;
		centerZ += data[i].mass*data[i].coords.Z;
	}
	centerX = centerX/totalMass;
	centerY = centerY/totalMass;
	centerZ = centerZ/totalMass;
	setCenter(); //sets the center of Mass to origin
}

void molecule2::rotate(double theta, double phi, double gamma)
{
	//NEW CODE

	//rotates the molecule by a given Euler angles theta, phi, gamma
	for (int i = 0; i < numAtoms; i ++)
	{
		data[i].coords.rotate(data[i].origCoords, theta, phi, gamma);
	}

}
void molecule2::makeXYZ(double x, double y, double z, double threadNum)
{

	///THIS IS FORM AN OLDER VERSION, USED FOR REFERENCE PURPOSES ONLY//
	if (threadNum != 0)
	{
		xyzFile.open(xyzFilePath_ext.c_str(), ::ios_base::out);
		xyzFile << (numAtoms + 1) << "\n";
		xyzFile << "OmegaCal rotated and centered molecule \n";

		for (int i = 0; i < numAtoms; i ++)
		{
			xyzFile << data[i].atom << "     " << data[i].coords.X << "       " << data[i].coords.Y << "       " << data[i].coords.Z << "\n";
		}
		xyzFile << "He     " << x << "       " << y << "       " << z << "\n";
		//cout << threadNum << "               " << y << "            " << z << endl;
		xyzFile.close();
	}
}

void molecule2::importTXT()
{
  //location of input file
    file.open("/network/rit/lab/ChenRNALab/awesomeSauce/drift/eclpise_WrkSpce/MobcalV3/molFiles/c180.txt");//open the file


    std::string temps; //temp stores first string in first colm of file
    float temp[3];//temp stores the data values in cols of file
    totalMass = 0;
    centerX = 0; centerY = 0; centerZ = 0;

    //used to grab data from file, then pops data into array
    queue<std::string> first;//from temps to front of first
    queue<float> xq, yq, zq;//from temp[n] to front of queue
    //********************************************************
    while (file >> temps >> temp[0] >> temp[1] >> temp[2])//takes data from data file into queues
    {
            first.push(temps); xq.push(temp[0]); yq.push(temp[1]); zq.push(temp[2]);
    }
    file.close();
    numAtoms = xq.size();

    data = new dataList[numAtoms];
    //elmName = lineElements.at(3);

    int stop = xq.size();
    for (int i = 0; i < stop; i++)//pops data from queues into arrays for easy re-use of data
    {
            data[i].atom = first.front();
            first.pop();

            data[i].coords.X = xq.front();
            xq.pop();

            data[i].coords.Y = yq.front();
            yq.pop();

            data[i].coords.Z = zq.front();
            zq.pop();

            //computes the total mass of the molecule
            //possibly impliment sql files as a library of
            //masses and charge values for other atoms
            switch(data[i].atom.at(0))
            {
                    case 'C':
                            data[i].mass = 12.0107;
                            //data[i].sigma = 3.043;
                            data[i].sigma = 2.771;
                            data[i].epsilon = .129;
                            break;
                    case 'N':
                            data[i].mass =  14.0067;
                            data[i].sigma = 2.771;
                            data[i].epsilon = 0.3;
                            break;
                    case 'O':
                            data[i].mass = 15.9994;
                            data[i].sigma = 2.771;
                            data[i].epsilon = 0.3;
                            break;
                    case 'P':
                            data[i].mass = 30.973762;
                            data[i].sigma = 2.771;
                            data[i].epsilon = 0.3;
                            break;
                    case 'H':
                            data[i].mass = 1.00794;
                            data[i].sigma = 2.771;
                            data[i].epsilon = 0.3;
            }
            centerX += data[i].mass*data[i].coords.X; //total center of massX
            centerY += data[i].mass*data[i].coords.Y; //total center of massY
            centerZ += data[i].mass*data[i].coords.Z; //total center of massZ
            totalMass += data[i].mass;

            centerX /= totalMass; centerY /= totalMass; centerZ /= totalMass;

    }

    setCenter(); //puts the molecule crtOfMass at origin
    cout << "done" << endl;
}

void molecule2::importPQR()
{
  file.open(molFilePath_ext.c_str(), ios::in);//open the file
  importFile();
  cout << "path to PDB FILE: " << molFilePath_ext.c_str() << endl;

  //queues to store atom info
  queue<double> xCoord, yCoord, zCoord, charge;
  queue <string>type;

  //reads line of file, used for extracting coordinates
  const int SIZE = 800;
  string line;
  stringstream sstream;

  //used to check if the line is an ATOM line
  char ATOM[78] = {'A', 'T', 'O', 'M'};
  char HETATM[78] = {'H', 'E', 'T', 'A', 'T', 'M'};

  //read the first line of the file
  //getline(file, line);
  //sstream.str(line);

  cout << "IN FILe " << line << endl;
  cout << "IN FILe " << sstream.str() << endl;
  vector<string> lineElements;

  bool chainIDPresent = false;
  int lineCount = 0;

  while(!file.eof())
  {
    //net next line in the file
    lineCount ++;
    getline(file, line);
    istringstream ss(line);
    cin.get();

    //parse each string in the line
    string s;
    while(ss >> s)
    {

      lineElements.push_back(s);
      cout << "LINE: " << s << endl;
    }

    //check if it is an atom
    if(lineElements.at(0).compare(ATOM) == 0 || lineElements.at(0).compare(HETATM) == 0 )
    {
      if(lineElements.size() != 10 || lineElements.size() != 11)
      {
        cout << " Error in reading molecule file: Number of elements in line " << lineCount
        << " exceeds the allowed number for PQR files" << endl;
        break;
      }

      lineElements.pop_back();
      charge.push(atof(lineElements.at(lineElements.size()).c_str()));

      lineElements.pop_back();
      charge.push(atof(lineElements.at(lineElements.size()).c_str()));

      lineElements.pop_back();
      charge.push(atof(lineElements.at(lineElements.size()).c_str()));

      lineElements.pop_back();
      charge.push(atof(lineElements.at(lineElements.size()).c_str()));
      lineElements.pop_back();

      type.push(lineElements.at(3));
    }
  }

  file.close();

  //for(int i = 0; i < (signed)xCoord.size(); i ++)
    //cout << "Elm: " << xCoord.[i] << "  " << yCoord.at(i) << "  " << zCoord.at(i) << "   " << type.at(i) << "\n";

}

bool molecule2::isInt(string inStr)
{

    //This function determins if a char array is an integer or not
    char *inChar = new char[inStr.length() + 1];
    strcpy(inChar, inStr.c_str());
    for (int i = 0; i < (signed)strlen(inChar); i++)
            if(!(isdigit(inChar[i]) || inChar[i] == ' '))
                    return false;
    return true;
}

bool molecule2::isDecimal(string inStr)
{
        //This function determines if a char array is a decimal or not
      char *inChar = new char[inStr.length() + 1];
      strcpy(inChar, inStr.c_str());
        for (int i = 0; i < (signed)strlen(inChar); i++)
                if (!isdigit(inChar[i]) && inChar[i] != '.')
                        return false;
        return true;
}

void molecule2::importMFJ(bool changedName, bool changedCharge)
{
  char buffer[200];
  char printOut[1000] = "";

  bool uniformCharge = false;
  bool calcCharge = false;
  string unit, chargeType, tempVal;
  string tempX, tempY, tempZ, tempMass, tempCharge;
  double scale, unitConv;
  queue<double> xCoord, yCoord, zCoord, pCharge;
  queue <string>type, fileLine;
  int tmpIntMass;


  assignDefaultFF();

  if(!atomFilePath_ext.empty())
                 importFile();

  file.open(molFilePath_ext.c_str(), ios::in);

  //determine if project name or charge has already been set form input file
  if(changedName)
  {
    sprintf(buffer, " Warning: Project name specified in both MFJ and input file.\n"
        " Program will use MFJ file name\n");
    strcat(printOut, buffer);
  }
  if(changedCharge)
    {
      sprintf(buffer, " Warning: Molecule charge specified in both MFJ and input file.\n"
          " Program will use MFJ molecule charge\n");
      strcat(printOut, buffer);
    }


  //get project name
  getline(file, tempVal);
  //remove spaces in line
  tempVal.erase(remove_if(tempVal.begin(), tempVal.end(), ::isspace), tempVal.end());
  if(tempVal.empty())
  {
    sprintf(buffer, " Warning: Project name in MFJ file is empty\n"
        " Using default name of %s\n", name_ext.c_str());
    strcat(printOut, buffer);
  }
  else
    name_ext = tempVal;


  //get number of coordinates (multiple coordinates
  //file >> tempVal;
  getline(file, tempVal);
  tempVal.erase(remove_if(tempVal.begin(), tempVal.end(), ::isspace), tempVal.end());
  if(atoi(tempVal.c_str()) > 1)
  {
    sprintf(buffer, " Warning: Multiple coordinate sets not yet supported.\n");
    strcat(printOut, buffer);
  }



  //get number of atoms in the molecule
  getline(file, tempVal);
  tempVal.erase(remove_if(tempVal.begin(), tempVal.end(), ::isspace), tempVal.end());
  if(tempVal.empty())
  {
    sprintf(buffer, " Error: Number of atoms not specified.\n");
    strcat(printOut, buffer);
  }
    numAtoms = atoi(tempVal.c_str());


  //get the units of the coordinates, either Angstrom ("ang")
  //or atomic units ("au")
  getline(file, tempVal);
  tempVal.erase(remove_if(tempVal.begin(), tempVal.end(), ::isspace), tempVal.end());
  transform(tempVal.begin(), tempVal.end(), tempVal.begin(), ::tolower);
  if(tempVal.empty() || !(tempVal.compare("ang") == 0 || tempVal.compare("au") == 0))
  {
    sprintf(buffer, " Warning: Coordinate units not specified; Assuming Angstroms.\n");
    strcat(printOut, buffer);
  }
  else
    unit = tempVal;


  //Get the type of charge distribution.
  //This program also allows integer values
  //to indicate a uniform charge
  getline(file, tempVal);
  tempVal.erase(remove_if(tempVal.begin(), tempVal.end(), ::isspace), tempVal.end());
  if(tempVal.empty())
  {
    sprintf(buffer, " Warning: Charge type not specified; Assuming no charge.\n");
    strcat(printOut, buffer);
    charge_ext = false;
  }
  else
    {
      chargeType = tempVal;
      if(chargeType.compare("none") == 0)
        charge_ext = false;
      else if(chargeType.compare("equal") == 0)
        {
          charge_ext = true;
          totalCharge_ext = 1.0;
          uniformCharge = true;
        }
      else if(isInt(tempVal.c_str()))
      {
        charge_ext = true;
        totalCharge_ext = atoi(tempVal.c_str());
        uniformCharge = true;
      }
      else if(chargeType.compare("calc") == 0)
      {
          calcCharge = true;
          charge_ext = true;
      }
      else
      {
        sprintf(buffer, " Warning: Charge type incorrectly specified; Assuming no charge.\n");
        strcat(printOut, buffer);
        charge_ext = false;
      }
    }


  //get correction factor
  scale = 1.0;
  getline(file, tempVal);
  tempVal.erase(remove_if(tempVal.begin(), tempVal.end(), ::isspace), tempVal.end());
  if(tempVal.empty())
  {
    sprintf(buffer, " Warning: Correction factor not given; assuming identity\n");
    strcat(printOut, buffer);
  }
  else
    if(isInt(tempVal) || isDecimal(tempVal))
      scale = atof(tempVal.c_str());
    else
    {
      sprintf(buffer, " Warning: Incorrect correction factor; assuming .\n");
      strcat(printOut, buffer);
    }

  cout << endl << printOut << endl;
  fileStream_ext.open(omegaFilePath_ext.c_str(), ::ios_base::app);
  fileStream_ext << endl << printOut << endl;
  fileStream_ext.close();

  unitConv = 1.0;
  if(unit.compare("au") == 0)
    unitConv = 0.529177210;

  //loop through each line and grap coordinate data
  double totalCharge = 0.0;
  double thisChg;
  while(file >> tempX >> tempY >> tempZ >> tempMass >> tempCharge)
    {
      //cout << tempX << "  " << tempY << "  " << tempZ << "  " << tempMass << "  " << tempCharge << endl;
      xCoord.push(atof(tempX.c_str())*scale*unitConv);
      yCoord.push(atof(tempY.c_str())*scale*unitConv);
      zCoord.push(atof(tempZ.c_str())*scale*unitConv);
      
	  fileLine.push("Line number " + to_string(xCoord.size()) + "\n");

      if(charge_ext)
      {
          if(uniformCharge)
          {
              thisChg = totalCharge_ext/(double)numAtoms;
              pCharge.push(thisChg);
          }
          else
          {
              thisChg = atof(tempCharge.c_str());
              pCharge.push(thisChg);
          }
          totalCharge += thisChg;
      }
      else
        pCharge.push(0.0);


      tmpIntMass = atoi(tempMass.c_str());

      for(int j = 0; j < (signed)fieldData.size(); j++)
      {
        if(tmpIntMass == fieldData[j].intMass)
        {
          type.push(fieldData[j].symbol);
          break;
        }
      }

    }


  file.close();
  if(calcCharge)
    totalCharge_ext = totalCharge;
  setMolProperties(xCoord, yCoord, zCoord, pCharge, type, fileLine);
}

void molecule2::importPDB()
{
	//ofstream molRead;
	//molRead.open("moleculeRead.txt");

        numH = 0; numO = 0; numN = 0; numC = 0; numP = 0;



        //location of input file
        assignDefaultFF();
        if(!atomFilePath_ext.empty())
               importFile();
        file.open(molFilePath_ext.c_str(), ios::in);//open the file


        std::string temps; //temp stores first string in first colm of file

        //queues to store atom info
        queue<double> xCoord, yCoord, zCoord;
        queue <string>type, fileLine;

        //reads line of file, used for extracting coordinates
        const int SIZE = 800;
        char tempRead[SIZE];
		int atomNameLength = 1;

        //used to check if the line is an ATOM line
        char ATOM[78] = {'A', 'T', 'O', 'M'};
		char HETATM[78] = { 'H', 'E', 'T', 'A', 'T', 'M' };

        //read the first line of the file
        file.getline(tempRead, SIZE, '\n');


        //temporarily stores the coordinate value of an atom
        //as a string, then converted to a float
        //atomType stores the atom type info
        char partX[9];
        char partY[9];
        char partZ[9];
        char atomType[10];


        while (!file.eof())
        {

          //checks if the line being read is an atom
          char checkATOM[10];
          char checkHETATM[10];

          strncpy(checkATOM, tempRead,6);
          strncpy(checkHETATM, tempRead, 6);

          checkATOM[4] = '\0';
          checkHETATM[6] = '\0';
          atomNameLength = 1;
		  //if this is an ATOM or a HETATM entry
          if (strcmp(checkATOM, ATOM) == 0 || strcmp(checkHETATM, HETATM) == 0)
          {
            temps = "";

            //if the line has enough space to hold a possible atom type
            if (strlen(tempRead) >= 78)
            {
                    memcpy(atomType, tempRead + 76, 2);
					atomType[2] = '\0';
					//molRead << tempRead << endl;
					//molRead << atomType[0] << "|" << atomType[1];
                    //if the atom type entry is empty, use the atom name instead
                    if (atomType[0] == ' ' && atomType[1] == ' ')
                    {
                            memcpy(atomType, tempRead + 12, 5);
                            atomType[5] = '\0';

                    }
                    else
                    {
                            atomType[2] = '\0';
                            atomNameLength = 2;

                    }
            }
            else
            {
                    memcpy(atomType, tempRead + 12, 5);
                    atomType[5] = '\0';
            }

            //Loop thorugh the segment of text and pick out the letters that
            //are used to describe that atom type. Note: If no atom type is
            //given in columns 76 and 77, then the program only uses the
            //first letter of the atom-name column as the atom type.
            //This obviously has limitations to the types of atoms found.
            for (int i = 0; i < (signed)strlen(atomType); i++)
                    if (isalpha(atomType[i]))
                    {

                            temps.append(1, atomType[i]);
                            if((signed)temps.length() == atomNameLength)
                                    break;

                    }

			//cout << "Temps: " << temps << endl;
            type.push(temps);


            //get what element the atom is
            //memcpy(atomType, tempRead + 12, 5);
            //atomType[5] = '\0';
            //type.push(atomType);

            //get x-coordinate
            memcpy(partX, tempRead+30, 8);
            partX[8] = '\0';
            xCoord.push(atof(partX));


            //get y-coordinate
            memcpy(partY, tempRead+38, 8);
            partY[8] = '\0';
            yCoord.push(atof(partY));

            //get z-coordinate
            memcpy(partZ, tempRead+46, 8);
            partZ[8] = '\0';
            zCoord.push(atof(partZ));

			fileLine.push(tempRead);
          }
          file.getline(tempRead, SIZE, '\n');
        }


		//molRead.close();

        file.close();

        queue<double> temp;
        setMolProperties(xCoord, yCoord, zCoord, temp, type, fileLine);
}

void molecule2::setMolProperties(queue<double> xCoord, queue<double> yCoord, queue<double> zCoord
	, queue<double> pCharge, queue<string> type, queue<string> fileLine)
{

	char buffer[1000];
	char printOut[1000] = "";
  numAtoms = xCoord.size();
  data = new dataList[numAtoms];
  string temp;
  bool typeFound;
                string atom;
                bool foundUpper = false;
                for (int i = 0; i < numAtoms; i ++)
                {

                    atom = type.front();
					typeFound = false;
                    for(int j = 0; j < (signed)fieldData.size(); j ++)
                    {
						
                        //transform(atom.begin(), atom.end(), atom.begin(), ::tolower);
                        if(fieldData[j].symbol.compare(atom) == 0)
                        {
						  typeFound = true;
                          data[i].mass = fieldData[j].mass;
                          data[i].sigma = fieldData[j].sigma;
                          data[i].epsilon = fieldData[j].epsilon;
                          data[i].sigAvg = fieldData[j].sigma;
                          data[i].epsAvg = fieldData[j].epsilon;
                          //convert from kcal/mol to u*ang^2*ps^-2
                          data[i].epsAvg *= calorie_ext / (avo_ext*protonMassKg_ext*10.0);
                          data[i].sigAvg6 = pow(data[i].sigAvg, 6);
                          data[i].sigAvg12 = data[i].sigAvg6*data[i].sigAvg6;
                          data[i].eps24_mass = data[i].epsAvg*24; //the mass part will be asssigned later
                          data[i].eps24 = data[i].epsAvg*24;
                          data[i].eps4 = data[i].epsAvg*4;

                          if(charge_ext)
                            {
                              if(pCharge.size() != 0)
                              {
                                data[i].charge = pCharge.front();
                                pCharge.pop();
                              }
                              else
                                data[i].charge = totalCharge_ext/(double)numAtoms;
                            }

                          else
                            data[i].charge = 0.0;

                          fieldData[j].atomCount++;
                          break;
                        }

                    }
					if (!typeFound)
					{
						fileStream_ext.open(omegaFilePath_ext.c_str(), ::ios_base::app);
						sprintf(buffer, "Warning: Force field does not contain atoms of type %s from molecule file line:\n%s\n"
							, atom.c_str(), fileLine.front().c_str());
						
						fileStream_ext << buffer;
						cout << buffer;
						fileStream_ext.close();
					}
					fileLine.pop();


                    data[i].atom = type.front();
                    type.pop();

                    data[i].coords.X = xCoord.front();
                    xCoord.pop();

                    data[i].coords.Y = yCoord.front();
                    yCoord.pop();

                    data[i].coords.Z = zCoord.front();
                    zCoord.pop();


                    centerX += data[i].mass*data[i].coords.X; //total center of massX
                    centerY += data[i].mass*data[i].coords.Y; //total center of mass Y
                    centerZ += data[i].mass*data[i].coords.Z; //total center of massZ
                    totalMass += data[i].mass;

                }

                            //finish the division of mass
                            mass_ext = (totalMass*heliumMassU_ext) / (totalMass + heliumMassU_ext);

                            centerX /= totalMass; centerY /= totalMass; centerZ /= totalMass;
                            for (int n = 1; n < numAtoms; n++)
                                    data[n].eps24_mass /= mass_ext;

                    setCenter(); //puts the molecule crtOfMass at origin

					if ((dispersionCutoff_ext && dispersion_radius_ext == 0) || (multipole_ext && multipole_radius_ext == 0))
					{
						double sigmaMax = 0;
						double minPct = 0.5;
						for (int n = 0; n < (signed)fieldData.size(); n++)
							if (fieldData[n].sigma > sigmaMax)
								sigmaMax = fieldData[n].sigma;

						//this is determined by approximating the integral of the LJ-potential
						//over the cut-off sphere for a contineous mass density with a VDW radius
						//of sigmaMax, and solving for the cut-off radius until a minimum 
						//percentage (100-minPct) of the total energy is accounted for
						

						if (dispersionCutoff_ext && dispersion_radius_ext == 0)
						{
							dispersion_radius_ext = sigmaMax*pow(150.0 / minPct, 1.0 / 3.0);
							sprintf(buffer, " Dispersion cut-off radius not specified in input file: \n\tUsing a cut-off radius of %6.2f\n"
								"\t Determined by largest VDW radius of %5.3f\n", dispersion_radius_ext, sigmaMax);
							strcat(printOut, buffer);
						}
						if (multipole_ext && multipole_radius_ext == 0)
						{
						        //multipole_radius_ext = sigmaMax*pow(150.0 / minPct, 1.0 / 3.0);
							multipole_radius_ext = 20.0;
							printf(buffer, " Multipole cut-off radius not specified in input file: \n\tUsing a cut-off radius of %6.2f\n"
								"\t Determined by largest VDW radius of %5.3f\n", multipole_radius_ext, sigmaMax);
							strcat(printOut, buffer);
						}

					}
					fileStream_ext.open(omegaFilePath_ext.c_str(), ::ios_base::app);
					fileStream_ext << printOut;
					fileStream_ext.close();
					cout << printOut;
                    //printMolStats();

}

void molecule2::printMolStats()
{
	char buffer[200];

  fileStream_ext.open(omegaFilePath_ext.c_str(), ios_base::app);
  stringstream output;

  output.precision(5);
  output << left << " Atom type breakdown:" << endl;
  output.fill(' ');
  //output.width(25);

  sprintf(buffer, "\t Name        Symbol    Mass    Sigma    Epsilon  IntMass  Count\n");
  output << buffer;



  //output << left << setw(15) << "\t Name" << setw(10) << "symbol"
      //<< setw(10) << "mass" << setw(10) << "count" << setw(10) << "sigma" << "epsilon" << endl;
  for(int i = 0; i < (signed)fieldData.size(); i ++)
  if(fieldData[i].atomCount != 0)
  {
	  sprintf(buffer, "\t %-12s  %2s  %8.2f  %8.4f  %8.5f    %3i    %5i\n"
		  , fieldData[i].name.c_str(), fieldData[i].symbol.c_str(), fieldData[i].mass
		  , fieldData[i].sigma, fieldData[i].epsilon, fieldData[i].intMass, fieldData[i].atomCount);
	  output << buffer;

	  /*
      output << "\t " << left << setw(15) << fieldData[i].name
      << setw(10) << fieldData[i].symbol
      << setw(10) << fieldData[i].mass
      << setw(10) << fieldData[i].atomCount
      << setw(10) << fieldData[i].sigma
      << setw(10) << fieldData[i].epsilon * calorie_ext / (avo_ext*protonMassKg_ext*10.0) << endl;
	  */
  }



  output << "\t Total number of atoms:   " << numAtoms << endl;
  output << "\t Total mass of molecule:  " << totalMass << " a.u." << endl;
  output << "\t Reduced mass:            " << mass_ext << " a.u." << endl;
  //output << "\tCenter of mass vector: (" << centerX/mass_ext << ", " << centerY/mass_ext << ", " << centerZ/mass_ext << " )\n";

  fileStream_ext << output.str();
  cout << output.str();
  fileStream_ext.close();
  setMultipole();

}

double molecule2::delta(int a, int b)
{
  if(a == b)
    return 1.0;
  else
    return 0.0;
}

void molecule2::setMultipole()
{
  ostringstream output;
  output << setw(20);
  output << setprecision(5);

  //find the dipole term
  setDipole();
  setQuadpole();
  setOctupole();

  //first find inertia tensor and rotate molecule along principle axes
  setInertia(true);



  //get maximum extend along each axis;
  double maxX, maxY, maxZ;
  double absX, absY, absZ;
  double radGy = 0;
  maxX = 0.0; maxY = 0.0; maxZ = 0.0;
  for(int i = 0; i < numAtoms; i ++)
  {
    absX = abs(data[i].coords.X);
    absY = abs(data[i].coords.Y);
    absZ = abs(data[i].coords.Z);
    if(absX > maxX) maxX = absX;
    if(absY > maxY) maxY = absY;
    if(absZ > maxZ) maxZ = absZ;
    radGy += data[i].coords.mag2();
  }
  radGy = sqrt(radGy/numAtoms);


  int width = 10;
  char buffer[200];

  output << "\n\n\t Total Mass of molecule (a.u.):  ";
  sprintf(buffer, "%8.4g\n", totalMass);
  output << buffer;

  output << "\t Radius of Gyration (ang):     ";
  sprintf(buffer, "%*.4f\n", width, radGy);
  output << buffer;

  output << "\n\t Center of Mass (a.u. Ang):\n";
  sprintf(buffer, "\t { %*.4f  %*.4f  %*.4f   }\n", width, centerX, width, centerY, width, centerZ);
  output << buffer;

  output << "\n\t Inertia Tensor (a.u. Ang^2):\n";
  sprintf(buffer, "\t { %*.4g  %*.4g  %*.4g   }\n", width, origInertia[0][0], width, origInertia[0][1], width, origInertia[0][2]);
  output << buffer;
  sprintf(buffer, "\t { %*.4g  %*.4g  %*.4g   }\n", width, origInertia[1][0], width, origInertia[1][1], width, origInertia[1][2]);
  output << buffer;
  sprintf(buffer, "\t { %*.4g  %*.4g  %*.4g   }\n", width, origInertia[2][0], width, origInertia[2][1], width, origInertia[2][2]);
  output << buffer;

  output << "\n\t Rotation matrix applied to molecule: \n";
  sprintf(buffer, "\t { %*.4f  %*.4f  %*.4f   }\n", width, rotMat[0][0], width, rotMat[0][1], width, rotMat[0][2]);
  output << buffer;
  sprintf(buffer, "\t { %*.4f  %*.4f  %*.4f   }\n", width, rotMat[1][0], width, rotMat[1][1], width, rotMat[1][2]);
  output << buffer;
  sprintf(buffer, "\t { %*.4f  %*.4f  %*.4f   }\n", width, rotMat[2][0], width, rotMat[2][1], width, rotMat[2][2]);
  output << buffer;

  if(charge_ext)
  {
    output << "\n\t Total Charge of molecule (e):   ";
    sprintf(buffer, "%*.4f\n", 8, totalCharge_ext);
    output << buffer;

    output << "\n\t Dipole Moment (e - Ang):\n";
    sprintf(buffer, "\t { %*.4f  %*.4f  %*.4f   }\n", width, dipole[0], width, dipole[1], width, dipole[2]);
    output << buffer;

    output <<    "\t Quadrupole Moment (e - Ang^2): \n";
	sprintf(buffer, "\t { %*.4f  %*.4f  %*.4f   }\n", width, quadpole[0][0], width, quadpole[0][1], width, quadpole[0][2]);
	output << buffer;
	sprintf(buffer, "\t { %*.4f  %*.4f  %*.4f   }\n", width, quadpole[1][0], width, quadpole[1][1], width, quadpole[1][2]);
	output << buffer;
	sprintf(buffer, "\t { %*.4f  %*.4f  %*.4f   }\n", width, quadpole[2][0], width, quadpole[2][1], width, quadpole[2][2]);
	output << buffer;


  }



  cout << output.str();
  fileStream_ext.open(omegaFilePath_ext.c_str(), ios_base::app);
  fileStream_ext << output.str();
  fileStream_ext.close();
}

void molecule2::setDipole()
{
  //make a temporary copy of the coordinates;
  double **coords = new double *[numAtoms];
  for (int i = 0; i < numAtoms; i++) coords[i] = new double[3];

  for(int i = 0; i < numAtoms; i ++)
  {
    coords[i][0] = data[i].coords.X;
    coords[i][1] = data[i].coords.Y;
    coords[i][2] = data[i].coords.Z;
  }

  //initialize the elements of the array
  for(int n = 0; n < 3; n ++)
        dipole[n] = 0.0;

  //perform the summation
  for(int i = 0; i < numAtoms; i ++)
  {
    for(int n = 0; n < 3; n ++)
      dipole[n] += data[i].charge*coords[i][n];
  }
}

void molecule2::setQuadpole()
{
  int i, a, b;

  //make a temporary copy of the coordinates;
  double **coords = new double *[numAtoms];
  for (int i = 0; i < numAtoms; i++) coords[i] = new double[3];
  double *r2 = new double[numAtoms];
  for(int i = 0; i < numAtoms; i ++)
  {
    coords[i][0] = data[i].coords.X;
    coords[i][1] = data[i].coords.Y;
    coords[i][2] = data[i].coords.Z;
    r2[i] = coords[i][0]*coords[i][0]
                       + coords[i][1]*coords[i][1]
                                     + coords[i][2]*coords[i][2];
  }

  //initialize the elements of the array
  for(a = 0; a < 3; a ++)
    for(b = 0; b <3; b ++)
        quadpole[a][b] = 0.0;

  //perform the summation
  for(i = 0; i < numAtoms; i ++)
  {
    for(a = 0; a < 3; a ++)
      for(b = 0; b < 3; b ++)
        quadpole[a][b] += data[i].charge*(3*coords[i][a]*coords[i][b] - r2[i]*delta(a, b));
  }
}

void molecule2::setOctupole()
{
  int i, a, b, c;

  //make a temporary copy of the coordinates;
  double **coords = new double *[numAtoms];
  for (int i = 0; i < numAtoms; i++) coords[i] = new double[3];
  double *r2 = new double[numAtoms];
  for(int i = 0; i < numAtoms; i ++)
  {
    coords[i][0] = data[i].coords.X;
    coords[i][1] = data[i].coords.Y;
    coords[i][2] = data[i].coords.Z;
    r2[i] = coords[i][0]*coords[i][0]
                       + coords[i][1]*coords[i][1]
                                     + coords[i][2]*coords[i][2];
  }

  //initialize the elements of the array
  for(a = 0; a < 3; a ++)
    for(b = 0; b < 3; b ++)
      for(int c = 0; c < 3; c ++)
        octupole[a][b][c] = 0.0;

  //perform the summation
  for(i = 0; i < numAtoms; i ++)
  {
    for(a = 0; a < 3; a ++)
      for(b = 0; b < 3; b ++)
        for(int c = 0; c < 3; c ++)
          octupole[a][b][c] += data[i].charge*(5*coords[i][a]*coords[i][b]*coords[i][c]
                                - 3*r2[i]*(coords[i][a]*delta(b, c) + coords[i][b]*delta(a, c)
                                    + coords[i][c]*delta(a, b)));
  }
}

void molecule2::setInertia(bool rotate)
{
  int i, a, b, maxcoord;
  vector3D<double> temp1, temp2, temp3;
  double dotProd, temp, maxMag, phi, theta, gamma;
  char buffer[300];

  //make a temporary copy of the coordinates;
  double **coords = new double *[numAtoms];
  for (int i = 0; i < numAtoms; i++) coords[i] = new double[3];
  double *r2 = new double[numAtoms];
  for(int i = 0; i < numAtoms; i ++)
  {
    coords[i][0] = data[i].coords.X;
    coords[i][1] = data[i].coords.Y;
    coords[i][2] = data[i].coords.Z;
    r2[i] = coords[i][0]*coords[i][0]
                       + coords[i][1]*coords[i][1]
                                     + coords[i][2]*coords[i][2];
  }

  //initialize the elements of the array
  for(a = 0; a < 3; a ++)
    for(b = 0; b <3; b ++)
        inertia[a][b] = 0.0;

  //perform the summation and szve original inertia atrix
  for(i = 0; i < numAtoms; i ++)
  {
    for(a = 0; a < 3; a ++)
      for(b = 0; b < 3; b ++)
      {
        inertia[a][b] += data[i].mass*(r2[i]*delta(a, b) - coords[i][a]*coords[i][b]);
        origInertia[a][b] = inertia[a][b];
      }
  }


  //diagonalize the inertia tensor
  getEigenvalues(inertia, inertiaValues);
  getEigenvectors(inertia, inertiaValues, inertiaVectors);



  temp1.setEqual(inertiaVectors[0][0], inertiaVectors[0][1], inertiaVectors[0][2]);
  temp2.setEqual(inertiaVectors[1][0], inertiaVectors[1][1], inertiaVectors[1][2]);
  temp3.setEqual(inertiaVectors[2][0], inertiaVectors[2][1], inertiaVectors[2][2]);
  dotProd = abs(temp1.dot(temp2));
  dotProd = max(abs(temp1.dot(temp3)), dotProd);
  dotProd = max(abs(temp2.dot(temp3)), dotProd);


  if(dotProd > 1E-6)
  {
    sprintf(buffer, "\nFailed to find proper eigenvectors of Inertia tensor:"
                    "\nAtom with longest vector magnitude will be aligned"
                    "\nwith z-axis and the atom with greatest XY-projection "
                    "\nwill then be aligned along the x-axis.\n");

    fileStream_ext.open(omegaFilePath_ext.c_str(), ::ios_base::app);
    cout << buffer;
    fileStream_ext << buffer;
    fileStream_ext.close();


    maxMag = 0;
    for(i = 0; i < numAtoms; i ++)
    {
      temp = data[i].coords.mag();
      if(temp > maxMag)
        {
          maxMag = temp;
          maxcoord = i;
        }
    }
    theta = acos(data[maxcoord].coords.Z/maxMag);
    phi = atan2(data[maxcoord].coords.Y, data[i].coords.X);

    for(i = 0; i < numAtoms; i ++)
    {
      data[i].coords.rotateZ(-phi);
      data[i].coords.rotateY(-theta);
    }

    maxMag = 0;
    for(i = 0; i < numAtoms; i ++)
    {
      temp = sqrt(pow(data[i].coords.X, 2) + pow(data[i].coords.Y, 2));
      if(temp > maxMag)
        {
          maxMag = temp;
          maxcoord = i;
        }
    }

    gamma = acos(data[maxcoord].coords.X/maxMag);
    for(i = 0; i < numAtoms; i ++)
      data[i].coords.rotateZ(-gamma);

  }
  else
  {

    //calculate new inertia matrix and save rotation matrix
    for(int n = 0; n < 3; n ++)
      for(int m = 0; m < 3; m ++)
      {
        if(m == n)
          inertia[n][m] = inertiaValues[n];
        else
          inertia[m][n] = 0.0;
        rotMat[n][m] = inertiaVectors[n][m];
      }

    if(rotate)
	for (int i = 0; i < numAtoms; i++)
		data[i].coords.rotate2(inertiaVectors);
  }


  if(rotate)
  {




        //find direction of longest extent
        maxX = 0; maxY = 0; maxZ = 0;
        double tmpCoord;
        for(int i = 0; i < numAtoms; i ++)
        {

          if(abs(data[i].coords.X) > abs(maxX))
                maxX = abs(data[i].coords.X);

          if(abs(data[i].coords.Y) > abs(maxY))
                maxY = abs(data[i].coords.Y);

          if(abs(data[i].coords.Z) > abs(maxZ))
                maxZ = abs(data[i].coords.Z);
        }


        //if longest extent is along x-axis,
        //rotate about y-axis -pi/2
        if(maxX > maxY && maxX > maxZ)
        {
          for(int i = 0; i < numAtoms; i ++)
          {
            tmpCoord = data[i].coords.X;
            data[i].coords.X = -data[i].coords.Z;
            data[i].coords.Z = tmpCoord;
          }
        }

        //if longest extent is along y-axis,
        //rotate about x-axis -pi/2
        else if(maxY > maxX && maxY > maxZ)
        {
          for(int i = 0; i < numAtoms; i ++)
          {
            tmpCoord = data[i].coords.Y;
            data[i].coords.Y = data[i].coords.Z;
            data[i].coords.Z = -tmpCoord;
          }
        }

        //reset maximum extent values
        maxX = 0; maxY = 0; maxZ = 0;
        for(int i = 0; i < numAtoms; i ++)
          {

            if(abs(data[i].coords.X) > abs(maxX))
              maxX = data[i].coords.X;

            if(abs(data[i].coords.Y) > abs(maxY))
                  maxY = data[i].coords.Y;

			if (abs(data[i].coords.Z) > abs(maxZ))
			{
				maxZ = data[i].coords.Z;
			}
          }
    }

  

}

/* This function finds the exact eigenvalues for a given 3x3 matrix A
 * The characteristic equation is solved exactly to find the eigenvalues,
 * and then these are used to find the eigenvectors of the vector.
 * A 3x3 matrix of orthogonal eigenvectors are returned.
 */
void molecule2::getEigenvalues(double A[3][3], double (&eigenvalues)[3])
{
    //double eigenvalues[3];
  //eigenvalues = new double[3];
    double a0, a1, a2;
    double p, q, R, Q, D, a, b, phase, mag;

    //roots of the characteristic equation
    a0 = A[0][0]*A[1][1]*A[2][2] -
        A[0][0]*A[1][2]*A[2][1] -
        A[1][1]*A[0][2]*A[2][0] -
        A[2][2]*A[0][1]*A[1][0] +
        A[0][1]*A[1][2]*A[2][0] +
        A[0][2]*A[2][1]*A[1][0];

    a1 = A[0][1]*A[1][0] + A[0][2]*A[2][0] + A[1][2]*A[2][1] - A[0][0]*A[1][1] - A[1][1]*A[2][2] -
         A[2][2]*A[0][0];

    a2 = A[0][0] + A[1][1] + A[2][2];
    a0 = -a0; a1 = -a1; a2 = -a2;

    //solve the cubic equation for 3 real roots
    p = (3*a1 - a2*a2)/3;
    q = (9*a1*a2 - 27*a0 - 2*a2*a2*a2)/27;
    Q = p/3; R = q/2; D = Q*Q*Q + R*R;

    if(D< 0)
    {
      a = R; b = sqrt(-D);
    }
    else
    {
      a = R + sqrt(D);
	  b = 0;
    }
    phase = atan2(b, a);
    mag = pow(a*a + b*b,1.0/6.0);

    eigenvalues[0] = -(a2/3) + mag*2*cos(phase/3.0);
    eigenvalues[1] =  -(a2/3) - mag*(cos(phase/3.0)
        +sqrt(3)*sin(phase/3.0) );
    eigenvalues[2] =  -(a2/3) - mag*(cos(phase/3.0)
        -sqrt(3)*sin(phase/3.0) );
}

void molecule2::getEigenvectors(double A[3][3], double lambda[3], double (&eigenvectors)[3][3])
{
    //double eigenvectors[3][3];
    double det, x, y, z, normFactor;
	double temp[3][3];

	det = lambda[0] * lambda[0] - lambda[0] * A[1][1] - A[1][2] * A[1][2] - lambda[0] * A[2][2] + A[1][1] * A[2][2];
	y = (lambda[0] * A[1][0] - A[2][2] * A[1][0] + A[1][2] * A[2][0]) / det;
	z = (A[2][1] * A[1][0] + lambda[0] * A[2][0] - A[1][1] * A[2][0]) / det;
	normFactor = sqrt(1 + y*y + z*z);
	eigenvectors[0][0] = 1.0 / normFactor;
	eigenvectors[0][1] = y / normFactor;
	eigenvectors[0][2] = z / normFactor;

	det = lambda[1] * lambda[1] - lambda[1] * A[0][0] - A[0][2] * A[0][2] - lambda[1] * A[2][2] + A[0][0] * A[2][2];
	x = (lambda[1] * A[0][1] + A[0][2] * A[2][1] - A[0][1] * A[2][2]) / det;
	z = (A[0][1] * A[0][2] + (lambda[1] - A[0][0])*A[2][1]) / det;
	normFactor = sqrt(1 + x*x + z*z);
	eigenvectors[1][0] = x / normFactor;
	eigenvectors[1][1] = 1 / normFactor;
	eigenvectors[1][2] = z / normFactor;

	eigenvectors[2][0] = eigenvectors[0][1] * eigenvectors[1][2] - eigenvectors[0][2] * eigenvectors[1][1];
	eigenvectors[2][1] = eigenvectors[0][2] * eigenvectors[1][0] - eigenvectors[0][0] * eigenvectors[1][2];
	eigenvectors[2][2] = eigenvectors[0][0] * eigenvectors[1][1] - eigenvectors[0][1] * eigenvectors[1][0];


	vector3D<double> temp1, temp2, temp3;
	temp1.setEqual(eigenvectors[0][0], eigenvectors[0][1], eigenvectors[0][2]);
	temp2.setEqual(eigenvectors[1][0], eigenvectors[1][1], eigenvectors[1][2]);
	temp3.setEqual(eigenvectors[2][0], eigenvectors[2][1], eigenvectors[2][2]);

	/*for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			temp[i][j] = eigenvectors[i][j];

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			eigenvectors[i][j] = temp[j][i];*/



	/*
    for(int n = 0; n < 3; n ++)
    {
      //det = lambda[n]*A[0][2] - A[0][2]*A[1][1] + A[0][1]*A[1][2];
      //y = (A[0][2]*A[1][0] + lambda[n]*A[1][2] - A[0][0]*A[1][2])/det;
      //z = (lambda[n]*lambda[n] - lambda[n]*(A[0][0] + A[1][1]) - A[0][1]*A[1][0] + A[0][0]*A[1][1])/det;

	  //det = (lambda[n] - A[1][1])*(lambda[n] - A[2][2]) - A[2][3]*A[2][3];
	  det = lambda[n] * lambda[n] - lambda[n] * A[1][1] - A[1][2] * A[1][2] - lambda[n] * A[2][2] + A[1][1] * A[2][2];
	  y = (lambda[n] * A[1][0] - A[2][2] * A[1][0] + A[1][2] * A[2][0]) / det;
	  z = (A[2][1] * A[1][0] + lambda[n] * A[2][0] - A[1][1] * A[2][0]) / det;

      normFactor = sqrt(1 + y*y + z*z);
      eigenvectors[0][n] = 1.0/normFactor;
      eigenvectors[1][n] = y/normFactor;
      eigenvectors[2][n] = z/normFactor;

	  //eigenvectors[n][0] = 1.0 / normFactor;
	  //eigenvectors[n][1] = y / normFactor;
	  //eigenvectors[n][2] = z / normFactor;

    }
	*/

    //diagonalize the quadrupole tensor
    //double *quadEigenValues = getEigenvalues(inertia);
    //quadEigenVectors = getEigenvectors(inertia, quadEigenValues);

}






















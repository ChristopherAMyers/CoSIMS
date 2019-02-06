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
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <queue>
#include <algorithm>
#include <locale>
#include <cctype>
#include "vector3D.h"
#include "myStructs.h"
using namespace std;

class molecule2
{
	ifstream file;
	int numH, numC, numO, numN, numP;// Added 06/28/15
        struct forceField
        {
          double mass;
          double epsilon;
          double sigma;
          string name;
          string symbol;
          int atomCount;
          int intMass;
        };

        initParams params;
	vector <forceField> fieldData; //stores forceField data from file

	//stores multipole and inertia options
	bool multipole;
	bool diagQuad, zeroDipole;
	double totalCharge;
	double dipole[3];
	double quadpole[3][3];
	double octupole[3][3][3];
	double inertia[3][3];
	double origInertia[3][3];
	double rotMat[3][3];
	double inertiaVectors[3][3];
	double inertiaValues[3];
	double quadVectors[3][3];
	double quadValues[3];

	const double avo_ext = 6.02214199e+23; //Avogadro's number'
	const double calorie_ext = 4.184; //number of jouels in 1 calorie
	const double protonMassKg_ext = 1.67262158e-27; // mss of a proton in kg
	const double heliumMassU_ext = 4.002602; //mass of helium in a.u.
	
	void assignDefaultFF();
	bool isInt(string);
	bool isDecimal(string);


public:
	ofstream xyzFile;
	dataList *data;

	molecule2();
	~molecule2(void);
	molecule2(const molecule2& mol_in) : data(new dataList(*mol_in.data) ) {}

	int numAtoms;
	void setCenter();
	void findCenter();
	void rotate(double, double, double);
	void importFile();
	double totalMass;
	double centerX, centerY, centerZ;
	double maxX, maxY, maxZ;
	void makeXYZ(double, double, double, double);
	void importTXT(); //Added 06/28/15
        void importPDB(); //Added 06/28/15
        void importPQR();
        void importMFJ(bool, bool);
        void setMolProperties(queue<double>, queue<double>, queue<double>, queue<double>, queue<string>, queue<string>);
        void printMolStats(); // Added 06/28/15

        //multipole annd inertia functions
            bool charge;
        void setMultipole();
        void setDipole();
        void setQuadpole();
        void setOctupole();
        void setInertia(bool);
        void setInertia(double[3][3]);
        void getEigenvalues(double[3][3], double(&)[3]);
        void getEigenvectors(double[3][3], double[3], double(&)[3][3]);
        void rotateLongVec();
        double delta(int, int);


};


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

#ifndef MYSTRUCTS_H
#define MYSTRUCTS_H
#include "vector3D.h"
#include <string>
#include <fstream>

//external (global) parameters
extern ofstream fileStream_ext;
extern string omegaFilePath_ext;
extern string fileLocaiton_ext; //TYPO
extern string molFilePath_ext;
extern string molFileName_ext;
extern string atomFilePath_ext;
extern string dataFilePath_ext;
extern string omegaFilePath_ext;
extern string xyzFilePath_ext;
extern string name_ext;
extern string subDir_ext;
extern string printMolName_ext;
extern double temp_ext;
extern double dt_ext;
extern double mass_ext;
extern double mu;
extern int trajNum_ext;
extern int iterations_ext;
extern int numThreads_ext;
extern long long int seed_ext;
extern int printRate_ext;
extern int maxClusterSize_ext;
extern int multipole_order_ext;
extern double cutOffRadius_ext; //remove
extern double tollerance_ext; //remove
extern double gMax_ext;
extern double gMin_ext;
extern double totalCharge_ext;
extern bool displayProg_ext;

extern double dispersion_radius_ext;
extern double multipole_radius_ext;

extern bool multipole_ext;
extern bool dispersionCutoff_ext;

extern bool charge_ext;
//extern bool orbitPi_ext;
extern bool projection_ext; //ellipsoid projection on/off
extern bool shorten_ext; //shorten trajectories
extern bool printTraj_ext; //print xyz files of trajectories
extern bool printTrajOnFail_ext; //print trajectory only if they fail = true, otherwise it prints all successful ones
extern bool printData_ext; //print starting data
extern bool printDataOnFail_ext; //print data only if they fail = true, otherwise it prints all successful ones
extern bool printMol_ext;
//extern double fixedVelocity_ext; //debugging only
extern bool printToColsole_ext;






//starting values for gas atom trajectories
//used for computing the cross section integral
typedef struct startVals
{
	double phi; //euler angle 1
	double theta; //euler angle2
	double gamma; //euler angle 3
	double g; //initial incoming velocity
	double b; //impact parameter
	vector3D <double>startPos; //starting position of gas atom
	vector3D <double>startVel;

	//these two vectors are just used to find
	//the angle between the incoming gas atom
	//and the deflected gas atom. startVec is computed
	//from the first two time steps of the trajectory,
	//while endVec is from the last two time steps. The
	//angle between these two vectors is the scattering angle "chi"
	vector3D <double>startVec;
	vector3D <double>endVec;
} startVals;

//holds min and max values for the integral;
typedef struct integralStats
{
	double phiMin;
	double phiMax;

	double thetaMin;
	double thetaMax;

	double gammaMin;
	double gammaMax;

	double bMin;
	double bMax;

	double gMin;
	double gMax;

	int numIntPts;

} integralStats;

typedef struct initParams
{
  string fileLocaiton; //TYPO
  string molFilePath;
  string molFileName;
  string atomFilePath;
  string dataFilePath;
  string omegaFilePath;
  string xyzFilePath;
  string iterFilePath;
  string name;
  string subDir;
  double temp;
  double dt;
  double mass;
  int trajNum;
  int iterations;
  int numThreads;
  int seed;
  int printRate;
  double resolution;
  double cutOffRadius;
  double tollerance;
  double gMax;
  double gMin;
  double totalCharge;

  bool multipole;
  bool charge;
  bool orbitPi;
  bool sampleDistr; //sample my MaxBoltz
  bool projection; //ellipsoid projection on/off
  bool shorten; //shorten trajectories
  bool printTraj; //print xyz files of trajectories
  bool printTrajOnFail; //print trajectory only if they fail = true, otherwise it prints all successful ones
  bool printData; //print starting data
  bool printDataOnFail; //print data only if they fail = true, otherwise it prints all successful ones



  /*
  initParams(const initParams &p)
  {
    fileLocaiton = p.fileLocaiton; filePath = p.filePath; molFileName = p.molFileName;



    printDataOnFail = p.printDataOnFail;
    iterations = p.iterations;
  }*/



}initParams;

typedef struct threadSharedParams
{
  double total;
  int trajCount;
  int iterCount;
  double tempOmega;
  double tempError;
}threadSharedParams;

struct localThreadParams
{
	double total;
	double total2;

	double prevTotal, prevTotal2;

	double currOmega;
	double currError;

	double g, b, phi, theta, gamma;

	vector3D<double> initPos;
	vector3D<double> initVel;

	int trajCount, prevTrajCount, missCount, orbitCount, hitCount
	, lowChiCount, stepCount, ranNumCount;

	double maxErgyPctDev;


};

struct dataList
       {
       public:
               vector3D <double>origCoords; //these are the inputed coords, only change by centering about CoM.
               vector3D <double>coords; //these are the ones the program will work with.
               double mass;
               double charge;
               double sigma;
               double epsilon;
               double epsAvg;
               double sigAvg;
               double sigAvg6;
               double sigAvg12;
               double eps24_mass;
               double eps4;
               double eps24;
               std::string atom;
               //struct dataList *next;
       };

struct userOptions
{
	bool name, dir, dt, temp, traj, iter, threads, seed, print_mol;
	bool print_traj, print_traj_fail, print_data, print_data_fail;
	bool dispersion_cutoff, dispersion_radius, multipole;
	bool multipole_radius, multipole_order, print_rate;
	bool max_cluster_size, charge, gMin, gMax, proj;
};




#endif

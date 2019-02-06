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

#ifndef MONTETRAJ_H_
#define MONTETRAJ_H_

#define _USE_MATH_DEFINES

#ifdef _WIN32
	#include <direct.h>
	#define GetCurrentDir _getcwd
	#define osParam 0
#elif defined _WIN64
	#include <direct.h>
	#define GetCurrentDir _getcwd
	#define osParam 0	
#else
	#include <unistd.h>
	#define GetCurrentDir getcwd
	#define osParam 1
#endif

#include <random>
#include <omp.h>
#include "molecule2.h"
#include "vector3D.h"
#include "diffEqSolver.h"
//#include "crosInt.h"
#include "maxBoltz.h"
#include "impactParam.h"
#include "VelocityDistribution.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <queue>
#include <sstream>



 //boltzmans constant in units of (A^2)(u)(ps^-2)(K^-1)
#define numGas 1 // number of gas atoms to simulate at one time


class monteTraj
{
  float maxDist; //maximum distance a gas atom can be away form molecule
  diffEqSolver gasAtom; //solves for the trajectories of the gas atoms
  molecule2 *molPtr; //points to the molecule data imported from main()
  ofstream xyzFile; // prints the trajectory in .xyz format, for DEBUGGING ONLY
  ofstream outMolFile; //prints the rotated molecule file, for DEBUGGING ONLY
  ofstream data; //prints trajectory stats, for DEBUGGING ONLY
  //ofstream omegaData; //prints the current value for the crossection, for DEBUGGING ONLY
  ofstream iterData; //prints a running value of omega DEBUG ONLY
  ofstream energyFile; //energy analysis, debugging only

  double *omega; //holds total total crossSection of molecule
  double *omegaErr; //holds the error for each cross-section
  double weightedMean;
  double weightedErr;
  double maxErgyPctDev;
  //double dOmega;
  maxBoltz<double> maxB; //maxwell-boltzmann distribution
  impactParam impactB; //class object to determine impact parameter

  int missCount; //number of missed trajectories
  int iterations; //number of sets of trajectory calculations
  double omegaAvg;
  double omegaStdDev;

 
  //std::thread threads;
  initParams params;
  threadSharedParams *shared;

  double const1;
  double const2;

  //printingVariables
  int currTrajCount, runningNumThreads, currSuccTraj;
  double currTotal, currTotal2;
  int prevCountMod, updateCond;
  

  struct threadParams
  {
    double bMax;
    diffEqSolver* gasAtom;
  };

  queue<string> xyzQueue;
  ostringstream xyzStream;
  //debugging purposes
  bool disp;
  double bCount;


  int distrPower;
  VelocityDistribution velDistr;
  double maxVelDistr;


  double areaRatio;
  int randNumSkipSize;
  unsigned int *restartGen;


  //double tempOmega, tempError, total;

  void updateCounts(localThreadParams&);



public:
  monteTraj(molecule2&, threadSharedParams&);
  virtual
  ~monteTraj();

  double gMean; //mean of maxBoltz distro
  double gMode; //mode of maxBoltz distro
  double gSigma; //stdDev of maxBOltz distro
  double gMin; //max velocity of gas atoms
  double gMax; //min velocity of gas atoms
  double temp; //temperature of the system
  double mu; //reduced mass of the gas atoms
  double heliumMass; // mass of helium atom
  double beta; // constant
  double  bMax; // maxium value of the impact parameter
  void findMaxDist();
  //void setUpTraj(int, int&);

  bool setUpTraj(int , int & , vector3D<double>& , vector3D<double>& 
	  , localThreadParams &, mt19937_64 &, uniform_real_distribution<double>&);
  void runTraj(); //compute trajectories of gas atoms
  void runThreads(int, bool, double, localThreadParams&);
  //static void *runThreads(void *args);
  //bool *threadComplete;//keep track which threads are still running
  //static void *callFunction(void *args);
  //debugging purposes only
  void printTraj(bool); //print the particle trajectory to XYX file
  void printMolecule();
  void storeTraj();
  void deleteStuff(); //for debugging purposes only
  void printData(int, double, localThreadParams&);
  void dispResults(int, int, double);
  void dispResults(localThreadParams&);

};

#endif /* MONTETRAJ_H_ */

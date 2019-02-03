/*
 * diffEqSolver.h
 *
 *  Created on: Mar 13, 2015
 *      Author: cmyers
 */

#ifndef DIFFEQSOLVER_H_
#define DIFFEQSOLVER_H_
#include <cmath>
#include <queue>
#include "vector3D.h"
#include "molecule2.h"
#include "forceAccel.h"
#include "myStructs.h"
#include <string>
#include <queue>
#include <fstream>
#include <sstream>

using namespace std;

class diffEqSolver {

	ofstream trajFile; //prints trajectory coordinates to file
	ofstream impact; ///DEBUG ONLY
	double dt; //sets integration step
	vector3D<double> position; //current position
	vector3D<double> velocity; //current velocity
	vector3D<double> accel; //current acceleration
	vector3D<double> oldPosition;
	vector3D<double> oldVelocity;

	//debug feature only
	ostringstream energyStream;
	bool printEnergy;

	int stepCount;
	double rDotOld, rDotNew;

	//for orbiting
	int maxOrbits;
	int orbits;

public:

	forceAccel forceFunct; //function to integrate
	diffEqSolver(); //constructor
	~diffEqSolver(); //destructor
	void setMol(molecule2&); //import molecular file
	void setPosition(vector3D<double>); //set initial velocity
	void setVelocity(vector3D<double>); //set initial position

	void setParams(initParams);
	vector3D<double> getPosition(); //retrieve current position
	vector3D<double> getVelocity(); //retrieve current velocity

	//verlet-velocity integrator
	//returns true if orbiting has not occurred
	bool verletVel(double, vector3D<double>&, vector3D<double>&, double, double, double, double&, localThreadParams&); //runge-kutta 4th order diffEq SOlver

	//Determines if orbiting has occurred by tracking the
	//projeciton of the position vector to the velocity vector
	bool trackRadius();

	//debugging features only
	void energySolver();
	void energyToFile(bool);

};

#endif /* DIFFEQSOLVER_H_ */

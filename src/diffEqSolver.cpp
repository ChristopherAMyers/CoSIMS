/*
 * diffEqSolver.cpp
 *
 *  Created on: Mar 13, 2015
 *      Author: cmyers
 */

#include "diffEqSolver.h"
#include <cmath>
#include <queue>
#include "vector3D.h"
#include "molecule2.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <queue>
#include <fstream>
#include <algorithm>

diffEqSolver::diffEqSolver() {

	//default integrator time-step
	dt = 0.01;

	stepCount = 0;
	orbits = 0;
	maxOrbits = 8;

    rDotNew = 0.0;
    rDotOld = 0.0;

    //prints kinetic and potential energy to test file
    //the writer is kept for debugging purposes, but it is
    //not called anywhere in the code. Reader must add own caller
    //and change 'printEnergy' to true.
    printEnergy = false;
    if(printEnergy)
    {
		impact.open("out/file/locaiton.txt", std::ios_base::out);
		impact <<"kin\tpoten\tsteps" << endl;
    }
}

diffEqSolver::~diffEqSolver() {

}

void diffEqSolver::setMol(molecule2 &mol)
{
	forceFunct.setMol(mol);
	forceFunct.setClusters();
}

void diffEqSolver::setPosition(vector3D <double>pos)
{position = pos;
stepCount = 0;}

void diffEqSolver::setVelocity(vector3D <double>vel)
{velocity = vel;}

vector3D<double> diffEqSolver::getPosition()
{return position;}

vector3D<double> diffEqSolver::getVelocity()
{return velocity;}


bool diffEqSolver::verletVel(double time, vector3D<double> &initPos, vector3D<double> &initVel
    , double a, double b, double c, double &chi, localThreadParams &p)
{

	/* Notes:
	 * It was found that orbiting of the gas particle around the ion
	 * could sometimes occur. This was inevitably due to energy not
	 * being conserved from a time-step that wasn't small enough.
	 *
	 * This was fixed by re-running the trajectory with a smaller
	 * time step. The orbiting logic originally used in this code
	 * is included for references only. This will be removed with
	 * a future update.
	 */


	  vector3D<double> tempVec;
	  vector3D<double> pos, vel, acc;
      double maxRad = max(a, max(b, c));
      double maxTries = 8, trajTries = 0;;
      double dt = time;
      double rDot_new, rDot_old, product, orbits;
      double startEnergy, endEnergy, energyPctDiff;
      int stepCount = 0;
      bool track;
      stepCount = 0;

	  startEnergy = forceFunct.getPotential(initPos) + 0.5*mass_ext*initVel.mag2();
      while(trajTries < maxTries && maxTries < 10)
      {
    	  //tasks to be performed prior to starting trajectory
          if(stepCount == 0)
          {
        	//set initial acceleration
            acc = forceFunct.getReducedVec(initPos);
            pos = initPos;
            vel = initVel;

            //initialize orbiting parameters
            rDot_old = 0;
            rDot_new = 0;
            if(printEnergy)
            	energyStream.str("");
			orbits = 0;

			//sets first orbit check
			rDot_new = pos.dot(vel);
          }


          //update position, velocities, and step counts
          pos = pos + vel*dt + acc*dt*dt*0.5;
          vel = vel + acc*dt*0.5;
          acc = forceFunct.getReducedVec(pos);
          vel = vel + acc*dt*0.5;
          stepCount++;


          //update orbiting
          track = true;
          rDot_old = rDot_new;
          rDot_new = pos.dot(vel);
          product = rDot_new*rDot_old;
          if(stepCount > 2)
          {
            if(product <0)
            {
              if(rDot_new > rDot_old)
                track =  true;
              else
              {
                orbits++;

                //purposely set to never happen
                if(orbits > maxOrbits)
                  track = false;
              }
            }
          }


          //this will never happen, logic is
          //kep for reference only
          if(!track )
          {
            trajTries++;
            dt = time/((trajTries+1));
            cout << "TRACKING\n";
          }
          else
          {

        	//multiply vector by elipsoid axis lengths
            tempVec.setEqual(pos.X/a, pos.Y/b, pos.Z/c );

            //check if (1) particle is outside of the ellipsoid
            //or, if a sphere is used, (2) if the particle is outside of the
            //sphere radius, still defined by the ellipsoid axes.
            if( ( (tempVec.mag2() > 1 && projection_ext) || (tempVec.mag2() > maxRad && (!projection_ext) ) ) && stepCount > 10)
            //if(tempVec.mag2() > 1 && stepCount > 10)
			{

            	//calculate total energy of gas particle
			    endEnergy = forceFunct.getPotential(pos) + 0.5*mass_ext*vel.mag2();
			    //endEnergy = startEnergy;

			    //calcualte percent difference form initial energy
			    energyPctDiff = 100 * abs(startEnergy - endEnergy) / startEnergy;
			    chi = initVel.Angle(vel);
			    p.stepCount = stepCount;

			    //if energy difference is greater than 0.5%,
			    //then re do the trajectory with a smaller time step
			    if(energyPctDiff > 0.5)
			    {
					   trajTries++;
					   maxTries++;
					   dt = time / ((trajTries + 1));
					   stepCount = 0;
			    }

			    //else, exit diff-eq solver indicating a successful trajectory
			    else
			    {
					   if (energyPctDiff > p.maxErgyPctDev)
							   p.maxErgyPctDev = energyPctDiff;
					   //printf("Energy: %6.3f\n", energyPctDiff);
					   return true;
			    }
			}
          }


      }

      chi = M_PI_2;
      p.stepCount = stepCount;
      return false;
}

void diffEqSolver::energySolver()
{
  double kinetic = 0.5*4*velocity.dot(velocity);
  double potential = forceFunct.getPotential(position);

  energyStream << setprecision(4);
  energyStream << kinetic << "\t" << potential << "\t"  << kinetic + potential << "\t"
      << position.X << "\t" << position.Y << "\t" << position.Z << "\t"
      << velocity.X << "\t" << velocity.Y << "\t" << velocity.Z <<
      "\n";
}

void diffEqSolver::energyToFile(bool print)
{
  if(print)
    impact << energyStream.str();
}


//this funciton is kept for reference only
bool diffEqSolver::trackRadius()
{
  double product;
  rDotOld = rDotNew;
  rDotNew = position.dot(velocity);
  product = rDotNew*rDotOld;
  if(stepCount > 2)
    if(product <0)
      if(rDotNew > rDotOld)
        return true;
              else
              {
                      orbits++;
                      if (orbits > maxOrbits)
                              return false;
                      else
                              return true;
              }
    else
      return true;
  else return true;


/*Logic Notes:
 * Check if rDot swaps from incoming (positive valued)
 * to outgoing (negative valued) by a sign swap.
 * If there is not a change in direction, then product
 * is positive, else it gives a negative number.
 * Then check if is scattering (rDotNew > rDotOld)
 * or bound to the molecule (rDotNew < rDotOld).
 */


	
}

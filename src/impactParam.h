/*
 * impactParam.h
 *
 *  Created on: Jul 14, 2015
 *      Author: cmyers
 */

#ifndef IMPACTPARAM_H_
#define IMPACTPARAM_H_
#define _USE_MATH_DEFINES
#include <vector>
#include <cstdlib>
#include <cfloat>
#include <fstream>
#include "molecule2.h"
#include "forceAccel.h"
#include "vector3D.h"

class impactParam
{

  molecule2 *mol; // points the molecular data for the atom coordinates
  forceAccel forceFunct;
  double *phi;
  double *theta;
  int numTestPoints;

  double aUncert, bUncert, cUncert, rmsd;
  double cov_aa, cov_bb, cov_cc, cov_ab, cov_ac, cov_bc;
  double avgAxes, avgAxesUncert;
  double surfArea, surfAreaUncert;
  double minPhi, minCosT, minSinT;
  double aMin, bMin, cMin;
  double aUncertMin, bUncertMin, cUncertMin;
  double enlargeAmount;

  vector<vector3D<double> > boundryPoints; //stores boundry points
  //vector3D<double> *boundryPoints;
  //int boundryPtSize;
  

  //ofstream omegaFile;
public:
  impactParam();
  virtual
  ~impactParam();
  double a, b, c; // axes for ellipsoid
  double aPlane, bPlane; // params for plane
  double phiPlane, cosThetaPlane; //angles for rotaiton on plane
  void setMol(molecule2&); //get molecular information
  void rotateMol(); // finds the atom with longest vector and points
                   // the z-axes along that vector
  vector3D<double> inverseM(vector3D<double>&, double, double);
  void findEllipsoid(); //finds the parameters 'a', 'b', and 'c'
  void findPlane();
  double getEllipMag(double, double, double);
  double getEllipMag2(double, double, double);
  double getMax();
  double getMaxAtomDist();
  double derivF(double);
  void refineEllip();
  void expandEllip();
  void getTestPoints();
  void rotateMol(double, double, double);
  void rotateMolInv(double, double, double);
  void rotateMolInvEuler(double, double, double, double);
  vector<vector3D<double> > genRandPoints(int);
  void multEllip(int);
  void printResults();

  void setBoundry(double);
  void setBoundry_old(double);
  void findEllipsoidBound();
  void printBoundryPoints(string); //Debuging-only
  void enlargeEllipBoundry();
  void printBoundryResults();
};

#endif /* IMPACTPARAM_H_ */

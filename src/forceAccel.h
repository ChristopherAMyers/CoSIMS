/*
 * forceAccel.h
 *
 *  Created on: Aug 17, 2017
 *      Author: cmyers
 */

#ifndef FORCEACCEL_H_
#define FORCEACCEL_H_

#include "vector3D.h"
#include "molecule2.h"
#include "PCM_Cluster.h"
#include <string>

namespace std
{

  class forceAccel
  {
    molecule2 *mol;
    double cutOffRadius;
    double resolution;
    double mass;
    double heliumMass;
    double dipoleScalar;


	double delta(int, int);

	PCM_Cluster cluster;
	vector<vector<int> > clusterGroups;
	vector<vector3D<double> > clusterPos;
	vector<double> clusterRadius;
	vector<bool> centerReplace;
	int totalCenterReplaced;
	double avgClusterSize;
	double avgClusterRad;

	vector<double> totalCharge;
	vector<vector<double> > dipole;
	vector<vector<vector<double> > > quadPole;

	void addMultipole(vector<int>, vector3D<double>);
	void addClusterPos(vector3D<double>, int);
	void addClusterRadius(double);
	void printClusters();

	vector3D<double> multipoleCutoff(vector3D<double>&);
	vector3D<double> multipoleOnly(vector3D<double>&);
	vector3D<double> dispersionOnly(vector3D<double>&);
	vector3D<double> chargeDipoleLJ(vector3D<double>&);
	vector3D<double> dispersionNoCharge(vector3D<double>&);
	vector3D<double> dispersionFull(vector3D<double>&);


	int forceType;


  public:
    forceAccel();
    virtual
    ~forceAccel();

    void setMol(molecule2&);
	void setClusters();
    vector3D<double> getReducedVec(vector3D<double>&);
    vector3D<double> getVec(vector3D<double>&);
    double getSecDerivMag(vector3D<double>&);
    double getPotential(vector3D<double>&);



  };

} /* namespace std */

#endif /* FORCEACCEL_H_ */

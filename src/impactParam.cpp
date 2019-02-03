/*
 * impactParam.cpp
 *
 *  Created on: Jul 14, 2015
 *      Author: cmyers
 */

#include "impactParam.h"
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

impactParam::impactParam()
{
  a = 0; b = 0; c = 0;
  // TODO Auto-generated constructor stub

}

impactParam::~impactParam()
{
  // TODO Auto-generated destructor stub
}

void impactParam::setMol(molecule2 &inMol)
{
  mol = &inMol;
  forceFunct.setMol(inMol);
}


void impactParam::rotateMol()
{
  //max vector magnitude and it's repective atom number
  double maxR = 0;
  int atomNum = 0;

  //check which atom has the margest vector magnitude from the origin
  for(int i = 0; i < mol->numAtoms; i++)
    if(mol->data[i].coords.mag() >= maxR)
    {
      maxR = mol->data[i].coords.mag();
      atomNum = i;
    }

  //get that atoms phi and theta components
  //THESE FUNCTIONS MUST BE ACTILVE IN VECTOR3d HEADER FILE!
  /*
  double phi = mol->data[atomNum].coords.getPhi();
  double theta = mol->data[atomNum].coords.getTheta();




  //rotate the molecule so that the recently chosen atom
  //points along the z-axis
  //rotateMol(-phi, cos(-theta), sin(-theta));
  minPhi = -phi; minCosT = cos(-theta); minSinT = sin(-theta);
  for(int i = 0; i <= mol->numAtoms; i++)
    mol->data[i].coords = inverseM(mol->data[i].coords, phi, theta);
*/
}

vector3D<double> impactParam::inverseM(vector3D<double> &vec, double phi, double theta)
{
  vector3D<double> temp;

  temp.X = vec.Y*cos(phi) - vec.X*sin(phi);
  temp.Y = -vec.X*cos(phi)*cos(theta) - vec.Y*cos(theta)*sin(phi) + vec.Z*sin(theta);
  temp.Z = vec.X*cos(phi)*sin(theta) + vec.Y*sin(phi)*sin(theta) + vec.Z*cos(theta);

  return temp;
}

void impactParam::rotateMol( double phi, double cTheta, double sTheta)
{
  // precalculate trig function values
  double cPhi = cos(phi); double sPhi = sin(phi);
  //double cTheta = cos(theta); double sTheta = sin(theta);

  //pre-calculate multiplication of trig values
  double cPcT = cPhi*cTheta;
  double cPsT = cPhi*sTheta;
  double cTsP = cTheta*sPhi;
  double sTsP = sTheta*sPhi;

  // temp values of the vector to rotate
  double x, y, z;

  //rotate molecule component by component
  //for(int n = 0; n < (signed)vec.size(); n ++)
  for(int n = 0; n < mol->numAtoms; n ++)
  {
    x = mol->data[n].origCoords.X;
    y = mol->data[n].origCoords.Y;
    z = mol->data[n].origCoords.Z;

    mol->data[n].coords.X = x*cPcT + z*cPsT - y*sPhi;
    mol->data[n].coords.Y = x*cTsP + y*cPhi + z*sTsP;
    mol->data[n].coords.Z = z*cTheta - x*sTheta;
  }
}

void impactParam::rotateMolInv(double phi, double cTheta, double sTheta)
{
	// precalculate trig function values
	double cPhi = cos(phi); double sPhi = sin(phi);
	//double cTheta = cos(theta); double sTheta = sin(theta);

	//pre-calculate multiplication of trig values
	double cPcT = cPhi*cTheta;
	double cPsT = cPhi*sTheta;
	double cTsP = cTheta*sPhi;
	double sTsP = sTheta*sPhi;

	// temp values of the vector to rotate
	double x, y, z;

	//rotate molecule component by component
	//for(int n = 0; n < (signed)vec.size(); n ++)
	for (int n = 0; n < mol->numAtoms; n++)
	{
		x = mol->data[n].origCoords.X;
		y = mol->data[n].origCoords.Y;
		z = mol->data[n].origCoords.Z;

		mol->data[n].coords.X = x*cPcT + y*cTsP - z*sTheta;
		mol->data[n].coords.Y = y*cPhi - x*sPhi;
		mol->data[n].coords.Z = x*cPsT + y*sTsP + z*cTheta;
	}
}

/* This functions uses a least squares calculation to determine the
 * axes components of an ellipsoid. The mathematical function
 *        sum{i = 1, atomNum} (1 - A*X^2 - B*Y^2 - C*Z^2)^2
 * is the least squares function, where it's derivatives w.r.t a, b, and c
 * are set to zero, and solved for here. Greek variables represent summations,
 * while Latin variables represent parameters we are solving for.
 *
 * We solve the matrix equation
 *
 * [  mu   delta  sigma  ][ A ]   [ alpha ]
 * [ delta   nu    rho   ][ B ] = [  beta ]
 * [ sigma  rho  epsilon ][ C ]   [ gamma ]
 *
 * using Cramer's Rule so linear systems of equations to solve for
 * A, B, and C. The x, y, and z axes of the ellipsoid are given as
 * one over the square-root of A, B, and C, respectively.
 *
 *
 *
 *
 */

double impactParam::getMaxAtomDist()
{
	//max vector magnitude and it's repective atom number
	double maxR = 0;
	int atomNum = 0;

	//check which atom has the margest vector magnitude from the origin
	for (int i = 0; i < mol->numAtoms; i++)
		if (mol->data[i].coords.mag() >= maxR)
		{
			maxR = mol->data[i].coords.mag();
			atomNum = i;
		}

	return maxR;
}

void impactParam::findEllipsoid()
{

  //stores summation information
  double mu = 0, nu = 0, epsilon = 0, delta = 0, sigma = 0, rho = 0;
  double alpha = 0, beta = 0, gamma = 0;
  double A = 0, B = 0, C = 0;
  rmsd = 0;

  //stores uncertainty information
  aUncert = 0, bUncert = 0, cUncert = 0;

  //double x2 = 0; double y2 = 0; double z2 = 0;
  //double x4 = 0; double y4 = 0; double z4 = 0;

  double x2 = 0, y2 = 0, z2 = 0, x4 = 0, y4 = 0,
      z4 = 0, x2y2 = 0, x2z2 = 0, y2z2 = 0;

  //calculate sumamtion over data coordinates
  for (int i = 0; i < mol->numAtoms; i ++)
  {

    x2 = mol->data[i].coords.X*mol->data[i].coords.X;
    y2 = mol->data[i].coords.Y*mol->data[i].coords.Y;
    z2 = mol->data[i].coords.Z*mol->data[i].coords.Z;



    x4 = x2*x2; y4 = y2*y2; z4 = z2*z2;
    x2y2 = x2*y2; x2z2 = x2*z2; y2z2 = y2*z2;

    mu = mu + x4;
    nu = nu + y4;
    epsilon = epsilon + z4;

    delta = delta + x2*y2;
    sigma = sigma + x2*z2;
    rho = rho + y2*z2;

    alpha = alpha + x2;
    beta = beta + y2;
    gamma = gamma + z2;

  }


  double det = epsilon*mu*nu - epsilon*delta*delta - mu*rho*rho + 2*delta*rho*sigma - nu*sigma*sigma;

  A = (alpha*(nu*epsilon - rho*rho)
      + delta*(rho*gamma - beta*epsilon)
      + sigma*(beta*rho - nu*gamma))/det;

  B = (mu*(beta*epsilon - gamma*rho)
      + alpha*(rho*sigma - delta*epsilon)
      + sigma*(delta*gamma - beta*sigma))/det;

  C = (mu*(nu*gamma - beta*rho)
      + delta*(beta*sigma - delta*gamma)
      + alpha*(delta*rho - nu*sigma))/det;


  a = 1/sqrt(A); b = 1/sqrt(B); c = 1/sqrt(C);

  rmsd = sqrt((A*A*mu + B*B*nu + C*C*epsilon +
      2*(A*B*delta + B*C*rho + A*C*sigma - A*alpha - B*beta - C*gamma ) + mol->numAtoms)/(mol->numAtoms - 1));


  aUncert = (a/A)*sqrt((epsilon*nu - rho*rho)*2/det)*rmsd;
  bUncert = (b/B)*sqrt((epsilon*mu - sigma*sigma)*2/det)*rmsd;
  cUncert = (c/C)*sqrt((mu*nu - delta*delta)*2/det)*rmsd;

}

void impactParam::printResults()
{
  cout << endl;
  cout << setprecision(4);
  cout << "Impact Parameter properties:" << endl;
  cout << "RMSD: " << rmsd << endl;
  cout << "Phi rotation: " << phiPlane << " Theta rotation: " << acos(cosThetaPlane) << endl;
  cout << "\t X \t Y \t Z " << endl;
  cout << "_______________________________" << endl;
  cout << "Axes \t " << a << "\t" << b << "\t" << c << endl;
  cout << "Uncert \t" << aUncert << "\t" << bUncert << "\t" << cUncert << endl;
  cout << "Percnt \t" << aUncert*100/a << "\t" << bUncert*100/b << "\t" << cUncert*100/c << endl;
}

void impactParam::findPlane()
{	double nu = 0, mu = 0, delta = 0, alpha = 0, beta = 0;
	double x2 = 0, y2 = 0, x = 0, y = 0, z = 0;
	double denom = 0;

	for (int i = 0; i < mol->numAtoms; i++)
	{
		x = mol->data[i].coords.X;
		y = mol->data[i].coords.Y;
		z = mol->data[i].coords.Z;
		x2 = x*x;
		y2 = y*y;

		nu = nu + x2;
		mu = mu + y2;
		delta = delta + x*y;
		alpha = alpha + x*z;
		beta = beta + y*z;
	}

	denom = nu*mu - delta*delta;
	aPlane = (mu*alpha - delta*beta) / denom;
	bPlane = (nu * beta - delta*alpha) / denom;

	phiPlane = atan2(bPlane, aPlane);
	cosThetaPlane = -1 / sqrt(aPlane*aPlane + bPlane*bPlane + 1);
}

vector< vector3D<double> > impactParam::genRandPoints(int numPoints)
{
	vector<vector3D<double> > pointList;
	pointList.reserve(numPoints);
	double bMax = max(max(a, b), c);
	double sinT = 0, cosT = 0, cosG = 0, sinG = 0, cosP = 0, sinP = 0;
	bool acceptN, acceptM;
	double dot = 0, discrim = 0, vMag2 = 0, posMag2 = 0, t1 = 0, t2 = 0;

	vector3D<double> tempPos, tempPos1, tempVel, tempVel1;

	double bCheck = 0, bb = 0, phi = 0, gamma = 0;

	for (int n = 0; n < numPoints; n += 2)
	{

		acceptN = false;
		while (!acceptN)
		{
		        acceptM = false;
			while (!acceptM)
			{
				bb = ((double)rand() / (double)RAND_MAX)*bMax;
				bCheck = ((double)rand() / (double)RAND_MAX) * 2 / bMax;
				phi = ((double)rand() / (double)RAND_MAX) * 2 * M_PI;
				cosT = ((double)rand() / (double)RAND_MAX)*2 - 1;
				gamma = ((double)rand() / (double)RAND_MAX) * 2 * M_PI;
				if (bCheck < (2 * bb / (bMax*bMax))) acceptM = true;
			}
			sinT = sqrt(1 - cosT*cosT);
			cosG = cos(gamma); sinG = sin(gamma);
			cosP = cos(phi); sinP = sin(phi);
			tempPos.setEqual(0, 0, 0);
			tempPos.X = bb*(cosG*cosP*cosT - sinG*sinP) + bMax*cosP*sinT;
			tempPos.Y = bb*(cosG*sinP*cosT + sinG*cosP) + bMax*sinP*sinT;
			tempPos.Z = bMax*cosT - bb*cosG*sinT;
			tempPos1.setEqual(tempPos.X / a, tempPos.Y / b, tempPos.Z / c);

			tempVel.setEqual(-cosP*sinT, -sinP*sinT, -cosT);
			tempVel1.setEqual(tempVel.X / a, tempVel.Y / b, tempVel.Z / c);

			dot = tempPos1.dot(tempVel1);
			vMag2 = tempVel1.mag2();
			posMag2 = tempPos1.mag2();
			discrim = dot*dot - vMag2*(posMag2 - 1);
			if(discrim > 0) acceptN = true;
		}

		t1 = (-dot - sqrt(discrim))/vMag2;
		t2 = (-dot + sqrt(discrim))/vMag2;

		pointList.push_back(tempPos + tempVel*t1);
		pointList.push_back(tempPos + tempVel*t2);

		vector3D<double> temp = (tempPos + tempVel*t1);
		temp.setEqual(temp.X/a, temp.Y/b, temp.Z/c);

	}
	
	return pointList;
}


void impactParam::multEllip(int num)
{
	cout << "NUMBER OF ATOMS: " << mol->numAtoms << endl;

	
	if (mol->numAtoms == 1)
	{
		cout << "Only one atom in molecule is present. Using a sphere instead\n";
		a = mol->data[0].sigAvg * 2.5;
		b = a; c = a;
	}
	else if (mol->numAtoms < 15)
	{
		a = 2 * getMaxAtomDist();
		b = a; c = a;
		cout << "maxDistance for 3 atoms: " << a << endl;
	}
	else
	{
		//plane method
		//rotateMol();
		//findPlane();
		//rotateMol(phiPlane, cosThetaPlane, sqrt(1 - cosThetaPlane*cosThetaPlane));
		//findEllipsoid();
	//}
		//printResults();

	//return;


		rotateMol();
		findEllipsoid();
		printResults();
		double randNum = 0.0, cosT = 1.0, sinT = 0.0; // random phi, cos(theta), sin(theta)
		double leastSquare = 0;
		double minLeastS = rmsd; // minimum RMSD from monte carlo search
		bool success; // tracks if successfully found a good ellipsoid

		//if a valid ellipsoid with percent error less than 50% is found
		if ((aUncert * 100 / a < 50) && (bUncert * 100 / b < 50) && (cUncert * 100 / c < 50))
			success = true;
		else
			success = false;

		cout << "SUCCESS: " << success << endl;
		cout << "SUCCESS: " << (aUncert * 100 / a < 50) << endl;
		cout << "SUCCESS: " << (bUncert * 100 / b < 50) << endl;
		cout << "SUCCESS: " << (cUncert * 100 / c < 50) << endl;

		minCosT = 1.0; minSinT = 0.0; minPhi = 0;

		//double minPhi, minCosT, minSinT;
		cout << setprecision(10);
		for (int n = 0; n < num; n++)
		{
			//sample theta by sin(theta)
			randNum = ((double)rand() / (double)RAND_MAX);
			cosT = 2 * randNum - 1;
			sinT = sqrt(1 - cosT*cosT);

			//sample phi uniformly
			randNum = ((double)rand() / (double)RAND_MAX) * 2 * M_PI;
			rotateMol(randNum, cosT, sinT);
			findEllipsoid();

			// if the calculated RMSD with random orientation is better than the last AND percent errors
			// in each axis lengths are less than 50%
			if ((rmsd < minLeastS) && (aUncert * 100 / a < 50) && (bUncert * 100 / b < 50) && (cUncert * 100 / c < 50))
			{
				if (!success) success = true;
				minLeastS = rmsd;
				minPhi = randNum;
				minCosT = cosT;
				minSinT = sinT;
				aMin = a; bMin = b; cMin = c;
				aUncertMin = aUncert; bUncertMin = bUncert; cUncertMin = cUncert;
			}
		}

		//is can't find an ellipsoid with appropriate percent error, use a sphere instead
		if (!success)
		{
			cout << "Warning: Failed to find a proper ellipsoid. \n"
				<< "Using a sphere instead. Number of atoms may be too small\n";
			a = mol->data[0].sigAvg * 2.5;
			b = a; c = a;
		}
		else
		{
			rotateMol(minPhi, minCosT, minSinT);
			findEllipsoid();
		}

		cout << "MIN:    " << minLeastS << endl;
		phiPlane = minPhi; cosThetaPlane = minCosT;
		printResults();
	}
	
}


double impactParam::getMax(){return max(a, max(b,c));}



double impactParam::getEllipMag(double phi, double theta, double gamma)
{
  return getEllipMag2(phi, theta, gamma);
  //for orthogonal vector to that of the ellipsoid
  vector3D<double> v(-b*c*sin(phi), a*c*cos(phi), 0.0);
  //form ellipsoid vector
  vector3D<double> ellipsoid(a*cos(phi)*sin(theta), b*sin(phi)*sin(theta), c*cos(theta));

  //form normal vector the the ellipsoid
  vector3D<double> k = ellipsoid/ellipsoid.mag();

  //use Rodrigue's rotation formula to rotate vector about gamma - Pi/2
  //v.rotateVec(k, gamma - M_PI_2);
  //v.rotateVec(k, gamma);
  v.rotateVec(k, gamma-M_PI_2);
  vector3D<double> rotation = v;

  //renormalize vector
  rotation.X = rotation.X/a; rotation.Y = rotation.Y/b; rotation.Z = rotation.Z/c;
  v = v/rotation.mag();
  return v.mag();

  //find new azimuthal and polar angles of this rotated vector
  double phiNew = atan2(rotation.Y, rotation.X);
  double thetaNew = acos(rotation.Z/rotation.mag());

  //form vector from these two angles that intersects the ellipsoids surface
  vector3D<double> newEllipsoid(a*cos(phiNew)*sin(thetaNew), b*sin(phiNew)*sin(thetaNew), c*cos(thetaNew));

  //return the magnitude of this vector
  return newEllipsoid.mag();

}

double impactParam::getEllipMag2(double phi, double theta, double gamma)
{
  //for orthogonal vector to that of the ellipsoid
  vector3D<double> v(-sin(phi), cos(phi), 0.0);

  //form normal vector
  vector3D<double> k(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));

  //use Rodrigue's rotation formula to rotate vector about gamma - Pi/2
  v.rotateVec(k, gamma-M_PI_2);
  vector3D<double> rotation = v;

  //renormalize vector
  rotation.X = rotation.X/a; rotation.Y = rotation.Y/b; rotation.Z = rotation.Z/c;
  v = v/rotation.mag();
  return v.mag();
}

void impactParam::refineEllip()
{
  cout << endl;
  double arg = 1.4;
  double tempArg = arg;
  double t = 0.01;
  double precision = 0.0000000001;
  double oldArg = arg + t;
  int count = 0;




  cout << setprecision(6);
  do
  {

    if(count % 10 == 0)
      cout << arg << "  " << "  " << oldArg << "  " << (oldArg - arg) << "  " << derivF(arg) << "  " << t << endl;
    oldArg = arg;

    tempArg = oldArg - t*derivF(arg);
    while(fabs(oldArg - tempArg) > 0.005)
    {t = t*0.8;
    cout << t << "\t"  << endl;
    tempArg = oldArg - t*derivF(arg);}


    arg = tempArg;
    count ++;

  }while(fabs(oldArg - arg) > precision);
  cout << "refine factor: " << arg << endl;
  a = a*arg;
  b = b*arg;
  c = c*arg;

}

void impactParam::expandEllip()
{

  forceFunct.setMol(*mol);
  double avgForce;
  vector3D<double> temp;
  bool accept = false;
  double aInc = a*0.01;
  double bInc = b*0.01;
  double cInc = c*0.01;

  do
  {
    avgForce = 0;

    //calculate the average force on the surface of the ellipsoid
    for(int i = 0; i < mol->numAtoms; i++)
    {
      temp.setEqual(mol->data[i].coords.X/a,mol->data[i].coords.Y/b,mol->data[i].coords.Z/c);
      temp = mol->data[i].coords/temp.mag();
      avgForce += forceFunct.getVec(temp).mag();
    }
    avgForce = avgForce/mol->numAtoms;

    //check if the average force is within the tollerance
    if(avgForce > 0.1)
    {
      cout << avgForce << "  " << a << endl;
      a = a + aInc;
      b = b + bInc;
      c = c + cInc;
      accept = false;
    }
    else accept = true;



  }while(!accept);
}

double impactParam::derivF(double x)
{
  getTestPoints();

  vector3D<double> u, v, v1,r;
  double force = 0;
  double forceDeriv;
  double epsilonAvg = 0;
  double sigmaAvg = 0;
  double total = 0;
  numTestPoints = mol->numAtoms;

  for(int j = 1; j < mol->numAtoms; j ++)
  {
    force = 0; forceDeriv = 0;
    v1.setEqual(a*cos(phi[j])*sin(theta[j]), b*sin(phi[j])*sin(theta[j]),c*cos(theta[j]));
    //cout << "here here" << endl;
    for(int i = 0; i < numTestPoints; i ++)
    {
      //cout << "here here" << endl;
      if(i != j)
      {
        u =mol->data[j].coords;
        v.setEqual(x*a*cos(phi[j])*sin(theta[j]), x*b*sin(phi[j])*sin(theta[j]),x*c*cos(theta[j]));
        r = u - v;

        epsilonAvg = sqrt(mol->data[i].epsilon*mol->data[j].epsilon);
        sigmaAvg = sqrt(mol->data[i].sigma*mol->data[j].sigma);

        force += 24*epsilonAvg*((pow(sigmaAvg,12)/pow(r.mag(), 13)) - (pow(sigmaAvg,6)/pow(r.mag(), 7)));
        forceDeriv += 24*epsilonAvg*(26*(pow(sigmaAvg,12)/pow(r.mag(), 14)) - 7*(pow(sigmaAvg,6)/pow(r.mag(), 8)))*r.dot(v1);

      }
    }
    //total += (0.01 + force)*forceDeriv;
    //total += (force/abs(force))*(0.01 + force)*forceDeriv;
   if(force > 0)
      total += (force - 0.01)*forceDeriv;
    else
      total += (0.01- force)*forceDeriv;


  }

return (total*2/mol->numAtoms);

}

void impactParam::getTestPoints()
{
  phi = new double[mol->numAtoms];
  theta = new double[mol->numAtoms];

  double x, y, z;

  for (int i = 0; i < mol->numAtoms; i ++)
  {
    x = mol->data[i].coords.X;
    y = mol->data[i].coords.Y;
    z = mol->data[i].coords.Z;
    phi[i] = atan2(a*y,b*x);
    theta[i] = atan2(c*y,b*z*sin(phi[i]));
  }
}

void impactParam::setBoundry(double minKin)
{
	/*
	This functions generates the potential energy boundry
	that will be fit to an ellipsoid. Any particle outside this
	boundry (and thus the ellipsoid) is defined to have a 
	potential energy deemed to be zero (by approximation).
	The zero potential energy is chosen such that the kinetic
	energy of the slowest possible moving gas particle is only
	lowered by 1%. The boundry thus defines the starting and 
	ending position	of all gas particle trajectories.
	*/

  vector<vector3D<double> > boundryPoints_temp;
    int numRotations = 20, count = 0;				//number of rotations about z-axis
    double theta;							//holds angle of rotation
    vector3D<double> vel, pos, prevPos, blank(0,0,0);		//vectors for determining test points
    double velMag = 0.1;					//incrimental vector 
    double posMag =
	    abs(max(abs(mol->maxX), abs(mol->maxY)));	//starting distance from z-axis
    double startX, startY;					//start positions
    double zPos = 0.0, zMax = 0.0;			//distance along z-axis, maximum extent along z-axis
    double energy;							//holds potential at pos
    //double minKinEnergy = 0.5*4*2*2;		//minimum kinetic energy value
    double slope = 1.0;						//determins search directions
    double zInc = 0.0; int maxZLevels = 25;		//z-axis incrimental values
    bool accept = false, hit = false;		//allows exit from while-loops

    //first, find maximum extent along z-axis
    accept = false;
    pos.setEqual(0, 0, abs(mol->maxZ) + 20);
    energy = forceFunct.getPotential(pos);
    minKin = 0.01*minKin;

	//printf("\ninitPos: %6.4f %6.4e %2f\n", pos.Z, energy + minKin, slope);
    if(energy + minKin < 0)
    {
      slope = -1;
    }
    while (!accept)
    {
      prevPos = pos;
      pos.Z -= 0.5*slope;
      energy = forceFunct.getPotential(pos);

	  //printf("initPos: %6.4f %6.4e %2f\n", pos.Z, energy + minKin, slope);

	  //found closest distance along z-axis
	  if ((energy + minKin)*slope < 0)
      //if(energy + minKin > 0)
        {
          zMax = max(pos.Z, prevPos.Z);
          accept = true;
        }
    }

	//start search from bottom-up
      zMax *=1.3;
	  posMag += abs(zMax - abs(mol->maxZ));
      maxZLevels = ceil(max(25, (int)ceil( zMax)));
      zPos = -zMax;
      zInc = 2*zMax/(maxZLevels);
      boundryPoints_temp.resize((int)(maxZLevels + 1)*numRotations, blank);
      //boundryPoints = new vector3D<double>[2*maxZLevels + 1];
      //while(zPos < zMax)
      omp_set_num_threads(numThreads_ext);
      #pragma omp parallel for reduction(+:count) private(zPos, slope, theta, startX, startY, pos, vel, hit, energy, prevPos)
      for(int i = 0; i < (maxZLevels + 1); i ++)
      {


        zPos = -zMax + i*zInc;
        //cin.get();
        //printf("Zpos: %6.4g %6.4f %4i %4.6f %4.6f %2f\n", zPos, zInc, i, zMax, abs(mol->maxZ), slope);
		//loop through each rotation about z-axis
        for(int n = 0; n < numRotations; n ++)
        {
          //create starting positions
          slope = 1;
          theta = n*2*M_PI/((double)numRotations);
          startX = posMag*cos(theta);
          startY = posMag*sin(theta);
          pos.setEqual(startX, startY, zPos);
          vel.setEqual(-velMag*cos(theta), -velMag*sin(theta), 0.0);
          hit = false;

          //get potential energy at pos
          energy = forceFunct.getPotential(pos);

          //if within the boundry that is being set up,
          //search from inside-out instead
          if(energy + minKin < 0)
            {
              slope = -1;
            }
          vel = vel*slope;

          while(!hit)
          {

            energy = forceFunct.getPotential(pos);

	    //potential energy vs kinetic energy condition
            if((energy + minKin)*slope <  0)
              {
                hit = true;
		//choose point with farthest distance
                if(prevPos.mag() > pos.mag())
                  boundryPoints_temp[i*numRotations + n] = prevPos;
                else
                  boundryPoints_temp[i*numRotations + n] = pos;
                count ++;
              }

            //if energy condition was not met,
            //but traversed across other side of molecule
            else if((startX == 0 && startY/pos.Y <= 0)
                || (startY == 0 && startX/pos.X <= 0)
                || (startX/pos.X <= 0 && startY/pos.Y <= 0))
            {
                hit = true;
                //boundryPoints_temp[i*numRotations + n] = pos;
                //boundryPoints[i*numRotations + n] = vector3D<double>(0,0,0);
				//cout << "Boundry: " << pos.X << "  " << pos.Y << "  " << pos.Z << endl;
            }
            //update position vectors
            else
            {
              prevPos = pos;
              pos = pos + vel;
            }
          }
        }

      }


      boundryPoints.reserve(count);
      count = 0;
      bool test;
      for(int i = 0; i < (signed)boundryPoints_temp.size(); i ++)
      {

        test = !(boundryPoints_temp[i].X == 0 && boundryPoints_temp[i].Y == 0 && boundryPoints_temp[i].Z == 0);
        if(!(boundryPoints_temp[i].X == 0 && boundryPoints_temp[i].Y == 0 && boundryPoints_temp[i].Z == 0))
          {
            boundryPoints.push_back(boundryPoints_temp[i]);
            count ++;
          }
      }





      //cout << endl;
      //for(int n = 0; n < boundryPoints.size(); n ++)
        //printf("Boundry: %8.4f  %8.4f  %8.4f\n", boundryPoints[n].X, boundryPoints[n].Y, boundryPoints[n].Z);
}

void impactParam::setBoundry_old(double minKin)
{
        /*
        This functions generates the potential energy boundry
        that will be fit to an ellipsoid. Any particle outside this
        boundry (and thus the ellipsoid) is defined to have a
        potential energy deemed to be zero (by approximation).
        The zero potential energy is chosen such that the kinetic
        energy of the slowest possible moving gas particle is only
        lowered by 1%. The boundry thus defines the starting and
        ending position of all gas particle trajectories.
        */

        int numRotations = 20;                                  //number of rotations about z-axis
    double theta;                                                       //holds angle of rotation
    vector3D<double> vel, pos, prevPos;         //vectors for determining test points
    double velMag = 0.1;                                        //incrimental vector
        double posMag =
                max(mol->maxX, mol->maxY) + 15.0;       //starting distance from z-axis
        double startX, startY;                                  //start positions
    double zPos = 0.0, zMax = 0.0;                      //distance along z-axis, maximum extent along z-axis
    double energy;                                                      //holds potential at pos
    double minKinEnergy = 0.5*4*2*2;            //minimum kinetic energy value
    double slope = 1.0;                                         //determins search directions
    double zInc = 0.0, maxZLevels = 25;         //z-axis incrimental values
    bool accept = false, hit = false;           //allows exit from while-loops

    //first, find maximum extent along z-axis
    accept = false;
    pos.setEqual(0, 0, abs(mol->maxZ) + 20);
    energy = forceFunct.getPotential(pos);

    if(energy + 0.01*minKinEnergy < 0)
    {
      slope = -1;
    }
    while (!accept)
    {
      prevPos = pos;
      pos.Z -= 0.5*slope;
      energy = forceFunct.getPotential(pos);
          //found closest distance along z-axis
      if(energy + 0.01*minKinEnergy > 0)
        {
          zMax = max(pos.Z, prevPos.Z);
          accept = true;
        }
    }

        //start search from bottom-up
      zPos = -zMax;
      zInc = 2*zMax/(maxZLevels);
      while(zPos < zMax)
      {
        zPos += 3;

                //loop through each rotation about z-axis
        for(int n = 0; n < numRotations; n ++)
        {
                        //create starting positions
          slope = 1;
          theta = n*2*M_PI/((double)numRotations);
          startX = posMag*cos(theta);
          startY = posMag*sin(theta);
          pos.setEqual(startX, startY, zPos);
          vel.setEqual(-velMag*cos(theta), -velMag*sin(theta), 0.0);
          hit = false;

                  //get potential energy at pos
          energy = forceFunct.getPotential(pos);

                  //if within the boundry that is being set up,
                  //search from inside-out instead
          if(energy + 0.01*minKinEnergy < 0)
            {
              slope = -1;
            }
          vel = vel*slope;

          while(!hit)
          {
            energy = forceFunct.getPotential(pos);



                        //potential energy vs kinetic energy condition
            if((energy + 0.01*minKinEnergy)*slope <  0)
              {
                hit = true;

                                //choose point with farthest distance
                if(prevPos.mag() > pos.mag())
                                        boundryPoints.push_back(prevPos);
                else
                                        boundryPoints.push_back(pos);
              }

                        //if energy condition was not met,
                        //but traversed across other side of molecule
            else if((startX == 0 && startY/pos.Y <= 0)
                || (startY == 0 && startX/pos.X <= 0)
                || (startX/pos.X <= 0 && startY/pos.Y <= 0))
            {
                hit = accept = true;
                                boundryPoints.push_back(pos);
            }

                        //update position vectors
            else
            {
              prevPos = pos;
              pos = pos + vel;
            }
          }
        }

      }
}

void impactParam::findEllipsoidBound()
{

	//stores summation information
	double mu = 0, nu = 0, epsilon = 0, delta = 0, sigma = 0, rho = 0;
	double alpha = 0, beta = 0, gamma = 0;
	double A = 0, B = 0, C = 0;
	rmsd = 0;

	//variables for ellipsoid surface area
	double p, f, deriv_a, deriv_b, deriv_c;

	//stores uncertainty information
	aUncert = 0, bUncert = 0, cUncert = 0;

	//double x2 = 0; double y2 = 0; double z2 = 0;
	//double x4 = 0; double y4 = 0; double z4 = 0;

	double x2 = 0, y2 = 0, z2 = 0, x4 = 0, y4 = 0,
		z4 = 0, x2y2 = 0, x2z2 = 0, y2z2 = 0;

	//calculate sumamtion over data coordinates
	//cout << endl;
	for (int i = 0; i < boundryPoints.size(); i++)
	{
		x2 = boundryPoints[i].X*boundryPoints[i].X;
		y2 = boundryPoints[i].Y*boundryPoints[i].Y;
		z2 = boundryPoints[i].Z*boundryPoints[i].Z;

		//cout << "Boundry: " << boundryPoints[i].X << "  " << boundryPoints[i].Y << "  " << boundryPoints[i].Z << endl;


		x4 = x2*x2; y4 = y2*y2; z4 = z2*z2;
		x2y2 = x2*y2; x2z2 = x2*z2; y2z2 = y2*z2;

		mu = mu + x4;
		nu = nu + y4;
		epsilon = epsilon + z4;

		delta = delta + x2*y2;
		sigma = sigma + x2*z2;
		rho = rho + y2*z2;

		alpha = alpha + x2;
		beta = beta + y2;
		gamma = gamma + z2;


	}


	double det = epsilon*mu*nu - epsilon*delta*delta - mu*rho*rho + 2 * delta*rho*sigma - nu*sigma*sigma;

	A = (alpha*(nu*epsilon - rho*rho)
		+ delta*(rho*gamma - beta*epsilon)
		+ sigma*(beta*rho - nu*gamma)) / det;

	B = (mu*(beta*epsilon - gamma*rho)
		+ alpha*(rho*sigma - delta*epsilon)
		+ sigma*(delta*gamma - beta*sigma)) / det;

	C = (mu*(nu*gamma - beta*rho)
		+ delta*(beta*sigma - delta*gamma)
		+ alpha*(delta*rho - nu*sigma)) / det;

	a = 1 / sqrt(A); b = 1 / sqrt(B); c = 1 / sqrt(C);

	rmsd = sqrt((A*A*mu + B*B*nu + C*C*epsilon +
		2 * (A*B*delta + B*C*rho + A*C*sigma - A*alpha - B*beta - C*gamma) + boundryPoints.size()) / (boundryPoints.size() - 1));

	cov_aa = 2 * rmsd*rmsd* (epsilon*nu - rho*rho)*a*a     / (4 * A*A * det);
	cov_bb = 2 * rmsd*rmsd* (epsilon*mu - sigma*sigma)*b*b / (4 * B*B * det);
	cov_cc = 2 * rmsd*rmsd* (mu*nu - delta*delta)*c*c      / (4 * C*C * det);

	cov_ab = 2 * rmsd*rmsd* (rho*sigma - delta*epsilon)*a*b / (4 * A*B * det);
	cov_ac = 2 * rmsd*rmsd* (delta*rho - nu*sigma)*a*c      / (4 * A*C * det);
	cov_bc = 2 * rmsd*rmsd* (delta*sigma - mu*rho)*b*c      / (4 * B*C * det);

	aUncert = (a / A)*sqrt((epsilon*nu - rho*rho) * 2 / det)*rmsd/2;
	bUncert = (b / B)*sqrt((epsilon*mu - sigma*sigma) * 2 / det)*rmsd/2;
	cUncert = (c / C)*sqrt((mu*nu - delta*delta) * 2 / det)*rmsd/2;

	avgAxes = (a + b + c) / 3.0;
	avgAxesUncert = sqrt(cov_aa + cov_bb + cov_cc + 2 * (cov_ab + cov_ac + cov_bc)) / 3.0;

	p = 1.6075;
	f = (pow(a*b, p) + pow(a*c, p) + pow(b*c, p)) / 3;
	surfArea = 4 * M_PI*pow(f, 1.0/p);
	deriv_a = surfArea*pow(a, p - 1)*(pow(b, p) + pow(c, p)) / (3 * f);
	deriv_b = surfArea*pow(b, p - 1)*(pow(a, p) + pow(c, p)) / (3 * f);
	deriv_c = surfArea*pow(c, p - 1)*(pow(a, p) + pow(b, p)) / (3 * f);

	surfAreaUncert = sqrt(deriv_a*deriv_a*cov_aa + deriv_b*deriv_b*cov_bb + deriv_c*deriv_c*cov_cc
		+ 2 * (deriv_a*deriv_b*cov_ab + deriv_a*deriv_c*cov_ac + deriv_b*deriv_c*cov_bc));

}

void impactParam::printBoundryPoints(string outFileName)
{
	ofstream energyPoints;
	energyPoints.open(outFileName.c_str());
	for (int n = 0; n < boundryPoints.size(); n++)
		energyPoints << boundryPoints[n].X << " " << boundryPoints[n].Y << "  " << boundryPoints[n].Z << endl;
	energyPoints.close();
}

void impactParam::enlargeEllipBoundry()
{
        double aOld = a, bOld = b, cOld = c;

	bool accept = false;
	//double aInc = a*0.01;
	//double bInc = b*0.01;
	//double cInc = c*0.01;
	double aInc = 0.01;
	double bInc = 0.01;
	double cInc = 0.01;
	int i = 0;

	double mag = 0;


	int counter = 0;

	// continue to enlarge ellipsoid until all
	//boundry points are inside
	do
	{
		counter++;
		//loop through all boundry points
		for (i = 0; i < (signed)boundryPoints.size() + 1; i++)
		{
			//incase running in an infinite loop
			//if (i == (signed)boundryPoints.size() || a > aInc * 300)
		        if(i == (signed)boundryPoints.size() || a > 3*aOld)
			{
				accept = true;
				break;
			}

			//definition of al ellipsoid: (x/z)^2 + (y/b)^2 + (z/c)^2 = 1
			mag = boundryPoints[i].X*boundryPoints[i].X / (a*a)
				+ boundryPoints[i].Y*boundryPoints[i].Y / (b*b)
				+ boundryPoints[i].Z*boundryPoints[i].Z / (c*c);
			if(mag > 1)
			{
				a = a + aInc;
				b = b + bInc;
				c = c + cInc;
				break;
			}
		}
	} while (!accept);

	//enlargePercentage = (abs(a - (aInc/0.01) )/a)*100;
	//enlargePercentage = ((a*b*c)/(aOld*bOld*cOld) - 1.0)*100.0;

	//enlargePercentage = 100*(a/aOld - 1);
	enlargeAmount = a - aOld;
}

void impactParam::printBoundryResults()
{
	char buffer[100];
	char output[1000];

	sprintf(buffer, "\n Impact parameter properties:\n");
	strcpy(output, buffer);

	sprintf(buffer, "                      X         Y         Z\n");
	strcat(output, buffer);

	sprintf(buffer, " -----------------------------------------------\n");
	strcat(output, buffer);

	sprintf(buffer, " Axis Length    : %9.5G %9.5G %9.5G\n", a, b, c);
	strcat(output, buffer);

	sprintf(buffer, " Uncertainty    : %9.5G %9.5G %9.5G\n", aUncert, bUncert, cUncert);
	strcat(output, buffer);

	sprintf(buffer, " Percent Uncert : %9.5G %9.5G %9.5G\n", aUncert * 100 / a, bUncert * 100 / b, cUncert * 100 / c);
	strcat(output, buffer);

	sprintf(buffer, " -----------------------------------------------\n");
	strcat(output, buffer);

        sprintf(buffer, "\t RMSD: %10.5G\n", rmsd);
        strcat(output, buffer);

	sprintf(buffer, " \t Covariance matrix:\n");
	strcat(output, buffer);

	sprintf(buffer, "\t { %*.4g  %*.4g  %*.4g   }\n", 10, cov_aa, 10, cov_ab, 10, cov_ac);
	strcat(output, buffer);

	sprintf(buffer, "\t { %*.4g  %*.4g  %*.4g   }\n", 10, cov_ab, 10, cov_bb, 10, cov_bc);
	strcat(output, buffer);

	sprintf(buffer, "\t { %*.4g  %*.4g  %*.4g   }\n", 10, cov_ac, 10, cov_bc, 10, cov_cc);
	strcat(output, buffer);


	sprintf(buffer, " Average ellipsoid axes: %9.5G +/- %9.5G Ang.\n", avgAxes, avgAxesUncert);
	strcat(output, buffer);

	sprintf(buffer, " Ellipsoid surface area: %9.5G +/- %9.5G Ang^2.\n", surfArea, surfAreaUncert);
	strcat(output, buffer);

	enlargeEllipBoundry();

	sprintf(buffer, " Ellipsoid axes further expanded by %6.2f Ang.\n", enlargeAmount);
	strcat(output, buffer);

	printf("\n%s\n", output);

	fileStream_ext.open(omegaFilePath_ext.c_str(), std::ios_base::app);
	fileStream_ext << endl << output << endl;
	fileStream_ext.close();

}


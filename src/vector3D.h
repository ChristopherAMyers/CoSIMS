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

#define _USE_MATH_DEFINE
#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#include <cmath>
#include <iostream>
using namespace std;


template <class type> class vector3D
{
public:
	vector3D(void);
	vector3D(type, type, type);
	vector3D(type[3]);
	~vector3D(void);

	//cartesian coordinates
	type X;
	type Y;
	type Z;

	//spherical coordinates
	//type R;
	//type theta;
	//type phi;

	//cylindrical coordinates
	//type S;
	//type phi;
	//type ZC;

	void setEqual(vector3D&); //sets this vecotr equal to another
	void setEqual(type, type, type);
	type mag();//returns the magnitude of this vector
	type mag2();//returns the magnitude of this vector squared
	
	double Angle(vector3D&);//returns the angle between this vector and another
	type Angle(type(&)[3]);
	double dot(vector3D&);//returns the dot-product between this vector and another;
	vector3D<type> cross(vector3D&);//returns the cross-product between this vector and another
	
	void carToSph();//sets current vectors spherical coordinates from its cartesian
	void sphToCar();//sets current vectors cartesian coordinates from its spherical
	type getPhi(); //gets the current vectors azimuthal angle
	type getS(); //gets the current vectors radial component (cylindrical coords)
	type getTheta(); //gets the current vectors polar angle
	type getR(); //gets the current vectors radial component (spherical coords)
	type at(int); //returns a specific component
	void rotate(vector3D&, type, type, type); //rotates the given vector about Euler angles
	void rotate(type, type, type); //rotates THIS vector about Euler angles
	void clear(); //clear vector ans set components to 0
	void rotateX(type); //rotate about x-axis
	void rotateY(type); //rotate about y-axis
	void rotateZ(type); //rotate about z-axis
	void rotateVec(vector3D&, type);
	void rotate(double[3][3]); //rotates vector with a given 3x3 roration matrix
	void rotate2(double[3][3]); //rotates vector with a given 3x3 roration matrix
	
	//overloaded operators
	void operator = (const vector3D&);
	vector3D<type> operator + (const vector3D&);
	vector3D<type> operator + (const type);
	vector3D<type> operator - (const vector3D&);
	vector3D<type> operator - (const type);
	vector3D<type> operator - (const type[3]);
	vector3D<type> operator * (const type);
	vector3D<type> operator / (const type);
	vector3D<type> operator / (const vector3D<type>&);
	

};




template <class type> vector3D<type>::vector3D()
{
	X = 0;
	Y = 0;
	Z = 0;

	//R = 0; theta = 0; phi = 0; S = 0;
}

template <class type> vector3D<type>::vector3D(type x, type y, type z)
{
	X = x;
	Y = y;
	Z = z;

	//R = 0; theta = 0;  phi = 0; S = 0;
}

template <class type> vector3D<type>::vector3D(type vec[3])
{
        X = vec[0];
        Y = vec[1];
        Z = vec[2];

        //R = 0; theta = 0;  phi = 0; S = 0;
}

template <class type> vector3D<type>::~vector3D()
{
  /*
	X = 0;
	Y = 0;
	Z = 0;

	R = 0; theta = 0; phi = 0; S = 0;
*/
}

template <class type> type vector3D<type>::mag()
{
	return sqrt((X*X) + (Y*Y) + (Z*Z));
}
template <class type> type vector3D<type>::mag2()
{
        return (X*X) + (Y*Y) + (Z*Z);
}

template <class type> void vector3D<type>::setEqual(vector3D &vec)
{
	X = vec.X;
	Y = vec.Y;
	Z = vec.Z;
}

template <class type> void vector3D<type>::setEqual(type x, type y, type z)
{
	X = x;
	Y = y;
	Z = z;
}

template <class type> double vector3D<type>::Angle(vector3D &vec)
{
	double term1 = dot(vec);
	double term2 = mag()*vec.mag();
	double term3 = term1 / term2;
	if (term3 > 1) return 0;
	else if (term3 < -1) return M_PI;

	return acos(term1 / term2);
}

template <class type> type vector3D<type>::Angle(type(&inVec)[3])
{
	type term1 = X*inVec[0] + Y*inVec[1] + Z*inVec[2];
	type term2 = mag()
		*sqrt(inVec[0] * inVec[0] + inVec[1] * inVec[1] + inVec[2] * inVec[2]);
	type term3 = term1 / term2;

	//accounts for numerical error
	if (term3 > 1) return 0;
	else if (term3 < -1) return M_PI;

	return acos(term1 / term2);
}

template <class type> vector3D<type> vector3D<type>::cross(vector3D &vec)
{
	double term1 = ((Y*vec.Z) - (Z*vec.Y));
	double term2 = ((Z*vec.X) - (X*vec.Z));
	double term3 = ((X*vec.Y) - (Y * vec.X));
	vector3D retVec(term1, term2, term3);
	return retVec;
}

template <class type> double vector3D<type>::dot(vector3D &vec)
{
	return ((X*vec.X) + (Y*vec.Y) + (Z*vec.Z));
}

/*
template <class type> void vector3D<type>::carToSph()
{
	R = sqrt((X*X) + (Y*Y) + (Z*Z));
	theta = acos(Z / R);
	phi = atan2(Y, X);
	if (phi < 0)
		phi = phi + 0;
}

template <class type> void vector3D<type>::sphToCar()
{
	X = R*sin(theta)*cos(phi);
	Y = R*sin(theta)*sin(phi);
	Z = R*cos(theta);
}


template <class type> type vector3D<type>::getPhi()
    {
      double tempPhi = atan2(Y,X);
      if (phi < 0)
        phi = phi + 0;
      return tempPhi;
    }


template <class type> type vector3D<type>::getTheta()
    {
      return acos(Z/mag());
    }

template <class type> type vector3D<type>::getR()
    {
      return sqrt((X*X) + (Y*Y) + (Z*Z));
    }

template <class type> type vector3D<type>::getS()
    {
      return sqrt((X*X) + (Y*Y));
    }

*/
template <class type> void vector3D<type>::rotateVec(vector3D& k, type theta)
{
    double x = X; double y = Y; double z = Z;

    vector3D<type> tempVec;

    tempVec.X = X;
    tempVec.Y = Y;
    tempVec.Z = Z;

    //tempVec.X = k.X*tempVec.dot(k)*(1 - cos(theta)) + X*cos(theta) + (k.Y*Z - k.Z*Y)*sin(theta);

    *this = *this*cos(theta) - this->cross(k)*sin(theta) + k*this->dot(k)*(1 - cos(theta));

    //X = tempVec.X; Y = tempVec.Y; Z = tempVec.Z;

}

/*template <class type> void vector3D<type>::rotateVec(vector3D& k, type theta)
{
    double x = X; double y = Y; double z = Z;

    x = k.X*(X*k.X + Y*k.Y + Z*k.Z)*(1 - cos(theta)) + X*cos(theta) + (k.Y*Z - k.Z*Y)*sin(theta);
    y = k.Y*(X*k.X + Y*k.Y + Z*k.Z)*(1 - cos(theta)) + Y*cos(theta) + (k.Z*X - k.X*Z)*sin(theta);
    z = k.Z*(X*k.X + Y*k.Y + Z*k.Z)*(1 - cos(theta)) + Z*cos(theta) + (k.X*Y - k.Y*X)*sin(theta);


    X = x; Y = y; Z = z;

}*/

template <class type> void vector3D<type>::rotate(type theta, type phi, type gamma)
{
	//rotates this vector about three Euler angles

	//holds original coordinates
	double x = X;
	double y = Y;
	double z = Z;

	//construct quaternion components
	double q0 = cos(0.5*theta)*cos(0.5*(phi + gamma));
	double q1 = sin(0.5*theta)*cos(0.5*(phi + gamma));
	double q2 = sin(0.5*theta)*sin(0.5*(phi + gamma));
	double q3 = cos(0.5*theta)*sin(0.5*(phi + gamma));

	X = ((q0*q0 + q1*q1 - q2*q2 - q3*q3)*x)
		+ (2 * (q1*q2 + q0*q3)*y)
		+ (2 * (q1*q3 - q0*q2)*z);

	Y = (2 * (q1*q2 - q0*q3)*x)
		+ ((q0*q0 - q1*q1 + q2*q2 - q3*q3)*y)
		+ (2 * (q2*q3 + q0*q1)*z);

	Z = (2 * (q1*q3 + q0*q2)*x)
		+ (2 * (q2*q3 - q0*q1)*y)
		+ ((q0*q0 - q1*q1 - q2*q2 + q3*q3)*z);

}

template <class type> void vector3D<type>::rotate(vector3D &vec, type theta, type phi, type gamma)
{
	//Takes in a vector and rotates it about three Euler angles


	double x = vec.X; double y = vec.Y; double z = vec.Z;

	//construct quaternion components
	double q0 = cos(0.5*(double)theta)*cos(0.5*((double)phi + (double)gamma));
	double q1 = sin(0.5*(double)theta)*cos(0.5*((double)phi + (double)gamma));
	double q2 = sin(0.5*(double)theta)*sin(0.5*((double)phi + (double)gamma));
	double q3 = cos(0.5*(double)theta)*sin(0.5*((double)phi + (double)gamma));

	//old
	X = ((q0*q0 + q1*q1 - q2*q2 - q3*q3)*x) + (2 * (q1*q2 + q0*q3)*y) + (2 * (q1*q3 - q0*q2)*z);
	Y = (2 * (q1*q2 - q0*q3)*x) + ((q0*q0 - q1*q1 + q2*q2 - q3*q3)*y) + (2 * (q2*q3 + q0*q1)*z);
	Z = (2 * (q1*q3 + q0*q2)*x) + (2 * (q2*q3 - q0*q1)*y) + ((q0*q0 - q1*q1 - q2*q2 + q3*q3)*z);

}

template <class type> void vector3D<type>::rotate(double matrix[3][3])
{

  //computes the matrix multiplication M*x
  //where M is a 3x3 matrix and x is a 1x3 vector
  double x = this->X; double y = this->Y; double z = this->Z;

  X = x*matrix[0][0] + y*matrix[0][1] + z*matrix[0][2];
  Y = x*matrix[1][0] + y*matrix[1][1] + z*matrix[1][2];
  Z = x*matrix[2][0] + y*matrix[2][1] + z*matrix[2][2];

}

template <class type> void vector3D<type>::rotate2(double matrix[3][3])
{

  //computes the matrix multiplication x*M
  //where M is a 3x3 matrix and x is a 1x3 vector
  double x = X; double y = Y; double z = Z;

  //X = x*matrix[0][0] + y*matrix[1][0] + z*matrix[2][0];
  //Y = x*matrix[0][1] + y*matrix[1][1] + z*matrix[2][1];
  //Z = x*matrix[0][2] + y*matrix[1][2] + z*matrix[2][2];

  X = x*matrix[0][0] + y*matrix[0][1] + z*matrix[0][2];
  Y = x*matrix[1][0] + y*matrix[1][1] + z*matrix[1][2];
  Z = x*matrix[2][0] + y*matrix[2][1] + z*matrix[2][2];

}

/*///////////////////////////////////////////////////////////
//// These next three functions
//// Rotates the molecule using the rotation matricies
////
////	  | 1     0      0    || X |   | X'|
//// Rx = | 0  cos(a) -sin(a) || Y | = | Y'|
////	  | 0  sin(a)  cos(a) || Z |   | Z'|
////
////	  | cos(b)  0  sin(b) || X |   | X'|
//// Ry = | 0       1    0    || Y | = | Y'|
////	  |-sin(b)  0  cos(b) || Z |   | Z'|
////
////	  | cos(c)  -sin(c)  0 || X |   | X'|
//// Rz = | sin(c)   cos(c)  0 || Y | = | Y'|
////	  |   0       0      1 || Z |   | Z'|
////
//// where a, b, and c is the amount of rotation about the axis,
//// on the intervals (-pi, pi), (-pi/2, pi/2), and (-pi, pi), respectively,
//// X, Y, and Z are the original coordinates,
//// and X', Y', and Z' are the rotaed coordinates
//////////////////////////////////////////////////////////*/

template <class type> void vector3D<type>::rotateX(type theta)
{
	double y = Y*cos(theta) - Z*sin(theta);
	double z = Y*sin(theta) + Z*cos(theta);
	Y = y; Z = z;
}

template <class type> void vector3D<type>::rotateY(type theta)
{
	double x = X*cos(theta) + Z*sin(theta);
	double z = -X*sin(theta) + Z*cos(theta);
	X = x; Z = z;
}

template <class type> void vector3D<type>::rotateZ(type theta)
{
	double x = X*cos(theta) - Y*sin(theta);
	double y = X*sin(theta) + Y*cos(theta);
	X = x; Y = y;
}

template <class type> void vector3D<type>::clear()
{
	X = 0; Y = 0; Z = 0;
}

template <class type> type vector3D<type>::at(int i)
{
	if (i == 0) return X;
	else if (i == 1) return Y;
	else return Z;
}

template <class type> void vector3D<type>::operator =(const vector3D &vec)
{
	X = vec.X;
	Y = vec.Y;
	Z = vec.Z;
}

template <class type> vector3D<type> vector3D<type>::operator+ (const vector3D &vec)
{
	vector3D temp;
	temp.X = X + vec.X;
	temp.Y = Y + vec.Y;
	temp.Z = Z + vec.Z;
	return temp;
}

template <class type> vector3D<type> vector3D<type>::operator +(const type myVar)
{
	vector3D temp;
	temp.X = X + (double)myVar;
	temp.Y = Y + (double)myVar;
	temp.Z = Z + (double)myVar;
	return temp;
}

template <class type> vector3D<type> vector3D<type>::operator -(const vector3D &vec)
{
	vector3D temp;
	temp.X = X - vec.X;
	temp.Y = Y - vec.Y;
	temp.Z = Z - vec.Z;
	return temp;
}

template <class type> vector3D<type> vector3D<type>::operator -(const type myVar)
{
	vector3D temp;
	temp.X = X - (double)myVar;
	temp.Y = Y - (double)myVar;
	temp.Z = Z - (double)myVar;
	return temp;
}

template <class type> vector3D<type> vector3D<type>::operator -(const type myVar[3])
{
        vector3D temp;
        temp.X = X - (double)myVar[1];
        temp.Y = Y - (double)myVar[2];
        temp.Z = Z - (double)myVar[3];
        return temp;
}

template <class type> vector3D<type> vector3D<type>::operator *(const type myVar)
{
	vector3D temp;
	temp.X = X * (double)myVar;
	temp.Y = Y * (double)myVar;
	temp.Z = Z * (double)myVar;
	return temp;
}

template <class type> vector3D<type> vector3D<type>::operator /(const type myVar)
{
	vector3D temp;
	temp.X = X / (double)myVar;
	temp.Y = Y / (double)myVar;
	temp.Z = Z / (double)myVar;
	return temp;
}

template <class type> vector3D<type> vector3D<type>::operator / (const vector3D<type> &vec)
{
  vector3D temp;
  temp.X = X / vec.X;
  temp.Y = Y / vec.Y;
  temp.Z = Z / vec.Z;
  return temp;

  std::cout << "hereInVec: " << temp.X << std::endl;
  return temp;
}


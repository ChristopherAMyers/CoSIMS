#pragma once
#include <math.h>
#include <iostream>
class eigen3D
{
public:
	eigen3D();
	~eigen3D();

	//provide a matrix and save it's eigenvalues
	void getEigenvalues(double[3][3], double(&)[3]);

	//provide a matrix and it's eigenvalues and save the eigenvectors
	void getEigenvectors(double[3][3], double[3], double(&)[3][3]);
};


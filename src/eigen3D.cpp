#include "eigen3D.h"



eigen3D::eigen3D()
{
}


eigen3D::~eigen3D()
{
}

void eigen3D::getEigenvalues(double A[3][3], double(&eigenvalues)[3])
{
	//This function saves the eigenvalues of 3x3 matrix 'A'
	//To the double array 'eigenvalues'.

	double a0, a1, a2;
	double p, q, R, Q, D, a, b, phase, mag;

	//roots of the characteristic equation
	a0 = A[0][0] * A[1][1] * A[2][2] -
		A[0][0] * A[1][2] * A[2][1] -
		A[1][1] * A[0][2] * A[2][0] -
		A[2][2] * A[0][1] * A[1][0] +
		A[0][1] * A[1][2] * A[2][0] +
		A[0][2] * A[2][1] * A[1][0];

	a1 = A[0][1] * A[1][0] + A[0][2] * A[2][0] + A[1][2] * A[2][1] - A[0][0] * A[1][1] - A[1][1] * A[2][2] -
		A[2][2] * A[0][0];

	a2 = A[0][0] + A[1][1] + A[2][2];
	a0 = -a0; a1 = -a1; a2 = -a2;
	
	//solve the cubic equation for 3 real roots
	p = (3 * a1 - a2*a2) / 3;
	q = (9 * a1*a2 - 27 * a0 - 2 * a2*a2*a2) / 27;
	Q = p / 3; R = q / 2; D = Q*Q*Q + R*R;

	if (D< 0)
	{
		a = R; b = sqrt(-D);
	}
	else
	{
		a = R + sqrt(D);
	}
	
	phase = atan2(b, a);

	mag = pow(a*a + b*b, 1.0 / 6.0);

	eigenvalues[0] = -(a2 / 3) + mag * 2 * cos(phase / 3.0);
	eigenvalues[1] = -(a2 / 3) - mag*(cos(phase / 3.0)
		+ sqrt(3)*sin(phase / 3.0));
	eigenvalues[2] = -(a2 / 3) - mag*(cos(phase / 3.0)
		- sqrt(3)*sin(phase / 3.0));
}

void eigen3D::getEigenvectors(double A[3][3], double lambda[3], double(&eigenvectors)[3][3])
{
	//This function saves the verctors of 3x3 matrix 'A'
	//To the 3x3 double array 'eigenvectors'.
	//WARNING: This function does not give orthogonal
	//eigenvectors for so similar eigenvalues, but is sufficient
	//for the PCM clustering implimented in the program.

	double det, y, z, normFactor;

	for (int n = 0; n < 3; n++)
	{
		det = lambda[n] * A[0][2] - A[0][2] * A[1][1] + A[0][1] * A[1][2];
		y = (A[0][2] * A[1][0] + lambda[n] * A[1][2] - A[0][0] * A[1][2]) / det;
		z = (lambda[n] * lambda[n] - lambda[n] * (A[0][0] + A[1][1]) - A[0][1] * A[1][0] + A[0][0] * A[1][1]) / det;

		normFactor = sqrt(1 + y*y + z*z);
		eigenvectors[0][n] = 1.0 / normFactor;
		eigenvectors[1][n] = y / normFactor;
		eigenvectors[2][n] = z / normFactor;

	}

	double transPose[3][3];
	for (int n = 0; n < 3; n++)
		for (int m = 0; m < 3; m++)
			transPose[n][m] = eigenvectors[m][n];

	for (int n = 0; n < 3; n++)
		for (int m = 0; m < 3; m++)
			eigenvectors[n][m] = transPose[n][m];

}

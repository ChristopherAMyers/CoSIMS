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


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
#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



   class VelocityDistribution
  {
    int n; //distribution power
    double a; //distribution parameter
    double mean_v;
    double mode_v;
    double stdDev_v;
    double variance_v;
    double norm_v;
    double range_a, range_b;

    int factorial(int);
    int d_factorial(int);

	

  public:
    VelocityDistribution();
    VelocityDistribution(int, double);
    VelocityDistribution(int, double, double, double);
    virtual ~VelocityDistribution();
    void setProperties(int, double);
    void renormalize(double range_A, double range_B);
    double normFactor();
    double normFactor(int, double);
    double normFactor(int, double, double, double);
	int power();
    double mean();
    double mean(int, double);
    double mode();
    double mode(int, double);
    double stdDev();
    double stdDev(int, double);
    double variance();
    double variance(int, double);
    double cdf(double);
    double cdf(int, double, double); 

    double operator()(double);

  };



/*
 * VelocityDistribution.h
 *
 *  Created on: Sep 20, 2017
 *      Author: cmyers
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



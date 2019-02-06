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

/*
 *  This template class defines the Maxwell-Boltzmann Distribution
 *  Note that this is build using the Boltzmann constant in units
 *  of (Angstrom^2)(atomic mass unit)/(pico-seconds^2)(Kelvin)
 */

#ifndef MAXBOLTZ_H_
#define MAXBOLTZ_H_
#define K 0.831446214563
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
using namespace std;

template <class type> class maxBoltz
{
public:
	type temp;
	type mass;
	type mean;
	type mode;
	type variance;
	type a; //scale parameter
	type normConst; //normalization constant
	type distMax; //maximum value of the distribution, evaluated at the mode


        maxBoltz();
	maxBoltz(type, type);
	virtual ~maxBoltz(void);

	//set temperature and mass parameters
	void setTemp(type);
	void setMass(type);

        //set the distro specials once mass or temp has changed
        void setSpecials();
        void setSpecials(type, type);

	//sets normalization cinstant over the interval [A,B];
	void setNormalization(type, type);

	//get values of temperature and mass parameters
	type getTemp();
	type getMass();
	type getMean();
	type getMode();
	type getvariance();
	type getStdDev();
	type getMaxDist();

	//evaluates the distribution at the current velocity
	//this distro is normalized over the interval [0, infinity]
	type distro(type);

	//returns the CDF as a function of it's ending value;
	type cdf(type);

	//evaluates the re-normalized distribution at the current velocity
	type distroN(type);

};

template <class type> maxBoltz<type>::maxBoltz()
    {
      temp = 298.0;
      mass = 1.0;
      setSpecials();
    }

template <class type> maxBoltz<type>::maxBoltz(type T, type m)
    {
      temp = T;
      mass = m;
      setSpecials();
    }

template <class type> maxBoltz<type>::~maxBoltz()
    {

    }

template <class type> type maxBoltz<type>::distro(type v)
    {
      return sqrt(M_2_PI)*v*v*exp(-v*v/(2*a*a))/(a*a*a);
    }

template <class type> type maxBoltz<type>::distroN(type v)
    {
      return distro(v)/normConst;
    }

template <class type> type maxBoltz<type>::cdf(type x)
    {
        return erf(x*M_SQRT1_2/a) - sqrt(M_2_PI)*x*exp(-x*x/(2*a*a))/a;
    }

template <class type> void maxBoltz<type>::setTemp(type T)
    {
      temp = T;
      setSpecials();
    }

template <class type> void maxBoltz<type>::setMass(type m)
    {
      mass = m;
      setSpecials();
    }

template <class type> void maxBoltz<type>::setSpecials()
    {

      a = sqrt(K*temp/mass);
      mean = 2*a*sqrt(M_2_PI);
      mode = M_SQRT2*a;
      variance = a*a*(3*M_PI - 8)/M_PI;
      distMax = sqrt(M_2_PI)*2/(a*M_E);
    }

template <class type> void maxBoltz<type>::setSpecials(type m, type T)
    {

      mass = m;
      temp = T;
      a = sqrt(K*temp/mass);
      mean = 2*a*sqrt(M_2_PI);
      mode = M_SQRT2*a;
      variance = a*a*(3*M_PI - 8)/M_PI;
      distMax = sqrt(M_2_PI)*2/(a*M_E);
    }

template <class type> void maxBoltz<type>::setNormalization(type A, type B)
    {
      normConst = cdf(B) - cdf(A);
    }

template <class type> type maxBoltz<type>::getMaxDist()
    {
      return distMax;
    }

template <class type> type maxBoltz<type>::getTemp()
    {
      return temp;
    }

template <class type> type maxBoltz<type>::getMass()
    {
      return mass;
    }

template <class type> type maxBoltz<type>::getMean()
    {
      return 2*a*sqrt(M_2_PI);
    }

template <class type> type maxBoltz<type>::getMode()
    {
      return M_SQRT2*a;
    }

template <class type> type maxBoltz<type>::getvariance()
    {
      return a*a*(3*M_PI - 8)*M_1_PI;
    }

template <class type> type maxBoltz<type>::getStdDev()
    {
      return a*sqrt((3*M_PI - 8)*M_1_PI);
    }



#endif /* MAXBOLTZ_H_ */














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

#include "VelocityDistribution.h"

VelocityDistribution::VelocityDistribution()
{
    setProperties(n, a);
}

VelocityDistribution::VelocityDistribution(int power, double a_in)
{
    setProperties(n, a);
}

VelocityDistribution::VelocityDistribution(int power, double a_int, double range_1, double range_2)
{
    setProperties(n, a);

    range_a = range_1; range_b = range_2;
    norm_v = normFactor(n, a, range_1, range_2);
}

VelocityDistribution::~VelocityDistribution()
{

}

void VelocityDistribution::setProperties(int power, double param)
{
    n = power;
    a = param;

    norm_v = normFactor(n, a);
    mode_v = mode(power, param);
    mean_v = mean(power, param);
    stdDev_v = stdDev(power, param);
    variance_v = variance(power, param);
}

double VelocityDistribution::normFactor()
{
    return norm_v;
}

void VelocityDistribution::renormalize(double range_A, double range_B)
{
    norm_v = norm_v / (cdf(range_B) - cdf(range_A));
}

double VelocityDistribution::normFactor(int power, double param)
{
    return 2 * pow(param, (power + 1) / 2.0) / tgamma((n + 1) / 2.0);
}

double VelocityDistribution::normFactor(int power, double param, double range_A, double range_B)
{
    return normFactor(power, param) / (cdf(power, param, range_B) - cdf(power, param, range_A));
}

int VelocityDistribution::power()
{
    return n;
}

double VelocityDistribution::mean()
{
    return mean_v;
}

double VelocityDistribution::mean(int power, double param)
{
    return tgamma((power + 2)/2.0)/(tgamma((power + 1)/2.0)*sqrt(param));
}

double VelocityDistribution::mode()
{
    return mode_v;
}

double VelocityDistribution::mode(int power, double param)
{
    return sqrt(power/(2*param));
}

double VelocityDistribution::variance()
{
    return variance_v;
}

double VelocityDistribution::variance(int power, double param)
{
    return ((power + 1)/2.0 - pow(tgamma((power + 2)/2.0)/tgamma((power + 1)/2.0), 2))/param;
}

double VelocityDistribution::stdDev()
{
    return stdDev_v;
}

double VelocityDistribution::stdDev(int power, double param)
{
    return sqrt(variance(power, param));
}

double VelocityDistribution::cdf(double x)
{
    return cdf(n, a, x);
}

double VelocityDistribution::cdf(int power, double param, double x)
{
    double sum = 0;

    if (power % 2 == 0)
    {
        for (int k = 0; k <= (power / 2.0 - 1); k++)
            sum += pow(2 * x*x*param, k) / (double)d_factorial(2 * k + 1);
        return erf(sqrt(param)*x) - sqrt(param / M_PI)*exp(-param*x*x) * 2 * x*sum;
    }
    else
    {
        for (int k = 0; k <= power; k++)
        {
            sum += pow(x*x*param, k) / (double)factorial(k);
        }
        return 1 - exp(-param*x*x)*sum;
    }
}

double VelocityDistribution::operator()(double v)
{
    return pow(v, n)*exp(-a*v*v)*norm_v;
}


int VelocityDistribution::factorial(int x)
{
    if(x <= 1)
        return 1;
    else
        return x*factorial(x - 1);
}

int VelocityDistribution::d_factorial(int x)
{
    if(x == 1 || x == 0)
        return 1;
    else
        return x*d_factorial(x - 2);
}





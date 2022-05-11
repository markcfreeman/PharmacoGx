// #pragma once

// #include "stdafx.h"
#include <algorithm>
#include <cmath>
#include <vector>

// #include "<Rcpp.h>"
// [[Rcpp::plugins("cpp11")]]

double effectiveConcentration2(const double myConcentration,
	const double otherConcentration,
	const double myLogD0,
	const double otherLogD0,
	const double logAOtherOnMine,
	const double logAMineOnOther)
{
	const double myD0{ exp(myLogD0) };
	const double otherD0{ exp(otherLogD0) };
	const double AMineOnOther{ exp(logAMineOnOther) };
	const double AOtherOnMine{ exp(logAOtherOnMine) };
	
	double A{ AMineOnOther * otherD0 };
	const double A1{ AOtherOnMine * otherConcentration };
	A += A1;

	double B{ myD0 * otherD0 };
	double B1{ AOtherOnMine * otherConcentration };
	B1 *= myD0;
	B += B1;

	double B2{ AMineOnOther * myConcentration };
	B2 *= otherD0;
	B -= B2;

	const double B3{ myConcentration * otherConcentration };
	B -= B3;

	double C{ myD0 * otherD0 };
	C *= -myConcentration;
	double C1{ myConcentration * otherConcentration };
	C1 *= myD0;
	C -= C1;

	double ans{ B * B };
	double product{ 4.0 * A };
	product *= C;
	ans -= product;
	ans = sqrt(ans);
	ans -= B;
	ans /= A;
	ans *= 0.5;

	return ans;
}

double hill(const double x, const double logD0, const double logN)
{
	double denom{ log(x) };
	denom -= logD0;

	double N{ exp(logN) };
	denom *= N;
	denom = exp(denom);
	++denom;

	double ans{ 1.0 };
	ans /= denom;

	return ans;
}

double katzir2(const double xi, const double xj, const std::vector<double>& params)
{
	double A{ effectiveConcentration2(xi, xj, params.at(0), params.at(1), params.at(4), params.at(5)) };
	double ans{ hill(A, params.at(0), params.at(2)) };

	double B{ effectiveConcentration2(xj, xi, params.at(1), params.at(0), params.at(5), params.at(4)) };
	ans *= hill(B, params.at(1), params.at(3));

	return ans;
}

double katzirResidual(const std::vector<double>& xi,
	const std::vector<double>& xj,
	const std::vector<double>& y,
	const std::vector<double>& params)
{
	double ans{ 0.0 };
	double term;

	for (size_t i = 0; i < xi.size(); ++i)
	{
		term = katzir2(xi.at(i), xj.at(i), params);
		term -= y.at(i);
		term *= term;
		ans += term;
	}

	return ans;
}

// [[Rcpp::export]]
std::vector<double> katzirMeshOptimization(const std::vector<double>& xi,
	const std::vector<double>& xj,
	const std::vector<double>& y,
	const std::vector<double>& start,
	const std::vector<double>& lower,
	const std::vector<double>& upper,
	std::vector<double>& guessIntervals,
	const std::vector<bool>& missingParams)
{
	// NB all vectors must have values in order logD0i, logD0j, logNi, logNj, logAij, logAji
	double bestResidual{ katzirResidual(xi, xj, y, start) };
	double newResidual;
	std::vector<double> bestParams{ start };
	std::vector<double> currentParams{ lower };
	
	for (currentParams.at(0) = lower.at(0); currentParams.at(0) <= upper.at(0); currentParams.at(0) += guessIntervals.at(0))
	{
		for (currentParams.at(1) = lower.at(1); currentParams.at(1) <= upper.at(1); currentParams.at(1) += guessIntervals.at(1))
		{
			for (currentParams.at(2) = lower.at(2); currentParams.at(2) <= upper.at(2); currentParams.at(2) += guessIntervals.at(2))
			{
				for (currentParams.at(3) = lower.at(3); currentParams.at(3) <= upper.at(3); currentParams.at(3) += guessIntervals.at(3))
				{
					for (currentParams.at(4) = lower.at(4); currentParams.at(4) <= upper.at(4); currentParams.at(4) += guessIntervals.at(4))
					{
						for (currentParams.at(5) = lower.at(5); currentParams.at(5) <= upper.at(5); currentParams.at(5) += guessIntervals.at(5))
						{
							newResidual = katzirResidual(xi, xj, y, currentParams);

							if (newResidual < bestResidual)
							{
								bestResidual = newResidual;
								bestParams = currentParams;
							}
						}
					}
				}
			}
		}
	}

	bool reachedLocalMinimum;
	std::vector<double> iterationStartParams;
	std::vector<double> newParams;
	std::vector<double> valuesToTest{ 0, 0 }; // dummy values to create vector of correct size

	for (size_t i = 0; i < guessIntervals.size(); ++i)
	{
		guessIntervals.at(i) *= 0.5;
	}

	do
	{
		reachedLocalMinimum = true;
		iterationStartParams = start;

		for (size_t param = 0; param < start.size(); ++param)
		{
			if (!(missingParams.at(param)))
			{
				newParams = iterationStartParams;
				valuesToTest.at(0) = iterationStartParams.at(param);
				valuesToTest.at(0) -= guessIntervals.at(param);
				valuesToTest.at(0) = std::max(valuesToTest.at(0), lower.at(param));

				valuesToTest.at(1) = iterationStartParams.at(param);
				valuesToTest.at(1) += guessIntervals.at(param);
				valuesToTest.at(1) = std::min(valuesToTest.at(1), upper.at(param));
				
				for (double &value : valuesToTest)
				{
					newParams.at(param) = value;
					newResidual = katzirResidual(xi, xj, y, newParams);

					if (newResidual < bestResidual)
					{
						bestParams = newParams;
						bestResidual = newResidual;
						reachedLocalMinimum = false;
						guessIntervals.at(param) *= 0.5;
					}
				}
			}
		}
	} while (!reachedLocalMinimum);

	return bestParams;
}
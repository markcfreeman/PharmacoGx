#pragma once

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

double effectiveConcentration2(const double,
	const double,
	const double,
	const double,
	const double,
	const double);
double hill(const double, const double, const double);
double katzir2(const double, const double, const std::vector<double>&);
std::vector<double> katzirMeshOptimization(const std::vector<double>&,
	const std::vector<double>&,
	const std::vector<double>&,
	std::vector<double>&,
	const std::vector<double>&,
	const std::vector<double>&,
	const std::vector<double>&,
	const std::vector<bool>&);
double katzirResidual(const std::vector<double>&,
	const std::vector<double>&,
	const std::vector<double>&,
	const std::vector<double>&);
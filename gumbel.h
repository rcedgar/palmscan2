#pragma once

double Gumbel_PDF(double mu, double beta, double x);
double Gumbel_CDF(double mu, double beta, double x);
double Gumbel_logCDF(double mu, double beta, double x);
void FitGumbel(const vector<double> &xs, double dx,
  double &mu, double &beta);

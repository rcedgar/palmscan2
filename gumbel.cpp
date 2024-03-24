#include "myutils.h"

double Gumbel_PDF(double mu, double beta, double x)
	{
	double z = float(x - mu)/beta;
	double e_z = exp(-z);
	double PDF = (1.0/beta)*exp(-(z + e_z));
	return PDF;
	}

double Gumbel_CDF(double mu, double beta, double x)
	{
	double a = -(x - mu)/beta;
	double exp1 = exp(-a);
	double exp2 = exp(-exp1);
	return exp2;
	}

double Gumbel_logCDF(double mu, double beta, double x)
	{
	double a = -(x - mu)/beta;
	double exp1 = exp(-a);
	return 1.0/exp1;
	}

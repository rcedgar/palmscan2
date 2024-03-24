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

void cmd_fit_gumbel()
	{
	FILE *f = OpenStdioFile(g_Arg1);
	string Line;
	vector<string> Fields;
	vector<double> xs;
	vector<double> ys;
	double Mode = 0;
	double Maxy = 0;
	uint N = 0;
	double Sumy = 0;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) >= 2);
		double x = StrToFloat(Fields[0]);
		double y = StrToFloat(Fields[1]);
		Sumy += y;
		xs.push_back(x);
		ys.push_back(y);
		if (y > Maxy)
			{
			Maxy = y;
			Mode = x;
			}
		}
	const uint BinCount = SIZE(xs);
	asserta(SIZE(ys) == BinCount);

	vector<double> Freqs;
	double SumFreq = 0;
	for (uint Bin = 0; Bin < BinCount; ++Bin)
		{
		double y = ys[Bin];
		double Freq = y/Sumy;
		SumFreq += Freq;
		Freqs.push_back(Freq);
		}
	asserta(feq(SumFreq, 1.0));

	double Sumydx2 = 0;
	for (uint Bin = 0; Bin < BinCount; ++Bin)
		{
		double x = xs[Bin];
		double y = ys[Bin];
		double Freq = Freqs[Bin];
		double dx = (x - Mode);
		Sumydx2 += Freq*dx*dx;
		}
	double Var = Sumydx2/Sumy;
// Var = (PI^2/6)Beta^2
// Var*6/(PI^2) = Beta^2
// Beta = sqrt(6 Var/PI)
	const double PI = 3.1416;
	double Beta = sqrt(Var/(PI*PI));
	double SumRatio = 0;
	for (uint Bin = 0; Bin < BinCount; ++Bin)
		{
		double x = xs[Bin];
		double y = ys[Bin];
		double G = Gumbel_PDF(Mode, Beta, x);
		}

	ProgressLog("Mode = %.3g, var = %.3g, beta = %.3g\n",
	  Mode, Var, Beta);
	Log("x\ty\tG\tC\tlogC\n");
	for (uint Bin = 0; Bin < BinCount; ++Bin)
		{
		double x = xs[Bin];
		double y = ys[Bin];
		double G = Gumbel_PDF(Mode, Beta, x);
		double C = Gumbel_CDF(Mode, Beta, x);
		double logC = Gumbel_logCDF(Mode, Beta, x);
		Log("%.3g\t%.3g\t%.3g\t%.3g\t%.3g\n", x, y, G, C, logC);
		}
	}

#include "myutils.h"
#include "gumbel.h"

const double PI = 3.1415926536;
const double PI2 = PI*PI;

void FitGumbel(const vector<double> &xs, double dx,
  double &mu, double &beta)
	{
	mu = DBL_MAX;
	beta = DBL_MAX;

	const uint N = SIZE(xs);
	vector<double> smoothed_xs;
	vector<double> smoothed_ys;

// mu is modal x value
	double sumx = 0;
	double lox = xs[0];
	uint y = 0;
	uint maxy = 0;
	for (uint i = 1; i < N; ++i)
		{
		double x = xs[i];
		if (x - lox >= dx)
			{
			double meanx = (xs[i-1] + lox)/2; 
			smoothed_xs.push_back(meanx);
			smoothed_ys.push_back(y);
			if (y > maxy)
				{
				maxy = y;
				mu = meanx;
				}
			Log("smoothed\t%.3g\t%u\n", meanx, y);
			lox = x;
			y = 0;
			}
		else
			++y;
		}
	double meanx = (xs[N-1] + lox)/2; 
	smoothed_xs.push_back(meanx);
	smoothed_ys.push_back(y);
	const uint M = SIZE(smoothed_xs);
	asserta(SIZE(smoothed_ys) == M);

	vector<double> Ps;
	double SumP = 0;
	for (uint i = 0; i < M; ++i)
		{
		double y = smoothed_ys[i];
		double P = y/N;
		SumP += P;
		Ps.push_back(P);
		}
	asserta(feq(SumP, 1.0));

// Variance is expected value of (x - mu)^2
	double sumPdx2 = 0;
	for (uint i = 0; i < M; ++i)
		{
		double x = smoothed_xs[i];
		double P = Ps[i];
		double dx = (x - mu);
		sumPdx2 += P*dx*dx;
		}
	double Var = sumPdx2;

// Var = (PI^2/6)Beta^2
// Var*6/(PI^2) = Beta^2
// Beta = sqrt(6 Var/PI^2)
	beta = sqrt(Var/PI2);

	for (uint i = 0; i < M; ++i)
		{
		double x = smoothed_xs[i];
		double y = smoothed_ys[i];
		double G = Gumbel_PDF(mu, beta, x);
		Log("G\t%.3g\t%.3g\t%.3g\n", x, y, G);
		}
	}

void cmd_fit_gumbel()
	{
	FILE *f = OpenStdioFile(g_Arg1);
	string Line;
	vector<string> Fields;
	vector<double> xs;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		double x = StrToFloat(Fields[0]);
		if (x != 0)
			xs.push_back(x);
		}

	double mu, beta;
	FitGumbel(xs, 5, mu, beta);
	ProgressLog("mu %.3g beta %3g\n", mu, beta);
	}

#include "myutils.h"
#include "ppcprofileparams.h"

double Mean_L_AB = 71.4;
double StdDev_L_AB = 13.9;
double Mean_d_AB = 25.5;
double StdDev_d_AB = 2.1;

double Mean_L_AC = 109.3;
double StdDev_L_AC = 16.8;
double Mean_d_AC = 12.5;
double StdDev_d_AC = 3.9;

double Mean_L_BC = 38.9;
double StdDev_L_BC = 8.9;
double Mean_d_BC = 17.4;
double StdDev_d_BC = 2.6;

static double GetNormal(double Mu, double Sigma, double x)
	{
	static double TWOPI = (2.0*3.1415926535);
	static double FACTOR = 1.0/sqrt(TWOPI);

	double a = (x - Mu)/Sigma;
	double y = (FACTOR/Sigma)*exp(-0.5*a*a);
	return y;
	}

static double GetScoreNormal(double Mu, double Sigma, double x)
	{
	double f = GetNormal(Mu, Sigma, x);
	double Max = GetNormal(Mu, Sigma, Mu);
	double Ratio = f/Max;
	return Ratio;
	}

static void TestNormal1(double Mu, double Sigma)
	{
	Log("\n");
	Log("Mu=%.3g, Sigma=%.3g\n", Mu, Sigma);
	for (int i = -10; i <= 10; ++i)
		{
		double x = Mu + i*Sigma/3;
		double f = GetNormal(Mu, Sigma, x);
		double g = GetScoreNormal(Mu, Sigma, x);
		Log("  x=%8.3g  f=%8.3g  g=%8.3g\n", x, f, g);
		}
	}

static void TestNormal()
	{
	TestNormal1(1.0, 1.0);
	TestNormal1(1.0, 0.5);
	TestNormal1(10.0, 2.0);
	}

double GetPPCProfileScore(
  uint LAB, uint LAC, uint LBC,
  double dAB, double dAC, double dBC)
	{
	double Sum = 0.0;

#define X(x)	\
	Sum += GetScoreNormal(Mean_L_##x, StdDev_L_##x, L##x); \
	Sum += GetScoreNormal(Mean_d_##x, StdDev_d_##x, d##x); \

	X(AB)
	X(AC)
	X(BC)

	double Score = Sum/6;
	return Score;
	}

#include "myutils.h"

static double g_PhiTable[]
	{
	  -0.30103, // stddevs=0.00, y=5.000000e-01
	  -0.33708, // stddevs=0.10, y=4.601722e-01
	 -0.375986, // stddevs=0.20, y=4.207403e-01
	 -0.417836, // stddevs=0.30, y=3.820886e-01
	 -0.462712, // stddevs=0.40, y=3.445783e-01
	 -0.510692, // stddevs=0.50, y=3.085375e-01
	 -0.561848, // stddevs=0.60, y=2.742531e-01
	  -0.61625, // stddevs=0.70, y=2.419637e-01
	  -0.67396, // stddevs=0.80, y=2.118554e-01
	  -0.73504, // stddevs=0.90, y=1.840601e-01
	 -0.799546, // stddevs=1.00, y=1.586553e-01
	 -0.867529, // stddevs=1.10, y=1.356661e-01
	 -0.939039, // stddevs=1.20, y=1.150697e-01
	  -1.01412, // stddevs=1.30, y=9.680048e-02
	  -1.09282, // stddevs=1.40, y=8.075666e-02
	  -1.17518, // stddevs=1.50, y=6.680720e-02
	  -1.26123, // stddevs=1.60, y=5.479929e-02
	    -1.351, // stddevs=1.70, y=4.456546e-02
	  -1.44454, // stddevs=1.80, y=3.593032e-02
	  -1.54187, // stddevs=1.90, y=2.871656e-02
	  -1.64302, // stddevs=2.00, y=2.275013e-02
	  -1.74801, // stddevs=2.10, y=1.786442e-02
	  -1.85688, // stddevs=2.20, y=1.390345e-02
	  -1.96964, // stddevs=2.30, y=1.072411e-02
	  -2.08632, // stddevs=2.40, y=8.197536e-03
	  -2.20693, // stddevs=2.50, y=6.209665e-03
	   -2.3315, // stddevs=2.60, y=4.661188e-03
	  -2.46005, // stddevs=2.70, y=3.466974e-03
	  -2.59259, // stddevs=2.80, y=2.555130e-03
	  -2.72913, // stddevs=2.90, y=1.865813e-03
	   -2.8697, // stddevs=3.00, y=1.349898e-03
	   -3.0143, // stddevs=3.10, y=9.676032e-04
	  -3.16296, // stddevs=3.20, y=6.871379e-04
	  -3.31567, // stddevs=3.30, y=4.834241e-04
	  -3.47246, // stddevs=3.40, y=3.369293e-04
	  -3.63334, // stddevs=3.50, y=2.326291e-04
	  -3.79831, // stddevs=3.60, y=1.591086e-04
	  -3.96738, // stddevs=3.70, y=1.077997e-04
	  -4.14057, // stddevs=3.80, y=7.234804e-05
	  -4.31789, // stddevs=3.90, y=4.809634e-05
	  -4.49933, // stddevs=4.00, y=3.167124e-05
	  -4.68492, // stddevs=4.10, y=2.065751e-05
	  -4.87466, // stddevs=4.20, y=1.334575e-05
	  -5.06855, // stddevs=4.30, y=8.539905e-06
	   -5.2666, // stddevs=4.40, y=5.412544e-06
	  -5.46882, // stddevs=4.50, y=3.397673e-06
	  -5.67521, // stddevs=4.60, y=2.112455e-06
	  -5.88579, // stddevs=4.70, y=1.300807e-06
	  -6.10055, // stddevs=4.80, y=7.933282e-07
	   -6.3195, // stddevs=4.90, y=4.791833e-07
	  -6.54265, // stddevs=5.00, y=2.866516e-07
	  -6.76999, // stddevs=5.10, y=1.698267e-07
	  -7.00155, // stddevs=5.20, y=9.964426e-08
	  -7.23731, // stddevs=5.30, y=5.790134e-08
	  -7.47729, // stddevs=5.40, y=3.332045e-08
	  -7.72149, // stddevs=5.50, y=1.898956e-08
	   -7.9699, // stddevs=5.60, y=1.071759e-08
	  -8.22255, // stddevs=5.70, y=5.990371e-09
	  -8.47942, // stddevs=5.80, y=3.315746e-09
	  -8.74052, // stddevs=5.90, y=1.817508e-09
	  -9.00586, // stddevs=6.00, y=9.865876e-10
	  -9.27544, // stddevs=6.10, y=5.303423e-10
	  -9.54926, // stddevs=6.20, y=2.823158e-10
	  -9.82733, // stddevs=6.30, y=1.488228e-10
	  -10.1096, // stddevs=6.40, y=7.768847e-11
	  -10.3962, // stddevs=6.50, y=4.015999e-11
	   -10.687, // stddevs=6.60, y=2.055789e-11
	  -10.9821, // stddevs=6.70, y=1.042100e-11
	  -11.2814, // stddevs=6.80, y=5.230982e-12
	   -11.585, // stddevs=6.90, y=2.600142e-12
	  -11.8929, // stddevs=7.00, y=1.279810e-12
	   -12.205, // stddevs=7.10, y=6.237788e-13
	  -12.5214, // stddevs=7.20, y=3.010370e-13
	   -12.842, // stddevs=7.30, y=1.438849e-13
	  -13.1668, // stddevs=7.40, y=6.811218e-14
	   -13.496, // stddevs=7.50, y=3.191891e-14
	  -13.8291, // stddevs=7.60, y=1.482148e-14
	  -14.1657, // stddevs=7.70, y=6.827872e-15
	  -14.5074, // stddevs=7.80, y=3.108624e-15
	  -14.8577, // stddevs=7.90, y=1.387779e-15
	};

// LogDBSize = log_10(N) where N=number of queries
//   e.g. LogDBSize=6 for N=one million.
double GetExpValue(double Value, double Mean, double StdDev,
  double LogDBSize)
	{
	if (Value < Mean)
		return 1;
	const int N = int(sizeof(g_PhiTable)/sizeof(g_PhiTable[0]));

	double StdDevs = (Value - Mean)/StdDev;
	int Lo = int(floor(StdDevs*10));
	if (Lo+1 < N)
		{
		double ValueLo = g_PhiTable[Lo];
		double ValueHi = g_PhiTable[Lo+1];
		double ExpValue = ValueLo + (ValueHi - ValueLo)*(StdDevs - Lo/10.0);
		return ExpValue + LogDBSize;
		}
	double ExpValue = g_PhiTable[N-1];
	return ExpValue + LogDBSize;
	}

void cmd_test_expvalue()
	{
	for (uint i = 50; i < 100; ++i)
		{
		double Score = i/10.0 + 0.03;
		double ExpValue = GetExpValue(Score, 5, 0.5, 1);
		Log("%7.3f  %.6g\n", Score, ExpValue);
		}
	}

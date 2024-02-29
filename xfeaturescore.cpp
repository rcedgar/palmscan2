#include "myutils.h"
#include "xprof.h"

vector<vector<double> > XProf::g_BinLos(XFEATS);
vector<vector<double> > XProf::g_Scores(XFEATS);

double XProf::GetDiff(uint FeatureIndex, double Value1, double Value2)
	{
// Angles Ang_m2_p2=0, Ang_m3_p3=1
	if (FeatureIndex <= 1)
		{
		double d1 = fabs(Value1 - Value2);
		double d2 = fabs(fabs(Value1 - Value2) - 180);
		double Diff = min(d1, d2);
		asserta(Diff <= 90);
		return Diff;
		}
	return fabs(Value1 - Value2);
	}

double XProf::GetScore(char Amino1, char Amino2, const vector<uint> &Bins)
	{
	asserta(SIZE(Bins) == XFEATS);
	double Score = 0;
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		{
		uint Bin = Bins[FeatureIndex];
		Score += g_Scores[FeatureIndex][Bin];
		}
	return Score;
	}

uint XProf::GetFeatureBin(uint FeatureIndex, double Value)
	{
	asserta(FeatureIndex < XFEATS);
	for (uint Bin = 0; Bin < XBINS-1; ++Bin)
		{
		if (Value < g_BinLos[FeatureIndex][Bin+1])
			return Bin;
		}
	return XBINS-1;
	}

void XProf::FeatureScoreBin(const string &FeatureName, uint Bin,
  double BinLo, double Score)
	{
	uint FeatureIndex = XProf::GetFeatureIndex(FeatureName);
	asserta(FeatureIndex < XFEATS);
	asserta(Bin < XBINS);
	asserta(g_BinLos[FeatureIndex][Bin] == DBL_MAX);
	asserta(g_Scores[FeatureIndex][Bin] == DBL_MAX);
	g_BinLos[FeatureIndex][Bin] = BinLo;
	g_Scores[FeatureIndex][Bin] = Score;
	}

void XProf::InitScoreTable()
	{
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		{
		g_BinLos[FeatureIndex].resize(XBINS, DBL_MAX);
		g_Scores[FeatureIndex].resize(XBINS, DBL_MAX);
		}

	FeatureScoreBin("Ang_m2_p2", 0, 0, 1.38);
	FeatureScoreBin("Ang_m2_p2", 1, 1, 1.23);
	FeatureScoreBin("Ang_m2_p2", 2, 2, 1.16);
	FeatureScoreBin("Ang_m2_p2", 3, 3.68, 1.12);
	FeatureScoreBin("Ang_m2_p2", 4, 5, 0.931);
	FeatureScoreBin("Ang_m2_p2", 5, 7.1, 0.825);
	FeatureScoreBin("Ang_m2_p2", 6, 10, 0.528);
	FeatureScoreBin("Ang_m2_p2", 7, 14.2, 0.235);
	FeatureScoreBin("Ang_m2_p2", 8, 20.5, -0.115);
	FeatureScoreBin("Ang_m2_p2", 9, 30.6, -1.17);

	FeatureScoreBin("Ang_m3_p3", 0, 0, 1.22);
	FeatureScoreBin("Ang_m3_p3", 1, 1.2, 1.01);
	FeatureScoreBin("Ang_m3_p3", 2, 2.7, 1.01);
	FeatureScoreBin("Ang_m3_p3", 3, 4.3, 0.875);
	FeatureScoreBin("Ang_m3_p3", 4, 6.2, 0.724);
	FeatureScoreBin("Ang_m3_p3", 5, 8.6, 0.523);
	FeatureScoreBin("Ang_m3_p3", 6, 11.7, 0.293);
	FeatureScoreBin("Ang_m3_p3", 7, 16.1, -0.0062);
	FeatureScoreBin("Ang_m3_p3", 8, 22.8, -0.317);
	FeatureScoreBin("Ang_m3_p3", 9, 33, -0.985);

	FeatureScoreBin("ED_p4", 0, 0, 1.37);
	FeatureScoreBin("ED_p4", 1, 0.07, 1.33);
	FeatureScoreBin("ED_p4", 2, 0.13, 1.2);
	FeatureScoreBin("ED_p4", 3, 0.2, 1.16);
	FeatureScoreBin("ED_p4", 4, 0.3, 1.08);
	FeatureScoreBin("ED_p4", 5, 0.44, 0.88);
	FeatureScoreBin("ED_p4", 6, 0.6, 0.623);
	FeatureScoreBin("ED_p4", 7, 0.9, 0.267);
	FeatureScoreBin("ED_p4", 8, 1.37, -0.143);
	FeatureScoreBin("ED_p4", 9, 2.19, -1.2);

	FeatureScoreBin("ED_m4", 0, 0, 1.37);
	FeatureScoreBin("ED_m4", 1, 0.07, 1.35);
	FeatureScoreBin("ED_m4", 2, 0.13, 1.21);
	FeatureScoreBin("ED_m4", 3, 0.2, 1.16);
	FeatureScoreBin("ED_m4", 4, 0.3, 1.09);
	FeatureScoreBin("ED_m4", 5, 0.44, 0.885);
	FeatureScoreBin("ED_m4", 6, 0.6, 0.626);
	FeatureScoreBin("ED_m4", 7, 0.9, 0.27);
	FeatureScoreBin("ED_m4", 8, 1.38, -0.153);
	FeatureScoreBin("ED_m4", 9, 2.2, -1.2);

	FeatureScoreBin("NU", 0, 0, 1.08);
	FeatureScoreBin("NU", 1, 1, 0.941);
	FeatureScoreBin("NU", 2, 2, 0.788);
	FeatureScoreBin("NU", 3, 3, 0.6);
	FeatureScoreBin("NU", 4, 4, 0.394);
	FeatureScoreBin("NU", 5, 5, 0.181);
	FeatureScoreBin("NU", 6, 6, -0.0334);
	FeatureScoreBin("NU", 7, 7, -0.242);
	FeatureScoreBin("NU", 8, 8, -0.436);
	FeatureScoreBin("NU", 9, 9, 0);

	FeatureScoreBin("ND", 0, 0, 0.776);
	FeatureScoreBin("ND", 1, 1, 0.692);
	FeatureScoreBin("ND", 2, 2, 0.538);
	FeatureScoreBin("ND", 3, 3, 0.35);
	FeatureScoreBin("ND", 4, 4, 0.119);
	FeatureScoreBin("ND", 5, 5, -0.103);
	FeatureScoreBin("ND", 6, 6, -0.33);
	FeatureScoreBin("ND", 7, 7, -0.553);
	FeatureScoreBin("ND", 8, 8, -0.767);
	FeatureScoreBin("ND", 9, 9, 0);

	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		{
		for (uint Bin = 0; Bin < XBINS; ++Bin)
			{
			asserta(g_BinLos[FeatureIndex][Bin] != DBL_MAX);
			asserta(g_Scores[FeatureIndex][Bin] != DBL_MAX);
			}
		}
	}

#include "myutils.h"
#include "alpha.h"
#include "xprof.h"
#include "xbinner.h"

static char g_Aminos[20];
static vector<vector<double> > g_FeatureValuesVec;

void XBinner::DefineCentroid(uint Idx, char aa,
	double v0, double v1, double v2,
	double v3, double v4, double v5)
	{
	if (g_FeatureValuesVec.empty())
		{
		g_FeatureValuesVec.resize(20);
		for (uint i = 0; i < 20; ++i)
			g_FeatureValuesVec[i].resize(XFEATS);
		}
	asserta(Idx < 20);
	asserta(g_Aminos[Idx] == 0);

	g_Aminos[Idx] = aa;
	g_FeatureValuesVec[Idx][0] = v0;
	g_FeatureValuesVec[Idx][1] = v1;
	g_FeatureValuesVec[Idx][2] = v2;
	g_FeatureValuesVec[Idx][3] = v3;
	g_FeatureValuesVec[Idx][4] = v4;
	g_FeatureValuesVec[Idx][5] = v5;
	}

void XBinner::InitCentroids()
	{
DefineCentroid(0, 'W', 115, 33, 6.47, 8.3, 3, 20); // 0.7078
DefineCentroid(1, 'V', 5.31, 13.6, 13.6, 13.1, 20, 15); // 0.7078
DefineCentroid(2, 'K', 103, 89.2, 9.3, 5.74, 3, 14); // 0.7078
DefineCentroid(3, 'T', 68.7, 113, 11.4, 10.3, 12, 17); // 0.7078
DefineCentroid(4, 'F', 49.7, 37.8, 12.9, 11.7, 7, 19); // 0.7078
DefineCentroid(5, 'A', 45.9, 80.8, 8, 9.09, 3, 18); // 0.7078
DefineCentroid(6, 'T', 55.3, 64.9, 13.1, 5.28, 19, 10); // 0.7078
DefineCentroid(7, 'A', 25.6, 39.4, 6.27, 13.9, 19, 23); // 0.7078
DefineCentroid(8, 'Y', 127, 98.4, 6.16, 11.4, 15, 13); // 0.7078
DefineCentroid(9, 'K', 67, 79.8, 5.88, 6.67, 17, 15); // 0.7078
DefineCentroid(10, 'R', 6.75, 29.3, 10.7, 11.5, 5, 22); // 0.7078
DefineCentroid(11, 'R', 97.2, 116, 10.4, 14, 0, 8); // 0.7078
DefineCentroid(12, 'T', 92.2, 141, 12.4, 6.83, 0, 15); // 0.7078
DefineCentroid(13, 'S', 65.6, 53.8, 10.9, 7.46, 8, 9); // 0.7078
DefineCentroid(14, 'L', 132, 110, 10.3, 5.92, 30, 19); // 0.7078
DefineCentroid(15, 'L', 75.3, 79.1, 11.3, 13, 10, 13); // 0.7078
DefineCentroid(16, 'P', 24.7, 57.7, 13.4, 8.32, 5, 17); // 0.7078
DefineCentroid(17, 'A', 125, 60.5, 4.84, 6.04, 17, 20); // 0.7078
DefineCentroid(18, 'N', 129, 142, 10.1, 13.2, 9, 14); // 0.7078
DefineCentroid(19, 'L', 107, 23.3, 5.99, 5.95, 33, 12); // 0.7078
	}

uint XBinner::GetLetter(uint AminoLetter, const vector<double> &FeatureValues) const
	{
	//return AminoLetter;
	if (AminoLetter >= 20)
		return 0;

	char AminoChar = g_LetterToCharAmino[AminoLetter];
	double BestScore = -999;
	uint BestLetter = UINT_MAX;
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		const vector<double> &v = g_FeatureValuesVec[Letter];
		double Score =
		  XProf::GetScore2(AminoChar, g_Aminos[Letter], FeatureValues, v);
		if (Score > BestScore)
			{
			BestScore = Score;
			BestLetter = Letter;
			}
		}
	asserta(BestLetter != UINT_MAX);
	return BestLetter;
	}

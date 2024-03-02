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
DefineCentroid(0, 'E', 110, 25.6, 6.03, 5.6, 17, 12);
DefineCentroid(1, 'N', 83.1, 62.8, 9.57, 11.7, 0, 14);
DefineCentroid(2, 'L', 72.9, 48.9, 5.6, 8.06, 19, 4);
DefineCentroid(3, 'L', 13.9, 22, 11.6, 0, 1, 19);
DefineCentroid(4, 'K', 125, 109, 6.2, 13.1, 1, 13);
DefineCentroid(5, 'T', 34.4, 74.9, 6.27, 9.21, 2, 16);
DefineCentroid(6, 'L', 18.5, 10.6, 5.74, 12.9, 21, 5);
DefineCentroid(7, 'Q', 117, 111, 11.5, 6.58, 3, 20);
DefineCentroid(8, 'R', 6.25, 26.3, 13, 8.64, 2, 26);
DefineCentroid(9, 'K', 42.9, 39.7, 28.3, 13.5, 1, 16);
DefineCentroid(10, 'G', 63.1, 92.2, 14, 5.95, 3, 12);
DefineCentroid(11, 'G', 78, 111, 6.42, 10.2, 0, 18);
DefineCentroid(12, 'P', 83.3, 75.5, 8.64, 5.31, 16, 7);
DefineCentroid(13, 'G', 64.3, 38.2, 12.8, 12.4, 22, 19);
DefineCentroid(14, 'L', 108, 145, 10.1, 6.15, 13, 11);
DefineCentroid(15, 'V', 27.3, 59, 10.6, 13.8, 19, 16);
DefineCentroid(16, 'Y', 6.12, 21.5, 8.8, 13.6, 8, 12);
DefineCentroid(17, 'S', 99.5, 144, 11.6, 9.04, 1, 15);
DefineCentroid(18, 'V', 29.5, 70.3, 13.4, 6.58, 9, 8);
DefineCentroid(19, 'R', 32.2, 0, 0, 5.43, 1, 12);
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

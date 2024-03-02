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
DefineCentroid(0, 'K', 109, 28, 6.23, 6.07, 15, 9);
DefineCentroid(1, 'V', 15.3, 5.99, 13.8, 13, 12, 15);
DefineCentroid(2, 'Q', 43.1, 41.8, 12.5, 12.8, 0, 19);
DefineCentroid(3, 'S', 99.2, 63.8, 6.11, 7.95, 3, 17);
DefineCentroid(4, 'L', 30.6, 28.6, 12.8, 7.36, 20, 9);
DefineCentroid(5, 'L', 114, 53.4, 9.56, 6.36, 19, 23);
DefineCentroid(6, 'A', 3.82, 7.06, 10.9, 12.3, 20, 7);
DefineCentroid(7, 'D', 8.58, 43.8, 12.3, 9.39, 16, 22);
DefineCentroid(8, 'D', 84.2, 111, 11, 8.21, 9, 16);
DefineCentroid(9, 'A', 115, 119, 9.57, 11.8, 0, 9);
DefineCentroid(10, 'D', 42.5, 95.4, 5.81, 11.4, 9, 13);
DefineCentroid(11, 'S', 29.4, 22.6, 8.21, 10.9, 17, 15);
DefineCentroid(12, 'G', 55, 85.2, 12.8, 4.75, 5, 20);
DefineCentroid(13, 'H', 69.5, 62.3, 8.98, 5.76, 5, 13);
DefineCentroid(14, 'G', 89.8, 102, 12.7, 11.2, 16, 25);
DefineCentroid(15, 'H', 20.9, 68.1, 6.52, 13.3, 18, 20);
DefineCentroid(16, 'V', 96.6, 83, 8.53, 12.9, 14, 12);
DefineCentroid(17, 'L', 65.3, 64.8, 6, 10.4, 20, 22);
DefineCentroid(18, 'E', 116, 67.3, 10.5, 11, 5, 17);
DefineCentroid(19, 'E', 52.7, 71.5, 9.84, 9.63, 13, 9);
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

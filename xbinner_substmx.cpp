#include "myutils.h"
#include "xprof.h"
#include "alpha.h"
#include "xbinner.h"
#include "x2data.h"

/***
palmscan2 \
  -xbinner_substmx d:/int/dave_grant/scop40/aligned_xfeatvecs_0.6_0.8.tsv \
  -log xbinner_substmx.log

Tests minimal XBinner::GetLetter() which returns aa letter.

Input is pairs of XProf feature vectors aligned by SCOP40 TM 0.6 to 0.8
=======================================================================

Output is substmx and entropy in log
====================================

static double g_AminoSubstMx[20][20] = {
//      A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
{  0.7892,  0.1351, -0.2798, -0.1031, -0.2414, -0.0349, -0.2369, -0.1676, -0.1370, -0.1271, -0.0086, -0.2398, -0.1653, -0.1084, -0.1764,  0.0487, -0.0833,  0.0265, -0.3102, -0.2366}, // A
...
{ -0.2366, -0.2110, -0.5238, -0.4107,  0.7346, -0.5309,  0.3342, -0.0746, -0.3487,  0.0090,  0.0489, -0.2360, -0.3784, -0.2277, -0.1899, -0.2676, -0.2150, -0.1044,  0.6562,  1.4538}, // Y
};
Entropy 5.63, expected score 0.147
***/

void LogCountsMx(const vector<vector<uint> > &CountMx)
	{
	assert(SIZE(CountMx) == 20);
	Log("static double Counts[20][20] = {\n");
	Log("//      ");
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		Log("        %c", c);
		}
	Log("\n");
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		Log("/* %c */ {", c);
		assert(SIZE(CountMx[i]) == 20);
		for (uint j = 0; j < 20; ++j)
			{
			uint Count = CountMx[i][j];
			Log(" %7u", Count);
			if (j != 19)
				Log(",");
			}
		Log("}, // %c\n", c);
		}
	Log("};\n");
	}

void LogScoreMx(const vector<vector<double> > &ScoreMx)
	{
	assert(SIZE(ScoreMx) == 20);
	Log("static double g_AminoSubstMx[20][20] = {\n");
	Log("//      ");
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		Log("        %c", c);
		}
	Log("\n");
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		Log("/* %c */ {", c);
		assert(SIZE(ScoreMx[i]) == 20);
		for (uint j = 0; j < 20; ++j)
			{
			double Score = ScoreMx[i][j];
			Log(" %7.4f", Score);
			if (j != 19)
				Log(",");
			}
		Log("}, // %c\n", c);
		}
	Log("};\n");
	}

static void GetHMMStrings(const vector<vector<double> > &FreqMx,
  vector<string> &Lines)
	{
	Lines.clear();
#define ADD_STR(s)	Lines.push_back(s);

ADD_STR("HMM	aa")
ADD_STR("T.START_M	0.6")
ADD_STR("T.START_IS	0.02")
ADD_STR("T.START_IL	0.18")
ADD_STR("T.M_M	0.96")
ADD_STR("T.M_IS	0.012")
ADD_STR("T.M_IL	0.008")
ADD_STR("T.IS_IS	0.35")
ADD_STR("T.IS_M	0.65")
ADD_STR("T.IL_IL	0.90")
ADD_STR("T.IL_M	0.10")
#undef ADD_STR

//ADD_STR("E.AA	0.023731")
//ADD_STR("E.CA	0.0014551")
//ADD_STR("E.CC	0.010135")
	for (uint Letter1 = 0; Letter1 < 20; ++Letter1)
		{
		char Char1 = g_LetterToCharAmino[Letter1];
		for (uint Letter2 = 0; Letter2 <= Letter1; ++Letter2)
			{
			char Char2 = g_LetterToCharAmino[Letter2];
			double P = FreqMx[Letter1][Letter2];

			string s;
			Ps(s, "E.%c%c	%.4g", Char1, Char2, P);
			Lines.push_back(s);
			}
		}
	}

void cmd_xbinner_substmx()
	{
	XBinner XB;
	XBinner::InitCentroids();

	XProf::InitScoreTable();

	X2Data X2;
	X2.FromTsv(opt_xbinner_substmx);

	vector<double> Freqs;
	vector<vector<double> > FreqMx;
	XB.GetFreqs(X2, Freqs, FreqMx);

	vector<vector<double> > ScoreMx;
	double ExpScore = XB.GetLogOddsMx(Freqs, FreqMx, ScoreMx);

	LogScoreMx(ScoreMx);
	ProgressLog("Expected score %.3f\n", ExpScore);

	vector<string> Lines;
	GetHMMStrings(FreqMx, Lines);
	for (uint i = 0; i < SIZE(Lines); ++i)
		Log("%s\n", Lines[i].c_str());
	}

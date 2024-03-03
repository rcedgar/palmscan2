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

void ReadAlignedValueVecs(const string &FileName,
  vector<char> &AminosQ,
  vector<char> &AminosR,
  vector<vector<double> > &ValuesQVec,
  vector<vector<double> > &ValuesRVec)
	{
	AminosQ.clear();
	AminosR.clear();
	ValuesQVec.clear();
	ValuesRVec.clear();

	FILE *f = OpenStdioFile(FileName);
	string HdrLine;
	vector<string> HdrFields;
	bool Ok = ReadLineStdioFile(f, HdrLine);
	asserta(Ok);
	Split(HdrLine, HdrFields, '\t');
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields[0]) == 1);
		asserta(SIZE(Fields[XFEATS+1]) == 1);
		char AminoQ = Fields[0][0];
		char AminoR = Fields[XFEATS+1][0];
		vector<double> ValuesQ;
		vector<double> ValuesR;
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			double ValueQ = StrToFloat(Fields[FeatureIndex+1]);
			double ValueR = StrToFloat(Fields[XFEATS+FeatureIndex+2]);
			ValuesQ.push_back(ValueQ);
			ValuesR.push_back(ValueR);
			}
		AminosQ.push_back(AminoQ);
		AminosR.push_back(AminoR);
		ValuesQVec.push_back(ValuesQ);
		ValuesRVec.push_back(ValuesR);
		}
	}

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

void cmd_xbinner_substmx()
	{
	XBinner XB;
	XBinner::InitCentroids();

	XProf::InitScoreTable();

	X2Data X2;
	X2.FromTsv(opt_xbinner_substmx);
	const vector<vector<double> > ValuesQVec = X2.m_FeatureValuesVec1;
	const vector<vector<double> > ValuesRVec = X2.m_FeatureValuesVec2;
	const vector<char> AminoQs = X2.m_Aminos1;
	const vector<char> AminoRs = X2.m_Aminos2;

	vector<double> Freqs;
	vector<vector<double> > FreqMx;
	XB.GetFreqs(X2, Freqs, FreqMx);

	vector<vector<double> > ScoreMx;
	double ExpScore = XB.GetLogOddsMx(Freqs, FreqMx, ScoreMx);

	LogScoreMx(ScoreMx);
	ProgressLog("Expected score %.3f\n", ExpScore);
	}

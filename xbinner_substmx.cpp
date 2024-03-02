#include "myutils.h"
#include "xprof.h"
#include "alpha.h"
#include "xbinner.h"
#undef zero
#include <chrono>

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

void GetFreqs(XBinner &XB, 
  const vector<char> &AminoQs,
  const vector<char> &AminoRs,
  const vector<vector<double> > &ValuesQVec,
  const vector<vector<double> > &ValuesRVec,
  vector<double> &Freqs,
  vector<vector<double> > &FreqMx)
	{
	Freqs.clear();
	FreqMx.clear();

	const uint PosPairCount = SIZE(ValuesQVec);
	assert(SIZE(ValuesRVec) == PosPairCount);

	uint LetterPairCount = 0;
	vector<uint> CountVec(20, 1);
	vector<vector<uint> > CountMx(20);
	for (uint i = 0; i < 20; ++i)
		CountMx[i].resize(20, 1);

	for (uint PairIndex = 0; PairIndex < PosPairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PosPairCount, "freqs");
		char AminoQ = AminoQs[PairIndex];
		char AminoR = AminoRs[PairIndex];

		uint AminoLetterQ = g_CharToLetterAmino[AminoQ];
		uint AminoLetterR = g_CharToLetterAmino[AminoR];

		const vector<double> &ValuesQ = ValuesQVec[PairIndex];
		const vector<double> &ValuesR = ValuesRVec[PairIndex];

		uint iq = XB.GetLetter(AminoLetterQ, ValuesQ);
		uint ir = XB.GetLetter(AminoLetterR, ValuesR);
		if (iq >= 20 || ir >= 20)
			continue;

		LetterPairCount += 2;
		CountVec[iq] += 1;
		CountVec[ir] += 1;
		CountMx[iq][ir] += 1;
		CountMx[ir][iq] += 1;
		}

	double SumFreq = 0;
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		double Freq = double(CountVec[i])/(LetterPairCount);
		SumFreq += Freq;
		Freqs.push_back(Freq);
		Log("Freqs[%c] = %8.6f\n", c, Freq);
		}
	LogCountsMx(CountMx);

	assert(feq(SumFreq, 1.0));

	double SumFreq2 = 0;
	FreqMx.resize(20);
	for (uint i = 0; i < 20; ++i)
		{
		for (uint j = 0; j < 20; ++j)
			{
			uint n = CountMx[i][j];
			double Freq2 = double(n)/double(LetterPairCount);
			FreqMx[i].push_back(Freq2);
			SumFreq2 += Freq2;
			}
		}
	assert(feq(SumFreq2, 1.0));
	}

void GetLogOddsMx(const vector<double> &Freqs,
  const vector<vector<double> > &FreqMx,
  vector<vector<double> > &ScoreMx,
  double &H, double &RH)
	{
	assert(SIZE(Freqs) == 20);
	assert(SIZE(FreqMx) == 20);
	ScoreMx.clear();
	ScoreMx.resize(20);

	H = 0;
	RH = 0;
	for (uint i = 0; i < 20; ++i)
		{
		double Freqi = Freqs[i];
		for (uint j = 0; j < 20; ++j)
			{
			double Freqj = Freqs[j];
			double ObsFreq = FreqMx[i][j];
			assert(feq(FreqMx[i][j], FreqMx[j][i]));

			double ExpFreq = double(Freqi*Freqj);
			double Score = log(ObsFreq/ExpFreq);
			H -= ObsFreq*log(ObsFreq);
			RH += ObsFreq*Score;
			ScoreMx[i].push_back(Score);
			}
		}
	}

void cmd_xbinner_substmx()
	{
	XBinner XB;
	XBinner::InitCentroids();

	XProf::InitScoreTable();

	vector<vector<double> > ValuesQVec;
	vector<vector<double> > ValuesRVec;
	vector<char> AminoQs;
	vector<char> AminoRs;
	Progress("Reading pairs...");
	ReadAlignedValueVecs(opt_xbinner_substmx, AminoQs, AminoRs,
	  ValuesQVec, ValuesRVec);
	Progress(" done.\n");

	auto start_time = std::chrono::steady_clock::now();
	vector<double> Freqs;
	vector<vector<double> > FreqMx;
	GetFreqs(XB, AminoQs, AminoRs, ValuesQVec, ValuesRVec, Freqs, FreqMx);

	double H, RH;
	vector<vector<double> > ScoreMx;
	GetLogOddsMx(Freqs, FreqMx, ScoreMx, H, RH);
	auto end_time = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
	double elapsedd = std::chrono::duration<double>(elapsed).count();

	LogScoreMx(ScoreMx);
	ProgressLog("Entropy %.3g, expected score %.3g, calc time %.3g ms\n",
	  H, RH, elapsedd);
	}

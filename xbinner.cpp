#include "myutils.h"
#include "alpha.h"
#include "xprof.h"
#include "xbinner.h"

static char g_Aminos[20];
static vector<vector<double> > g_FeatureValuesVec;

void XBinner::SetCentroid(uint Idx, char aa, const vector<double> &Values)
	{
	asserta(Idx < 20);
	if (g_FeatureValuesVec.empty())
		{
		g_FeatureValuesVec.resize(20);
		for (uint i = 0; i < 20; ++i)
			g_FeatureValuesVec[i].resize(XFEATS);
		}
	g_Aminos[Idx] = aa;
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		g_FeatureValuesVec[Idx][FeatureIndex] = Values[FeatureIndex];
	}

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

void XBinner::LogCentroids()
	{
	Log("XBinner::LogCentroids\n");
	for (uint i = 0; i < 20; ++i)
		{
		Log("%c", g_Aminos[i]);
		for (uint j = 0; j < XFEATS; ++j)
			Log(" %10.3g", g_FeatureValuesVec[i][j]);
		Log("\n");
		}
	}

void XBinner::InitCentroids()
	{
DefineCentroid(0, 'Y', 105, 21.8, 6.02, 6.29, 18, 28); // 0.7065
DefineCentroid(1, 'L', 2.56, 21.6, 13.6, 12.3, 13, 24); // 0.7065
DefineCentroid(2, 'H', 114, 80.2, 6.42, 13.4, 13, 14); // 0.7065
DefineCentroid(3, 'K', 87.1, 96, 12.7, 11.9, 9, 11); // 0.7065
DefineCentroid(4, 'G', 32.7, 54.4, 11.4, 7.66, 0, 17); // 0.7065
DefineCentroid(5, 'A', 23.1, 73.2, 6.56, 10.8, 6, 18); // 0.7065
DefineCentroid(6, 'T', 34.7, 60, 8.99, 13, 19, 18); // 0.7065
DefineCentroid(7, 'V', 69.1, 29, 9.35, 8.26, 22, 18); // 0.7065
DefineCentroid(8, 'G', 116, 145, 12.2, 9.55, 20, 13); // 0.7065
DefineCentroid(9, 'L', 115, 40.7, 8.33, 6.88, 2, 18); // 0.7065
DefineCentroid(10, 'E', 90.5, 119, 5.99, 8.92, 1, 15); // 0.7065
DefineCentroid(11, 'H', 72.5, 58.4, 13.1, 10.5, 3, 21); // 0.7065
DefineCentroid(12, 'K', 133, 112, 10.6, 6.58, 2, 12); // 0.7065
DefineCentroid(13, 'Y', 32.9, 66, 9.69, 5.03, 13, 12); // 0.7065
DefineCentroid(14, 'V', 98, 67, 13.2, 6.59, 21, 23); // 0.7065
DefineCentroid(15, 'Y', 20.2, 29.2, 13.1, 8.98, 25, 16); // 0.7065
DefineCentroid(16, 'G', 119, 107, 8.25, 11.2, 0, 19); // 0.7065
DefineCentroid(17, 'M', 60.1, 20.6, 10.9, 11.5, 18, 8); // 0.7065
DefineCentroid(18, 'G', 84.8, 129, 10.7, 6.22, 8, 16); // 0.7065
DefineCentroid(19, 'F', 15.8, 13, 6.25, 13.1, 23, 11); // 0.7065
	}

uint XBinner::GetLetter(char AminoChar,
  const vector<double> &FeatureValues) const
	{
	extern bool g_Trace;
	uint AminoLetter = g_CharToLetterAmino[AminoChar];
	if (AminoLetter >= 20)
		return 0;

	double BestScore = -999;
	uint BestLetter = UINT_MAX;
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		const vector<double> &v = g_FeatureValuesVec[Letter];
		char AminoChar = g_Aminos[Letter];
		uint AminoLetter2 = g_CharToLetterAmino[AminoChar];
		double Score =
		  XProf::GetScore_Letters2(AminoLetter, AminoLetter2, FeatureValues, v);
		if (g_Trace && Letter == 0)
			{
			Log("XBinner::GetLetter(%c)", AminoChar);
			Log(" [%2u]=%.3g\n", Letter, Score);
			Log("  Ref=");
			for (uint k = 0; k < XFEATS; ++k)
				Log(" %.3g", v[k]);
			Log("\n");
			Log("  Q=");
			for (uint k = 0; k < XFEATS; ++k)
				Log(" %.3g", FeatureValues[k]);
			Log("\n");
			}
		if (Score > BestScore)
			{
			BestScore = Score;
			BestLetter = Letter;
			}
		}
	asserta(BestLetter != UINT_MAX);
	if (g_Trace)
		{
		Log(">>>  XBinner::GetLetter(%c", AminoChar);
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			Log(" %.3g", FeatureValues[FeatureIndex]);
		Log(") = %u\n", BestLetter);
		}
	return BestLetter;
	}

bool g_Trace;

void XBinner::GetFreqs(
  const X2Data &X2,
  vector<double> &Freqs,
  vector<vector<double> > &FreqMx) const
	{
	Freqs.clear();
	FreqMx.clear();

	const vector<char> &AminoQs = X2.m_Aminos1;
	const vector<char> &AminoRs = X2.m_Aminos2;
	const vector<vector<double> > &ValuesQVec = X2.m_FeatureValuesVec1;
	const vector<vector<double> > &ValuesRVec = X2.m_FeatureValuesVec2;

	const uint PosPairCount = SIZE(ValuesQVec);
	assert(SIZE(ValuesRVec) == PosPairCount);

	uint LetterPairCount = 0;
	vector<uint> CountVec(20, 1);
	vector<vector<uint> > CountMx(20);
	for (uint i = 0; i < 20; ++i)
		CountMx[i].resize(20, 1);

	uint Counter = 0;
//#pragma omp parallel for
	for (int iPairIndex = 0; iPairIndex < (int) PosPairCount; ++iPairIndex)
		{
		g_Trace = (Counter < 10);
		uint PairIndex = (uint) iPairIndex;
		char AminoQ = AminoQs[PairIndex];
		char AminoR = AminoRs[PairIndex];

		//uint AminoLetterQ = g_CharToLetterAmino[AminoQ];
		//uint AminoLetterR = g_CharToLetterAmino[AminoR];

		const vector<double> &ValuesQ = ValuesQVec[PairIndex];
		const vector<double> &ValuesR = ValuesRVec[PairIndex];

		uint iq = GetLetter(AminoQ, ValuesQ);
		uint ir = GetLetter(AminoR, ValuesR);
		if (iq >= 20 || ir >= 20)
			continue;
#pragma omp critical
		{
		if (Counter < 100)
			Log("TRACE %u %c %c %u %u\n", Counter, AminoQ, AminoR, iq, ir);
		ProgressStep(Counter++, PosPairCount, "freqs");
		LetterPairCount += 2;
		CountVec[iq] += 1;
		CountVec[ir] += 1;
		CountMx[iq][ir] += 1;
		CountMx[ir][iq] += 1;
		}
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
	//LogCountsMx(CountMx);

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

double XBinner::GetLogOddsMx(const vector<double> &Freqs,
  const vector<vector<double> > &FreqMx,
  vector<vector<double> > &ScoreMx)
	{
	assert(SIZE(Freqs) == 20);
	assert(SIZE(FreqMx) == 20);
	ScoreMx.clear();
	ScoreMx.resize(20);

	double ExpScore = 0;
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
			ExpScore += ObsFreq*Score;
			ScoreMx[i].push_back(Score);
			}
		}
	return ExpScore;
	}

uint XBinnerC::GetLetter(char AminoChar,
  const vector<double> &FeatureValues) const
	{
	asserta(m_ptrX2 != 0);
	asserta(SIZE(m_Idxs) == 20);
	uint AminoLetter = g_CharToLetterAmino[AminoChar];
	if (AminoLetter >= 20)
		return 0;

	extern bool g_Trace;
	double BestScore = -999;
	uint BestLetter = UINT_MAX;
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		uint Idx = m_Idxs[Letter];
		const vector<double> &v = m_ptrX2->m_FeatureValuesVec[Idx];
		char AminoChar2 = m_ptrX2->m_Aminos[Idx];
		uint AminoLetter2 = g_CharToLetterAmino[AminoChar2];
		double Score =
		  XProf::GetScore_Letters2(AminoLetter, AminoLetter2, FeatureValues, v);
		if (g_Trace && Letter == 0)
			{
			Log("XBinnerC::GetLetter(%c)", AminoChar);
			Log(" [%2u]=%.3g\n", Letter, Score);
			Log("  Ref=");
			for (uint k = 0; k < XFEATS; ++k)
				Log(" %.3g", v[k]);
			Log("\n");
			Log("  Q=");
			for (uint k = 0; k < XFEATS; ++k)
				Log(" %.3g", FeatureValues[k]);
			Log("\n");
			}
		if (Score > BestScore)
			{
			BestScore = Score;
			BestLetter = Letter;
			}
		}
	asserta(BestLetter != UINT_MAX);
	if (g_Trace)
		{
		Log(">>>  XBinnerC::GetLetter(%c", AminoChar);
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			Log(" %.3g", FeatureValues[FeatureIndex]);
		Log(") = %u\n", BestLetter);
		}
	return BestLetter;
	}

void XBinnerC::LogMyCentroids() const
	{
	Log("XBinnerC::LogMyCentroids\n");
	for (uint i = 0; i < 20; ++i)
		{
		uint Idx = m_Idxs[i];
		const vector<double> &v = m_ptrX2->m_FeatureValuesVec[Idx];
		Log("%c", g_Aminos[i]);
		for (uint j = 0; j < XFEATS; ++j)
			Log(" %10.3g", g_FeatureValuesVec[i][j]);
		Log("\n");
		}
	}


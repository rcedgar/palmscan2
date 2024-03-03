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
DefineCentroid(0, 'L', 108, 40.4, 7.99, 6.15, 26, 15);
DefineCentroid(1, 'V', 16.8, 22.4, 12.8, 11.4, 18, 31);
DefineCentroid(2, 'K', 38.8, 29.4, 7.13, 11.8, 15, 20);
DefineCentroid(3, 'G', 62, 108, 11, 5.93, 0, 14);
DefineCentroid(4, 'S', 72.7, 61.5, 13.6, 12.9, 3, 19);
DefineCentroid(5, 'E', 98.4, 116, 6.13, 11.1, 3, 22);
DefineCentroid(6, 'H', 66.3, 46.4, 10.2, 7.86, 19, 22);
DefineCentroid(7, 'S', 19.6, 47.2, 8.5, 13.7, 13, 14);
DefineCentroid(8, 'L', 129, 103, 9.73, 8.33, 1, 17);
DefineCentroid(9, 'D', 128, 145, 12, 13.1, 6, 13);
DefineCentroid(10, 'V', 112, 77.9, 6.43, 9.18, 21, 24);
DefineCentroid(11, 'N', 69.1, 73.3, 5.65, 6.26, 6, 14);
DefineCentroid(12, 'I', 52, 69.6, 13.3, 5.13, 18, 15);
DefineCentroid(13, 'N', 97.3, 95.9, 12.2, 9.24, 0, 18);
DefineCentroid(14, 'I', 80.1, 91.7, 7.99, 13.7, 7, 18);
DefineCentroid(15, 'K', 39.8, 40.2, 11.3, 6.96, 10, 15);
DefineCentroid(16, 'V', 119, 136, 11.1, 6.33, 13, 21);
DefineCentroid(17, 'V', 92.6, 127, 13.6, 7.08, 1, 13);
DefineCentroid(18, 'I', 6.8, 19.8, 11.8, 13.8, 11, 25);
DefineCentroid(19, 'D', 110, 27.5, 6.24, 6.65, 10, 22);
	}

uint XBinner::GetLetter(char AminoChar,
  const vector<double> &FeatureValues) const
	{
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
		if (Score > BestScore)
			{
			BestScore = Score;
			BestLetter = Letter;
			}
		}
	asserta(BestLetter != UINT_MAX);
	return BestLetter;
	}

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
#pragma omp parallel for
	for (int iPairIndex = 0; iPairIndex < (int) PosPairCount; ++iPairIndex)
		{
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
		//Log("Freqs[%c] = %8.6f\n", c, Freq);
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
		if (Score > BestScore)
			{
			BestScore = Score;
			BestLetter = Letter;
			}
		}
	asserta(BestLetter != UINT_MAX);
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

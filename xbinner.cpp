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
DefineCentroid(0, 'T', 112, 35.3, 6.37, 5.19, 16, 11); // 0.7196
DefineCentroid(1, 'H', 2.83, 26, 13.5, 11, 10, 13); // 0.7196
DefineCentroid(2, 'M', 30.6, 29.1, 11.5, 13, 21, 14); // 0.7196
DefineCentroid(3, 'N', 110, 85.1, 6.29, 9.18, 0, 19); // 0.7196
DefineCentroid(4, 'S', 58.9, 114, 11.6, 5.75, 0, 16); // 0.7196
DefineCentroid(5, 'E', 67.9, 65, 13.3, 8.6, 7, 19); // 0.7196
DefineCentroid(6, 'V', 29.8, 57.6, 12.2, 7.94, 17, 7); // 0.7196
DefineCentroid(7, 'P', 78.1, 78.5, 7.73, 5.51, 10, 11); // 0.7196
DefineCentroid(8, 'N', 66.4, 90.4, 6.12, 11.6, 5, 9); // 0.7196
DefineCentroid(9, 'P', 92.1, 63.1, 11.4, 11.8, 2, 11); // 0.7196
DefineCentroid(10, 'R', 98.8, 34.2, 9.49, 7.85, 9, 20); // 0.7196
DefineCentroid(11, 'F', 134, 139, 10.6, 7.56, 13, 10); // 0.7196
DefineCentroid(12, 'G', 116, 134, 8.94, 12.5, 20, 21); // 0.7196
DefineCentroid(13, 'P', 45.6, 25.8, 6.24, 12.3, 12, 18); // 0.7196
DefineCentroid(14, 'D', 52.6, 56.3, 9.34, 13.2, 4, 22); // 0.7196
DefineCentroid(15, 'D', 41, 59.3, 5.5, 9.04, 14, 9); // 0.7196
DefineCentroid(16, 'V', 83.4, 136, 12.6, 10.6, 15, 9); // 0.7196
DefineCentroid(17, 'V', 70, 60.5, 7.54, 9.78, 24, 15); // 0.7196
DefineCentroid(18, 'D', 114, 149, 12.5, 8.91, 7, 11); // 0.7196
DefineCentroid(19, 'Q', 112, 52.4, 0, 6.21, 12, 21); // 0.7196
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

		uint AminoLetterQ = g_CharToLetterAmino[AminoQ];
		uint AminoLetterR = g_CharToLetterAmino[AminoR];

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

	double BestScore = -999;
	uint BestLetter = UINT_MAX;
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		uint Idx = m_Idxs[Letter];
		const vector<double> &v = m_ptrX2->m_FeatureValuesVec[Idx];
		uint AminoLetter2 = g_CharToLetterAmino[g_Aminos[Letter]];
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

#include "myutils.h"
#include "alpha.h"
#include "xprof.h"
#include "xbinner.h"

vector<vector<double> > XBinner::m_ValuesVec;

double XBinner::GetValue(uint Letter, uint FeatureIndex)
	{
	asserta(Letter < 20);
	asserta(FeatureIndex < XFEATS);
	return m_ValuesVec[Letter][FeatureIndex];
	}

void XBinner::DeltaValue(uint Letter, uint FeatureIndex, double Fract)
	{
	asserta(Letter < 20);
	asserta(FeatureIndex < XFEATS);
	double Value = m_ValuesVec[Letter][FeatureIndex];
	if (fabs(Value) < 0.01)
		m_ValuesVec[Letter][FeatureIndex] += Fract;
	else
		m_ValuesVec[Letter][FeatureIndex] *= (1 + Fract);
	}

void XBinner::SetCentroid(uint Idx, const vector<double> &Values)
	{
	asserta(Idx < 20);
	if (m_ValuesVec.empty())
		{
		m_ValuesVec.resize(20);
		for (uint i = 0; i < 20; ++i)
			m_ValuesVec[i].resize(XFEATS);
		}
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		m_ValuesVec[Idx][FeatureIndex] = Values[FeatureIndex];
	}

void XBinner::DefineCentroid(uint Idx,
	double v0, double v1, double v2,
	double v3, double v4, double v5)
	{
	if (m_ValuesVec.empty())
		{
		m_ValuesVec.resize(20);
		for (uint i = 0; i < 20; ++i)
			m_ValuesVec[i].resize(XFEATS);
		}
	asserta(Idx < 20);

	m_ValuesVec[Idx][0] = v0;
	m_ValuesVec[Idx][1] = v1;
	m_ValuesVec[Idx][2] = v2;
	m_ValuesVec[Idx][3] = v3;
	m_ValuesVec[Idx][4] = v4;
	m_ValuesVec[Idx][5] = v5;
	}

void XBinner::LogCentroids()
	{
	Log("XBinner::LogCentroids\n");
	for (uint i = 0; i < 20; ++i)
		{
		for (uint j = 0; j < XFEATS; ++j)
			Log(" %10.3g", m_ValuesVec[i][j]);
		Log("\n");
		}
	}

void XBinner::InitCentroids()
	{
//DefineCentroid(0, 115, 30.2, 6.21, 5.84, 23, 21); // 0.7614
//DefineCentroid(1, 7.98, 7.17, 13.9, 13, 19, 19); // 0.7614
//DefineCentroid(2, 42.1, 52, 8.41, 12.8, 11, 12); // 0.7614
//DefineCentroid(3, 89.5, 103, 12.3, 13, 7, 15); // 0.7614
//DefineCentroid(4, 30, 43, 11, 10.8, 0, 19); // 0.7614
//DefineCentroid(5, 116, 123, 6.3, 11.9, 15, 13); // 0.7614
//DefineCentroid(6, 74.5, 81.2, 10.4, 5.81, 13, 13); // 0.7614
//DefineCentroid(7, 110, 89.9, 8.54, 9.96, 13, 9); // 0.7614
//DefineCentroid(8, 138, 111, 11, 6.47, 1, 15); // 0.7614
//DefineCentroid(9, 59.7, 59.9, 12.9, 11.5, 20, 21); // 0.7614
//DefineCentroid(10, 53.2, 59.1, 9.1, 9.05, 4, 25); // 0.7614
//DefineCentroid(11, 80.1, 43.7, 8.22, 7.58, 25, 14); // 0.7614
//DefineCentroid(12, 46.5, 87.2, 6.22, 13.9, 18, 24); // 0.7614
//DefineCentroid(13, 27.1, 76.4, 6.58, 10.3, 5, 15); // 0.7614
//DefineCentroid(14, 16.2, 17.8, 12.9, 10.3, 21, 13); // 0.7614
//DefineCentroid(15, 102, 126, 10.5, 9.3, 1, 11); // 0.7614
//DefineCentroid(16, 29.8, 57.6, 13.5, 6.69, 14, 14); // 0.7614
//DefineCentroid(17, 76, 132, 12.9, 5.68, 3, 10); // 0.7614
//DefineCentroid(18, 9.98, 31.4, 10.4, 13.7, 23, 31); // 0.7614
//DefineCentroid(19, 127, 147, 9.12, 12.7, 0, 8); // 0.7614
DefineCentroid( 0,    46.50,    87.20,     6.22,    13.90,    18.00,    24.00); // xpermute[12]
DefineCentroid( 1,   116.00,   123.00,     6.30,    11.90,    15.00,    13.00); // xpermute[5]
DefineCentroid( 2,    59.70,    59.90,    12.90,    11.50,    20.00,    21.00); // xpermute[9]
DefineCentroid( 3,     9.98,    31.40,    10.40,    13.70,    23.00,    31.00); // xpermute[18]
DefineCentroid( 4,    53.20,    59.10,     9.10,     9.05,     4.00,    25.00); // xpermute[10]
DefineCentroid( 5,    42.10,    52.00,     8.41,    12.80,    11.00,    12.00); // xpermute[2]
DefineCentroid( 6,    29.80,    57.60,    13.50,     6.69,    14.00,    14.00); // xpermute[16]
DefineCentroid( 7,    89.50,   103.00,    12.30,    13.00,     7.00,    15.00); // xpermute[3]
DefineCentroid( 8,    27.10,    76.40,     6.58,    10.30,     5.00,    15.00); // xpermute[13]
DefineCentroid( 9,   110.00,    89.90,     8.54,     9.96,    13.00,     9.00); // xpermute[7]
DefineCentroid(10,    16.20,    17.80,    12.90,    10.30,    21.00,    13.00); // xpermute[14]
DefineCentroid(11,     7.98,     7.17,    13.90,    13.00,    19.00,    19.00); // xpermute[1]
DefineCentroid(12,   102.00,   126.00,    10.50,     9.30,     1.00,    11.00); // xpermute[15]
DefineCentroid(13,   127.00,   147.00,     9.12,    12.70,     0.00,     8.00); // xpermute[19]
DefineCentroid(14,    76.00,   132.00,    12.90,     5.68,     3.00,    10.00); // xpermute[17]
DefineCentroid(15,    74.50,    81.20,    10.40,     5.81,    13.00,    13.00); // xpermute[6]
DefineCentroid(16,    30.00,    43.00,    11.00,    10.80,     0.00,    19.00); // xpermute[4]
DefineCentroid(17,    80.10,    43.70,     8.22,     7.58,    25.00,    14.00); // xpermute[11]
DefineCentroid(18,   138.00,   111.00,    11.00,     6.47,     1.00,    15.00); // xpermute[8]
DefineCentroid(19,   115.00,    30.20,     6.21,     5.84,    23.00,    21.00); // xpermute[0]
	}

uint XBinner::GetLetter(const vector<double> &FeatureValues) const
	{
	double BestScore = -999;
	uint BestLetter = UINT_MAX;
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		const vector<double> &v = m_ValuesVec[Letter];
		double Score = XProf::GetScore(FeatureValues, v);
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

		const vector<double> &ValuesQ = ValuesQVec[PairIndex];
		const vector<double> &ValuesR = ValuesRVec[PairIndex];

		uint iq = GetLetter(ValuesQ);
		uint ir = GetLetter(ValuesR);
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
		double Freq = double(CountVec[i])/(LetterPairCount);
		SumFreq += Freq;
		Freqs.push_back(Freq);
		}

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

uint XBinnerC::GetLetter(const vector<double> &FeatureValues) const
	{
	asserta(m_ptrX2 != 0);
	asserta(SIZE(m_Idxs) == 20);

	double BestScore = -999;
	uint BestLetter = UINT_MAX;
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		uint Idx = m_Idxs[Letter];
		const vector<double> &v = m_ptrX2->m_FeatureValuesVec[Idx];
		double Score = XProf::GetScore(FeatureValues, v);
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
		for (uint j = 0; j < XFEATS; ++j)
			Log(" %10.3g", m_ValuesVec[i][j]);
		Log("\n");
		}
	}

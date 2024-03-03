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
DefineCentroid(0, 109, 28.2, 6.17, 6.27, 10, 19); // 0.7553
DefineCentroid(1, 109, 86.7, 6.36, 11.3, 18, 26); // 0.7553
DefineCentroid(2, 23.1, 46.2, 13.3, 6.55, 16, 14); // 0.7553
DefineCentroid(3, 7.14, 8.73, 13.3, 11.7, 17, 19); // 0.7553
DefineCentroid(4, 71.6, 55.2, 8.16, 6.24, 8, 12); // 0.7553
DefineCentroid(5, 133, 103, 9.29, 5.84, 4, 19); // 0.7553
DefineCentroid(6, 52.8, 29.7, 8.8, 13.8, 12, 18); // 0.7553
DefineCentroid(7, 49.8, 51.3, 13.8, 12.6, 6, 24); // 0.7553
DefineCentroid(8, 101, 132, 11.8, 8.64, 0, 12); // 0.7553
DefineCentroid(9, 112, 125, 13.3, 11.3, 4, 21); // 0.7553
DefineCentroid(10, 68.8, 87.4, 9.49, 12.3, 9, 31); // 0.7553
DefineCentroid(11, 33.8, 33, 7.74, 8.71, 4, 19); // 0.7553
DefineCentroid(12, 53.8, 66, 11.6, 5.8, 21, 8); // 0.7553
DefineCentroid(13, 16.6, 71.6, 8.68, 12.1, 3, 25); // 0.7553
DefineCentroid(14, 9.11, 30.3, 11.6, 13.1, 13, 26); // 0.7553
DefineCentroid(15, 80.1, 85.6, 13.8, 9.79, 8, 24); // 0.7553
DefineCentroid(16, 66.5, 51.3, 5.86, 9.02, 10, 9); // 0.7553
DefineCentroid(17, 47.2, 85.2, 7.36, 13.7, 4, 12); // 0.7553
DefineCentroid(18, 132, 135, 10, 11.9, 0, 13); // 0.7553
DefineCentroid(19, 104, 56.2, 8.85, 8.92, 12, 19); // 0.7553
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

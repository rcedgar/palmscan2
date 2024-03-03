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
DefineCentroid(0, 114, 36, 6.44, 6.03, 17, 16); // 0.7407
DefineCentroid(1, 9.95, 24.5, 12.2, 10, 3, 22); // 0.7407
DefineCentroid(2, 42.3, 40.7, 12.3, 12.6, 1, 7); // 0.7407
DefineCentroid(3, 84.2, 68.3, 9.16, 12.5, 20, 18); // 0.7407
DefineCentroid(4, 61.7, 76.5, 12.3, 9.39, 0, 11); // 0.7407
DefineCentroid(5, 103, 74, 8.62, 9.63, 6, 15); // 0.7407
DefineCentroid(6, 107, 126, 10.5, 6.05, 7, 11); // 0.7407
DefineCentroid(7, 60.2, 112, 7.02, 12.7, 5, 19); // 0.7407
DefineCentroid(8, 78, 17.3, 8.09, 6.31, 15, 9); // 0.7407
DefineCentroid(9, 34, 22, 9.05, 10.8, 26, 18); // 0.7407
DefineCentroid(10, 36.4, 54.6, 13, 6.22, 24, 13); // 0.7407
DefineCentroid(11, 95.4, 104, 13.7, 9.28, 7, 23); // 0.7407
DefineCentroid(12, 83, 50.2, 11.2, 7.5, 5, 20); // 0.7407
DefineCentroid(13, 3.56, 0.492, 13.8, 11.4, 21, 16); // 0.7407
DefineCentroid(14, 119, 145, 13.3, 9.61, 0, 11); // 0.7407
DefineCentroid(15, 15.2, 17.8, 12.9, 13.3, 30, 27); // 0.7407
DefineCentroid(16, 121, 127, 7.22, 10.8, 7, 10); // 0.7407
DefineCentroid(17, 102, 72.2, 6.2, 12.7, 1, 12); // 0.7407
DefineCentroid(18, 97.6, 125, 11.5, 11.7, 0, 15); // 0.7407
DefineCentroid(19, 75.7, 63.6, 5.91, 7.31, 20, 13); // 0.7407
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

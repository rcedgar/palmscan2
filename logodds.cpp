#include "myutils.h"
#include "logodds.h"
#include "alpha.h"

void LogOdds::GetSymbol(uint Letter, string &s) const
	{
	s.clear();
	if (m_AlphaSize <= 20)
		{
		s += char(g_LetterToCharAmino[Letter]);
		return;
		}
	else if (m_AlphaSize <= 26*2)
		{
		if (Letter < 26)
			{
			s += 'A' + Letter;
			return;
			}
		else if (Letter < 26*2)
			{
			s += 'a' + Letter;
			return;
			}
		else
			asserta(false);
		}
	Ps(s, "%02x", Letter);
	}

void LogOdds::Init(uint AlphaSize)
	{
	m_AlphaSize = AlphaSize;
	m_CountMx.clear();
	m_Freqs.clear();
	m_FreqMx.clear();
	m_CountMx.resize(m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		m_CountMx[i].resize(m_AlphaSize);
	}

void LogOdds::AddPair(uint Letter1, uint Letter2)
	{
	assert(m_AlphaSize > 0);
	if (Letter1 >= m_AlphaSize || Letter2 >= m_AlphaSize)
		return;
	m_CountMx[Letter1][Letter2] += 1;
	m_CountMx[Letter2][Letter1] += 1;
	}

void LogOdds::GetCounts(vector<uint> &Counts) const
	{
	Counts.clear();
	Counts.resize(m_AlphaSize);
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			uint Count = m_CountMx[Letter1][Letter2];
			Counts[Letter1] += Count;
			Counts[Letter2] += Count;
			}
		}
	}

uint LogOdds::GetTotal() const
	{
	uint Total = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			uint Count = m_CountMx[Letter1][Letter2];
			Total += Count;
			}
		}
	return Total;
	}

void LogOdds::GetFreqs(vector<double> &Freqs) const
	{
	Freqs.clear();
	vector<uint> Counts(m_AlphaSize);
	GetCounts(Counts);
	uint Total = GetTotal();

	double Sum = 0;
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		double Freq = Counts[Letter]/double(2*Total);
		Freqs.push_back(Freq);
		Sum += Freq;
		}
	asserta(feq(Sum, 1.0));
	}

void LogOdds::GetFreqMx(vector<vector<double> > &Mx) const
	{
	Mx.clear();
	Mx.resize(m_AlphaSize);
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		Mx[Letter].resize(m_AlphaSize);

	uint Total = GetTotal();
	uint SumCount = 0;
	double SumFreq = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			uint Count = m_CountMx[Letter1][Letter2];
			double Freq = Count/double(Total);
			Mx[Letter1][Letter2] = Freq;
			SumCount += Count;
			SumFreq += Freq;
			}
		}
	asserta(SumCount == Total);
	asserta(feq(SumFreq, 1.0));
	}

double LogOdds::GetLogOddsMx(vector<vector<double> > &Mx) const
	{
	Mx.clear();
	Mx.resize(m_AlphaSize);
	vector<double> Freqs;
	GetFreqs(Freqs);
	vector<vector<double> > FreqMx;
	GetFreqMx(FreqMx);
	uint Total = GetTotal();
	double SumFreq = 0;
	double ExpectedScore = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		Mx[Letter1].resize(m_AlphaSize);
		double f1 = Freqs[Letter1];
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			double f2 = Freqs[Letter2];
			double ObsFreq = FreqMx[Letter1][Letter2];
			double ExpectedFreq = double(f1*f2);
			if (ObsFreq == 0 || ExpectedFreq == 0)
				continue;
			double Ratio = ObsFreq/ExpectedFreq;
			double Score = log(Ratio);
			Mx[Letter1][Letter2] = Score;
			ExpectedScore += ObsFreq*Score;
			SumFreq += ObsFreq;
			}
		}
	asserta(feq(SumFreq, 1.0));
	return ExpectedScore;
	}

void LogOdds::VecToSrc(FILE *f, const string &Name, 
  const vector<double> &v) const
	{
	asserta(SIZE(v) == m_AlphaSize);
	fprintf(f, "static double %s[%u] = {\n",
	  Name.c_str(), m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		string si;
		GetSymbol(i, si);
		fprintf(f, "\t%.4g, // %s\n", v[i], si.c_str());
		}
	fprintf(f, "};\n");
	}

void LogOdds::MxToSrc(FILE *f, const string &Name, 
  const vector<vector<double> > &Mx) const
	{
	if (f == 0)
		return;
	asserta(SIZE(Mx) == m_AlphaSize);

	fprintf(f, "static double %s[%u][%u] = {\n",
	  Name.c_str(), m_AlphaSize, m_AlphaSize);
	fprintf(f, "//        ");
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		string si;
		GetSymbol(i, si);
		fprintf(f, "         %3.3s", si.c_str());
		}
	fprintf(f, "\n");
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		string si;
		GetSymbol(i, si);
		asserta(SIZE(Mx[i]) == m_AlphaSize);
		fprintf(f, "/* %3.3s */ {", si.c_str());
		for (uint j = 0; j < m_AlphaSize; ++j)
			{
			fprintf(f, " %10.4g", Mx[i][j]);
			if (j+1 != m_AlphaSize)
				fprintf(f, ",");
			}
		fprintf(f, "}, // %s\n", si.c_str());
		}
	fprintf(f, "};\n");
	}

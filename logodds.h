#pragma once

class LogOdds
	{
public:
	uint m_AlphaSize = 0;
	vector<vector<uint> > m_CountMx;
	vector<double> m_Freqs;
	vector<double> m_FreqMx;

public:
	void Init(uint AlphaSize);
	void AddPair(uint Letter1, uint Letter2);
	void GetCounts(vector<uint> &Counts) const;
	void GetFreqs(vector<double> &Freqs) const;
	void GetFreqMx(vector<vector<double> > &Mx) const;
	double GetLogOddsMx(vector<vector<double> > &Mx) const;
	uint GetTotal() const;
	void MxToSrc(FILE *f, const string &Name, 
	  const vector<vector<double> > &Mx) const;
	void VecToSrc(FILE *f, const string &Name, 
	  const vector<double> &v) const;
	void GetSymbol(uint Letter, string &s) const;
	};

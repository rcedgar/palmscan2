#pragma once

#include "x2data.h"

class XBinner
	{
public:
	virtual uint GetLetter(char AminoChar,
	  const vector<double> &FeatureValues) const;

	void GetFreqs(
		const X2Data &X2,
		vector<double> &Freqs,
		vector<vector<double> > &FreqMx) const;

public:
	static void InitCentroids();
	static void DefineCentroid(uint Idx, char aa,
	  double v0, double v1, double v2,
	  double v3, double v4, double v5);

	static double GetLogOddsMx(const vector<double> &Freqs,
	  const vector<vector<double> > &FreqMx,
	  vector<vector<double> > &ScoreMx);
	};

class XBinnerC : public XBinner
	{
public:
	const X2Data *m_ptrX2 = 0;
	vector<uint> m_Idxs;

public:
	virtual uint GetLetter(char AminoChar,
	  const vector<double> &FeatureValues) const;
	};

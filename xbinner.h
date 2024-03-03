#pragma once

#include "x2data.h"

class XBinner
	{
public:
	static vector<vector<double> > m_ValuesVec;

public:
	virtual uint GetLetter(const vector<double> &FeatureValues) const;

// Cannot be static coz calls virtual GetLetter()
	void GetFreqs(
		const X2Data &X2,
		vector<double> &Freqs,
		vector<vector<double> > &FreqMx) const;

public:
	static const vector<vector<double> > &GetValuesVec()
		{
		return m_ValuesVec;
		}

	static double GetValue(uint Letter, uint FeatureIndex);
	static void DeltaValue(uint Letter, 
	  uint FeatureIndex, double Fract);

	static void InitCentroids();
	static void LogCentroids();

	static void DefineCentroid(uint Idx,
	  double v0, double v1, double v2,
	  double v3, double v4, double v5);

	static void SetCentroid(uint Idx, const vector<double> &Values);

	static double GetLogOddsMx(const vector<double> &Freqs,
	  const vector<vector<double> > &FreqMx,
	  vector<vector<double> > &ScoreMx);
	};

class XBinnerC : public XBinner
	{
public:
	const X2Data *m_ptrX2 = 0;
	vector<uint> m_Idxs;
	void LogMyCentroids() const;

public:
	virtual uint GetLetter(const vector<double> &FeatureValues) const;
	};

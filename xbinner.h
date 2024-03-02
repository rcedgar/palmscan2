#pragma once

class XBinner
	{
public:
	virtual uint GetLetter(uint AminoLetter,
	  const vector<double> &FeatureValues) const;

public:
	static void InitCentroids();
	static void DefineCentroid(uint Idx, char aa,
	  double v0, double v1, double v2,
	  double v3, double v4, double v5);
	};

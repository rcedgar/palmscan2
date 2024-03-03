#pragma once

#include "x2data.h"
#include <map>

class XCluster
	{
public:
	X2Data m_X2;

public:
	void ReadFeatureTsv(const string &FileName);
	void Cluster(const vector<uint> &InputIdxs,
		double MinScore,
		vector<uint> &CentroidIdxs,
		map<uint, uint> &IdxToCentroidIdx,
		map<uint, vector<uint> > &CentroidIdxToMemberIdxs);
	void GetClusterSizes(
		const map<uint, vector<uint> > &CentroidIdxToMemberIdxs,
		const vector<uint> &CentroidIdxs,
		vector<uint> &Order,
		vector<uint> &Sizes) const;
	double GetScore(uint Idx1, uint Idx2) const;
	void LogDCs(const vector<uint> AlphaIdxs, double ExpScore) const;
	void DefineBinnerLetters(const vector<uint> AlphaIdxs) const;
	};

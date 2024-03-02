#pragma once

#include <map>

class XCluster
	{
public:
// Uniques
	vector<char> m_Aminos;
	vector<vector<double> > m_FeatureValuesVec;
	vector<string> m_Doms;
	vector<uint> m_Coords;

// Aligned pairs
	vector<char> m_Aminos1;
	vector<char> m_Aminos2;
	vector<vector<double> > m_FeatureValuesVec1;
	vector<vector<double> > m_FeatureValuesVec2;

	vector<uint> m_CountVec;
	vector<vector<uint> > m_CountMx;
	vector<double> m_Freqs;
	vector<vector<double> > m_FreqMx;

public:
	void ReadFeatureTsv(const string &FileName);
	void Cluster(const vector<uint> &InputIdxs,
		double MinScore,
		vector<uint> &CentroidIdxs,
		map<uint, uint> &IdxToCentroidIdx,
		map<uint, vector<uint> > &CentroidIdxToMemberIdxs);
	void HitToTsv(double Score, uint Idx, uint CentroidIdx) const;
	void CentroidToTsv(uint CentroidIdx) const;
	void VecToTsv(FILE *f, uint Idx) const;
	void ClustersToTsv(
	  const vector<uint> &CentroidIdxs,
	  const map<uint, vector<uint> > &CentroidIdxToMemberIdxs) const;
	void GetClusterSizes(
		const map<uint, vector<uint> > &CentroidIdxToMemberIdxs,
		const vector<uint> &CentroidIdxs,
		vector<uint> &Order,
		vector<uint> &Sizes) const;
	double GetScore(uint Idx1, uint Idx2) const;
	void GetScoreMx(const vector<uint> &Idxs,
	  vector<vector<double> > &ScoreMx) const;
	void ScoreMxToDistMx(
		const vector<vector<double> > &ScoreMx,
		vector<vector<double> > &DistMx);
	void MxToTsv(FILE *f, const vector<vector<double> > &Mx) const;
	uint GetTopIdx(uint QueryIdx, const vector<uint> &RefIdxs,
	  double *ptrScore = 0, uint *ptrIndex = 0) const;
	uint GetMy3DLetter(char Amino, const vector<double> &Values,
	  const vector<uint> &RefIdxs) const;
	void TopsToTsv(FILE *f, const vector<uint> &RefIdxs) const;
	void GetFreqs(const vector<uint> &RefIdxs);
	uint GetFirstRandomHit(uint Letter, const vector<uint> &RefIdxs) const;
	void LogDCs(const vector<uint> SelectedCentroidIdxs,
	  double ExpScore) const;
	};

#pragma once

class XCluster
	{
public:
	vector<char> m_Aminos;
	vector<vector<double> > m_FeatureValuesVec;
	vector<string> m_Doms;
	vector<uint> m_Coords;

public:
	void ReadFeatureTsv(const string &FileName);
	void Cluster(
		const vector<uint> &InputIdxs,
		double MinScore,
		vector<uint> &CentroidIdxs,
		vector<uint> &IdxToCentroidIdx,
		vector<vector<uint> > &CentroidIdxToMemberIdxs);
	void HitToTsv(double Score, uint Idx, uint CentroidIdx) const;
	void CentroidToTsv(uint CentroidIdx) const;
	void VecToTsv(FILE *f, uint Idx) const;
	void ClustersToTsv(
	  const vector<uint> &CentroidIdxs,
	  const vector<vector<uint> > &CentroidIdxToMemberIdxs) const;
	void GetClusterSizes(
		const vector<vector<uint> > &CentroidIdxToMemberIdxs,
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
	void TopsToTsv(FILE *f, const vector<uint> &RefIdxs) const;
	};

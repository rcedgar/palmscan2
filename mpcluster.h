#pragma once

#include <set>
#include "abcxyz.h"

const uint MPL = AL + BL + CL;

class MotifProfile
	{
public:
	string m_Seq;
	vector<vector<float> > m_FreqVec;

public:
	MotifProfile()
		{
		m_FreqVec.resize(MPL);
		for (uint i = 0; i < MPL; ++i)
			m_FreqVec[i].resize(20, 0.0);
		}

	void ValidateFreqs() const
		{
		asserta(m_FreqVec.size() == MPL);
		for (uint i = 0; i < MPL; ++i)
			{
			const vector<float> &v = m_FreqVec[i];
			asserta(v.size() == 20);
			double Sum = 0;
			for (uint j = 0; j < 20; ++j)
				Sum += v[j];
			asserta(feq(Sum, 1.0));
			}
		}

	void Clear()
		{
		m_Seq.clear();
		asserta(m_FreqVec.size() == MPL);
		for (uint i = 0; i < MPL; ++i)
			{
			for (uint j = 0; j < 20; ++j)
				m_FreqVec[i][j] = 0;
			}
		}

	void FromXxxSeq(const string &Seq);
	void FromSeqs(const vector<string> &Seqs);

	void LogMe() const;

	char PosToMotif(uint i) const
		{
		if (i <= AL)
			return 'A';
		if (i <= AL + BL)
			return 'B';
		return 'C';
		}

	uint PosToOffset(uint i) const
		{
		if (i < AL)
			return i;
		if (i < AL + BL)
			return i - AL;
		return i - (AL + BL);
		}

	void GetLogo(string &Logo) const;
	void GetMaxLetter(uint i, uint &Letter, float &Freq) const;

public:
	static void GetLettersFromXxxSeq(const string &Seq,
	  vector<uint> &Letters);
	static void GetLettersFromSeq(const string &Seq,
	  vector<uint> &Letters);
	};

// Motif profile clustering
class MPCluster
	{
public:
	const vector<MotifProfile *> *m_Input = 0;

	set<uint> m_PendingIndexes;

// Greedy cluster
	float m_MinScore = FLT_MAX;
	vector<MotifProfile *> m_Centroids;
	vector<uint> m_CentroidIndexes;
	vector<vector<uint> > m_CentroidIndexToMemberIndexes;
	vector<uint> m_ClusterSizeOrder;

// NN cluster
	vector<MotifProfile *> m_MPs;
	vector<uint> m_Parents;
	vector<uint> m_Lefts;
	vector<uint> m_Rights;
	vector<uint> m_Sizes;

public:
	void Clear()
		{
		m_Input = 0;
		m_Centroids.clear();
		m_MinScore = FLT_MAX;
		m_PendingIndexes.clear();
		m_CentroidIndexes.clear();
		m_CentroidIndexToMemberIndexes.clear();
		m_ClusterSizeOrder.clear();
		m_MPs.clear();
		m_Parents.clear();
		m_Lefts.clear();
		m_Rights.clear();
		m_Sizes.clear();
		}

	void NNCluster(const vector<MotifProfile *> &Input,
	  float MinScore);
	void FindNN(uint &Index1, uint &Index2) const;
	void Join(uint Index1, uint Index2);
	MotifProfile &CreateProfileNN(uint Index1, uint Index2) const;

	void GreedyCluster(const vector<MotifProfile *> &Input,
	  float MinScore);
	float GetScore(const MotifProfile &MP1,
	  const MotifProfile &MP2) const;
	float GetScoreNNPair(uint i1, uint i2) const;
	void LogPair(const MotifProfile &MP1,
	  const MotifProfile &MP2) const;
	MotifProfile &GetProfile(uint i) const
		{
		const vector<MotifProfile *> &v = *m_Input;
		asserta(i < SIZE(v));
		return *v[i];
		}
	void LogClusters() const;
	void LogCluster(uint i) const;
	void ClustersToTsv(FILE *f) const;
	void ClusterToTsv(FILE *f, uint ClusterIndex) const;

public:
	static void ReadMPs(const string &FileName,
	  vector<MotifProfile *> &MPs);
	static void ReadSeqsVec(const string &TsvFileName,
	  vector<vector<string> > &SeqsVec);
	static void LogLogos(const vector<MotifProfile *> &MPs);
	};

#pragma once

#include "pdbchain.h"
#include "cddata.h"
#include "cdtemplate.h"

class CDSearcher
	{
public:
	const PDBChain *m_Query = 0;
	const CDInfo *m_Info = 0;
	const CDData *m_Dists = 0;
	const CDData *m_StdDevs = 0;
	const CDTemplate *m_Template = 0;

	vector<uint> m_MotifStarts;
	vector<double> m_MotifScores;
	vector<vector<uint> > m_MotifIndexToHits;
	vector<vector<double> > m_MotifIndexToScores;

public:
	void Clear()
		{
		m_Query = 0;
		m_Info = 0;
		m_Dists = 0;
		m_StdDevs = 0;
		m_MotifStarts.clear();
		m_MotifScores.clear();
		m_MotifIndexToHits.clear();
		m_MotifIndexToScores.clear();
		}

	void LogHits() const;
	void ClearSearch();
	void Init(const CDInfo &Info, const CDData &Dists,
	  const CDData &StdDevs);
	double GetScore(uint SeqPos1, uint SeqPos2,
	  uint Ix1, uint Ix2,
	  uint L1, uint L2) const;
	void Search(const PDBChain &Q);
	void GetRangeNext(uint NextMotifIndex, uint &Lo, uint &Hi) const;
	void Search1(uint MotifIndex, uint Lo, uint Hi,
	  vector<uint> &Hits, vector<double> &Scores);
	};

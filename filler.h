#pragma once

#include "rdrpsearcher.h"
#include "pdbchain.h"

class Filler
	{
public:
	const int m_MinX = -45;
	const int m_MaxX = 34;
	const int m_MinY = -27;
	const int m_MaxY = 41;
	const int m_MinZ = -10;
	const int m_MaxZ = 46;
	int m_RangeX = INT_MAX;
	int m_RangeY = INT_MAX;
	int m_RangeZ = INT_MAX;
	vector<vector<vector<uint> > > m_CountMx;
	vector<vector<vector<double> > > m_ScoreMx;
	vector<double> m_PositiveScores;

	RdRpModel m_RdRpModel;
	RdRpSearcher m_RS;

	uint m_TrainCount = 0;
	uint m_NoPSSMHitCount = 0;
	uint m_PermutedCount = 0;

	vector<PDBChain *> m_TrainPPXs;

public:
	Filler()
		{
		m_RangeX = m_MaxX - m_MinX + 1;
		m_RangeY = m_MaxY - m_MinY + 1;
		m_RangeZ = m_MaxZ - m_MinZ + 1;
		}

	void InitMxs();
	void ScoreMxToFile(FILE *f) const;
	void FromLines(const vector<string> &Lines);
	double CalcScoreFromCounts(int ix, int iy, int iz) const;
	double CalcZeroScoreFromCounts(double Base,
	  int ix, int iy, int iz) const;
	void CalcScoreMx();
	void UpdateCountMx(const PDBChain &PPX);
	void UpdateCountMx1(int xn, int yn, int zn);
	double GetPPScore(const PDBChain &PPX) const;
	void Train(const vector<PDBChain *> &Chains);
	void Train1(const PDBChain &Chain);
	void GetPPX(const PDBChain &Chain,
	  uint PosA, uint PosB, uint PosC,
	  PDBChain &PPX) const;
	};

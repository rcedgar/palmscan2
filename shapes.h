#pragma once

#include "pdbchain.h"

typedef vector<vector<double> > t_Mx2;
typedef vector<vector<vector<vector<double> > > > t_Mx3;

// Includes ABC
class Shapes
	{
	const double m_MinContract = 0.9;
	const double m_MaxExpand = 1.1;

public:
// [ShapeIndex]
	vector<uint> m_Lengths;
	vector<string> m_Names;

	vector<uint> m_MinNeighborDists;
	vector<uint> m_MaxNeighborDists;

// Pair-wise residue distances in 3D
//   [ShapeIndex1][ShapeIndex2][ResidueIndex1][ResidueIndex2}
	t_Mx3 m_MeanDistMx3;
	t_Mx3 m_StdDevMx3;

	uint m_Offset_D_MotifA = UINT_MAX;
	uint m_Offset_G_MotifB = UINT_MAX;
	uint m_Offset_D_MotifC = UINT_MAX;

public:
	void Clear()
		{
		m_Names.clear();
		m_Lengths.clear();
		m_MeanDistMx3.clear();
		m_StdDevMx3.clear();
		}

	uint GetShapeCount() const { return SIZE(m_Names); }
	uint GetShapeIndex(const string &Name) const;
	void Init(const vector<string> &Names,
	  const vector<uint> &Lengths);
	void Train(const vector<PDBChain *> &Chains,
	  const vector<vector<string> > &SeqsVec);
	void InitMx2(t_Mx2 &Mx) const;
	void InitMx3(t_Mx3 &Mx) const;
	void TrainGetPosVec(const PDBChain &Chain, const vector<string> &Seqs,
	  vector<uint> &PosVec) const;

	void ToFile(const string &FileName) const;
	void Mx2ToFile(FILE *f, const string &MxName, const t_Mx2 &Mx) const;
	void Mx3ToFile(FILE *f, const string &MxName, const t_Mx3 &Mx) const;

	void FromFile(const string &FileName);
	void FromLines(const vector<string> &Lines);
	uint Mx2FromLines(const vector<string> &Lines, uint LineNr,
	  const string &MxName, t_Mx2 &Mx) const;
	uint Mx3FromLines(const vector<string> &Lines, uint LineNr,
	  const string &MxName, t_Mx3 &Mx) const;

	double GetMeanDist3(uint ShapeIndex1, uint ShapeIndex2,
	  uint Offset1, uint Offset2) const;
	double GetStdDev3(uint ShapeIndex1, uint ShapeIndex2,
	  uint Offset1, uint Offset2) const;
	void GetMinMaxDist(uint ShapeIndex, const vector<uint> &NeighborDists,
	  uint &MinDist, uint &MaxDist) const;
	};

void GetTrainingMotifs(const string &FileName,
  const vector<PDBChain *> &Chains, vector<string> &ChainLabels,
  vector<string> &MotifNames, vector<uint> &MotifLengths,
  vector<vector<string> > &MotifSeqsVec);

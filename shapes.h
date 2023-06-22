#pragma once

#include "pdbchain.h"

typedef vector<vector<double> > t_Mx2;
typedef vector<vector<vector<vector<double> > > > t_Mx3;

// Includes ABC
class Shapes
	{
public:
// [ShapeIndex]
	vector<uint> m_Lengths;
	vector<string> m_Names;

// Pair-wise start distances in primary sequence
//   [ShapeIndex1][ShapeIndex2]
	t_Mx2 m_MeanDistMx2;
	t_Mx2 m_StdDevMx2;

// Pair-wise residue distances in 3D
//   [ShapeIndex1][ShapeIndex2][ResidueIndex1][ResidueIndex2}
	t_Mx3 m_MeanDistMx3;
	t_Mx3 m_StdDevMx3;

public:
	void Clear()
		{
		m_Names.clear();
		m_Lengths.clear();
		m_MeanDistMx2.clear();
		m_StdDevMx2.clear();
		m_MeanDistMx3.clear();
		m_StdDevMx3.clear();
		}

	uint GetShapeCount() const { return SIZE(m_Names); }
	void Init(const vector<string> &Names,
	  const vector<uint> &Lengths);
	void Train(const vector<PDBChain *> &Chains,
	  const vector<vector<string> > &SeqsVec);
	void InitMx2(t_Mx2 &Mx) const;
	void InitMx3(t_Mx3 &Mx) const;
	void TrainGetPosVec(const PDBChain &Chain, const vector<string> &Seqs,
	  vector<uint> &PosVec) const;
	void Mx2ToFile(FILE *f, const string &MxName, const t_Mx2 &Mx) const;
	void Mx3ToFile(FILE *f, const string &MxName, const t_Mx3 &Mx) const;
	void Mx2FromFile(FILE *f, const string &MxName, t_Mx2 &Mx) const;
	void Mx3FromFile(FILE *f, const string &MxName, t_Mx3 &Mx) const;
	void ToFile(const string &FileName) const;
	void FromFile(const string &FileName);
	double GetMeanDist2(uint ShapeIndex1, uint ShapeIndex2) const;
	double GetStdDev2(uint ShapeIndex1, uint ShapeIndex2) const;
	double GetMeanDist3(uint ShapeIndex1, uint ShapeIndex2,
	  uint Offset1, uint Offset2) const;
	double GetStdDev3(uint ShapeIndex1, uint ShapeIndex2,
	  uint Offset1, uint Offset2) const;
	};

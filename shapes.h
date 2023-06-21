#pragma once

#include "pdbchain.h"

// Includes ABC
class Shapes
	{
public:
// [ShapeIndex]
	vector<uint> m_Lengths;
	vector<string> m_Names;

// Pair-wise start distances in primary sequence
//   [ShapeIndex1][ShapeIndex2]
	vector<vector<double> > m_Mean2DistMx;
	vector<vector<double> > m_Mean2StdDevMx;

// Pair-wise residue distances in 3D
//   [ShapeIndex1][ShapeIndex2][ResidueIndex1][ResidueIndex2}
	vector<vector<vector<vector<double> > > > m_Mean3DistMx;
	vector<vector<vector<vector<double> > > > m_Mean3StdDevMx;

public:
	void Clear()
		{
		m_Names.clear();
		m_Lengths.clear();
		m_Mean2DistMx.clear();
		m_Mean2StdDevMx.clear();
		m_Mean3DistMx.clear();
		m_Mean3StdDevMx.clear();
		}

	uint GetShapeCount() const { return SIZE(m_Names); }
	void Init(const vector<string> &Names,
	  const vector<uint> &Lengths);
	void Train(const vector<PDBChain *> &Chains,
	  const vector<vector<string> > &SeqsVec);
	void Init2Mx(vector<vector<double> > &Mx) const;
	void Init3Mx(vector<vector<vector<vector<double> > > > &Mx) const;
	void TrainGetPosVec(const PDBChain &Chain, const vector<string> &Seqs,
	  vector<uint> &PosVec) const;
	void TrainAdd1(const PDBChain &Chain, vector<uint> &PosVec,
	  vector<vector<double> > &CountMx,  
	  vector<vector<double> > &Sum2DistMx,  
	  vector<vector<vector<vector<double> > > > &Sum3DistMx) const;
	};

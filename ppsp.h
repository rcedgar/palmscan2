#pragma once

#include "pdbchain.h"
#include "abcxyz.h"

const uint PPSPL = 34;
const uint AIX = 0;
const uint BIX = AL;
const uint CIX = AL + BL;

class PPSP
	{
public:
	vector<vector<double> > m_Means;
	vector<vector<double> > m_StdDevs;
	
// For training
	vector<vector<vector<double> > > m_DistMxVec;

public:
	PPSP()
		{
		Clear();
		};

	void Clear()
		{
		m_Means.clear();
		m_StdDevs.clear();
		m_Means.resize(PPSPL);
		m_StdDevs.resize(PPSPL);
		for (uint i = 0; i < PPSPL; ++i)
			{
			m_Means[i].resize(PPSPL, DBL_MAX);
			m_StdDevs[i].resize(PPSPL, DBL_MAX);
			}
		}

	void ToFile(const string &FileName) const;
	void FromFile(const string &FileName);
	double GetScore3(const PDBChain &Chain,
	  uint PosA, uint PosB, uint PosC) const;
	double GetScoreA(const PDBChain &Chain, uint PosA) const;
	double GetScoreB(const PDBChain &Chain, uint PosB) const;
	double GetScoreC(const PDBChain &Chain, uint PosC) const;
	double GetScoreAB(const PDBChain &Chain, uint PosA, uint PosB) const;
	double GetScoreBC(const PDBChain &Chain, uint PosB, uint PosC) const;
	double GetScoreAC(const PDBChain &Chain, uint PosA, uint PosC) const;
	double GetScore(const PDBChain &Chain,
	  uint SeqPos, uint Ix, uint L, bool Trace = false) const;
	double GetScore2(const PDBChain &Chain,
	  uint SeqPos1, uint SeqPos2,
	  uint Ix1, uint Ix2,
	  uint L1, uint L2, bool Trace = false) const;
	char GetMotifChar(uint Ix) const;

// For training
public:
	void InitTrain()
		{
		m_DistMxVec.clear();
		}
	void FinalizeTrain();
	void TrainChain(const PDBChain &Chain,
	  uint APos, uint BPos, uint CPos);

	void GetDistMx(const PDBChain &Chain, 
	  uint APos, uint BPos, uint CPos,
	  vector<vector<double> > &DistMx);

	void AppendTrainDistMx(const vector<vector<double> > &DistMx)
		{
		asserta(SIZE(DistMx) == PPSPL);
		asserta(SIZE(DistMx[0]) == PPSPL);
		m_DistMxVec.push_back(DistMx);
		}

	void GetMeanStdDev(uint i, uint j,
	  double &Mean, double &StdDev) const;

public:
	static uint GetSeqPos(uint i, uint APos, uint BPos, uint CPos);
	};

double GetNormal(double Mu, double Sigma, double x);

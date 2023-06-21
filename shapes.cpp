#include "myutils.h"
#include "shapes.h"

/***
ABC found by CMP.
Two passes:
	1. For each shape, find bext fit vs. ABC
	2. For each shape, hold all other shapes fixed and re-optimize.
***/

void Shapes::Init(const vector<string> &Names,
  const vector<uint> &Lengths)
	{
	Clear();
	const uint N = SIZE(Names);
	asserta(SIZE(Lengths) == N);
	m_Names = Names;
	m_Lengths = Lengths;
	}

void Shapes::Init2Mx(vector<vector<double> > &Mx) const
	{
	uint ShapeCount = GetShapeCount();
	Mx.resize(ShapeCount);
	for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
		Mx[ShapeIndex1].resize(ShapeCount, DBL_MAX);
	}

void Shapes::Init3Mx(vector<vector<vector<vector<double> > > > &Mx) const
	{
	uint ShapeCount = GetShapeCount();
	Mx.resize(ShapeCount);
	for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
		{
		uint L1 = m_Lengths[ShapeIndex1];
		Mx[ShapeIndex1].resize(ShapeCount);
		for (uint ShapeIndex2 = 0; ShapeIndex2 < ShapeCount; ++ShapeIndex2)
			{
			uint L2 = m_Lengths[ShapeIndex2];
			Mx[ShapeIndex1][ShapeIndex2].resize(L1);
			for (uint ResidueIndex1 = 0; ResidueIndex1 < L1; ++ResidueIndex1)
				Mx[ShapeIndex1][ShapeIndex2][ResidueIndex1].resize(L2, DBL_MAX);
			}
		}
	}

void Shapes::TrainGetPosVec(const PDBChain &Chain, const vector<string> &Seqs,
  vector<uint> &PosVec) const
	{
	uint ShapeCount = GetShapeCount();
	asserta(SIZE(Seqs) == ShapeCount);
	PosVec.clear();
	for (uint ShapeIndex = 0; ShapeIndex < ShapeCount; ++ShapeIndex)
		{
		uint Pos = UINT_MAX;
		const string &Seq = Seqs[ShapeIndex];
		if (Seq != "" && Seq != ".")
			{
			size_t k = Chain.m_Seq.find(Seq);
			asserta(Chain.m_Seq.rfind(Seq) == k);
			if (k != string::npos)
				Pos = uint(k);
			}
		PosVec.push_back(Pos);
		}
	}

void Shapes::TrainAdd1(const PDBChain &Chain, vector<uint> &PosVec,
  vector<vector<double> > &CountMx,  
  vector<vector<double> > &Sum2DistMx,  
  vector<vector<vector<vector<double> > > > &Sum3DistMx) const
	{
	uint ShapeCount = GetShapeCount();
	asserta(SIZE(PosVec) == ShapeCount);
	for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
		{
		uint Pos1 = PosVec[ShapeIndex1];
		if (Pos1 == UINT_MAX)
			continue;
		uint L1 = m_Lengths[ShapeIndex1];
		for (uint ShapeIndex2 = ShapeIndex1; ShapeIndex2 < ShapeCount; ++ShapeIndex2)
			{
			uint Pos2 = PosVec[ShapeIndex2];
			if (Pos2 == UINT_MAX)
				continue;
			uint L2 = m_Lengths[ShapeIndex2];

			CountMx[ShapeIndex1][ShapeIndex2] += 1;

			double Dist2 = double(Pos2) - double(Pos1);
			Sum2DistMx[ShapeIndex1][ShapeIndex2] += Dist2;

			for (uint Offset1 = 0; Offset1 < L1; ++Offset1)
				{
				uint SeqPos1 = Pos1 + Offset1;
				for (uint Offset2 = 0; Offset2 < L2; ++Offset2)
					{
					uint SeqPos2 = Pos2 + Offset2;
					double d = Chain.GetDist(SeqPos1, SeqPos2);
					Sum2DistMx[Offset1][Offset2] = d;
					}
				}
			}
		}
	}

void Shapes::Train(const vector<PDBChain *> &Chains,
  const vector<vector<string> > &SeqsVec)
	{
	const uint ChainCount = SIZE(Chains);
	asserta(SIZE(SeqsVec) == ChainCount);

	Init2Mx(m_Mean2DistMx);
	Init2Mx(m_Mean2StdDevMx);

	Init3Mx(m_Mean3DistMx);
	Init3Mx(m_Mean3StdDevMx);

	vector<vector<double> > CountMx;
	vector<vector<double> > Sum2DistMx;
	Init2Mx(CountMx);
	Init2Mx(Sum2DistMx);

	vector<vector<vector<vector<double> > > > Sum3DistMx;
	Init3Mx(Sum3DistMx);

	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const vector<string> &Seqs = SeqsVec[ChainIndex];
		vector<uint> PosVec;
		TrainGetPosVec(Chain, Seqs, PosVec);
		TrainAdd1(Chain, PosVec, CountMx, Sum2DistMx, Sum3DistMx);
		}
	}

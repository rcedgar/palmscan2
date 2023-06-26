#include "myutils.h"
#include "shapes.h"
#include "quarts.h"

uint Shapes::GetShapeIndex(const string &Name) const
	{
	for (uint i = 0; i < SIZE(m_Names); ++i)
		if (Name == m_Names[i])
			return i;
	Die("GetShapeIndex(%s)", Name.c_str());
	return UINT_MAX;
	}

void Shapes::Init(const vector<string> &Names,
  const vector<uint> &Lengths)
	{
	Clear();
	const uint N = SIZE(Names);
	asserta(SIZE(Lengths) == N);
	m_Names = Names;
	m_Lengths = Lengths;
	}

void Shapes::InitMx2(t_Mx2 &Mx) const
	{
	uint ShapeCount = GetShapeCount();
	Mx.resize(ShapeCount);
	for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
		Mx[ShapeIndex1].resize(ShapeCount, 0);
	}

void Shapes::InitMx3(t_Mx3 &Mx) const
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
				Mx[ShapeIndex1][ShapeIndex2][ResidueIndex1].resize(L2, 0);
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

void Shapes::GetMinMaxDist(uint ShapeIndex, const vector<uint> &NeighborDists,
  uint &MinDist, uint &MaxDist) const
	{
	Quarts Q;
	GetQuarts(NeighborDists, Q);
	MinDist = uint(Q.Min*m_MinContract);
	MaxDist = uint(Q.Max*m_MaxExpand);
	}

void Shapes::Train(const vector<PDBChain *> &Chains,
  const vector<vector<string> > &SeqsVec, bool ExtendABC)
	{
	//if (ExtendABC)
	//	{
	//	m_Offset_D_MotifA = 5;
	//	m_Offset_G_MotifB = 1;
	//	m_Offset_D_MotifC = 5;
	//	}
	//else
	//	{
	//	m_Offset_D_MotifA = 3;
	//	m_Offset_G_MotifB = 1;
	//	m_Offset_D_MotifC = 3;
	//	}
	m_Offset_D_MotifA = 5;
	m_Offset_G_MotifB = 1;
	m_Offset_D_MotifC = 7;

	const uint ChainCount = SIZE(Chains);
	asserta(SIZE(SeqsVec) == ChainCount);
	const uint ShapeCount = GetShapeCount();

	vector<vector<uint> > NeighborDistVec(ShapeCount-1);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Training");
		const PDBChain &Chain = *Chains[ChainIndex];
		const vector<string> &Seqs = SeqsVec[ChainIndex];
		vector<uint> PosVec;
		TrainGetPosVec(Chain, Seqs, PosVec);
		asserta(SIZE(PosVec) == ShapeCount);
		for (uint ShapeIndex = 0; ShapeIndex + 1 < ShapeCount; ++ShapeIndex)
			{
			uint Pos = PosVec[ShapeIndex];
			if (Pos == UINT_MAX)
				continue;
			uint NextPos = PosVec[ShapeIndex + 1];
			if (NextPos == UINT_MAX)
				continue;
		// permuted?
			if (NextPos <= Pos)
				continue;
			uint Dist = NextPos - Pos;
			NeighborDistVec[ShapeIndex].push_back(Dist);
			}
		}

	m_MinNeighborDists.clear();
	m_MaxNeighborDists.clear();
	for (uint ShapeIndex = 0; ShapeIndex + 1 < ShapeCount; ++ShapeIndex)
		{
		const vector<uint> &NeighborDists = NeighborDistVec[ShapeIndex];

		uint MinDist, MaxDist;
		GetMinMaxDist(ShapeIndex, NeighborDists, MinDist, MaxDist);

		m_MinNeighborDists.push_back(MinDist);
		m_MaxNeighborDists.push_back(MaxDist);
		}
	m_MinNeighborDists.push_back(0);
	m_MaxNeighborDists.push_back(0);

	InitMx3(m_MeanDistMx3);
	InitMx3(m_StdDevMx3);

	t_Mx2 CountMx;
	InitMx2(CountMx);

	t_Mx3 SumDistMx3;
	InitMx3(SumDistMx3);

////////////////////////////////////////////////////////////////////
// Sums (3D)
////////////////////////////////////////////////////////////////////
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const vector<string> &Seqs = SeqsVec[ChainIndex];
		vector<uint> PosVec;
		TrainGetPosVec(Chain, Seqs, PosVec);
		asserta(SIZE(PosVec) == ShapeCount);
		for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
			{
			uint Pos1 = PosVec[ShapeIndex1];
			if (Pos1 == UINT_MAX)
				continue;
			uint L1 = m_Lengths[ShapeIndex1];
			for (uint ShapeIndex2 = 0; ShapeIndex2 <= ShapeIndex1; ++ShapeIndex2)
				{
				uint Pos2 = PosVec[ShapeIndex2];
				if (Pos2 == UINT_MAX)
					continue;
				uint L2 = m_Lengths[ShapeIndex2];

				CountMx[ShapeIndex1][ShapeIndex2] += 1;

				for (uint Offset1 = 0; Offset1 < L1; ++Offset1)
					{
					uint SeqPos1 = Pos1 + Offset1;
					for (uint Offset2 = 0; Offset2 < L2; ++Offset2)
						{
						uint SeqPos2 = Pos2 + Offset2;
						double d = Chain.GetDist(SeqPos1, SeqPos2);
						SumDistMx3[ShapeIndex1][ShapeIndex2][Offset1][Offset2] += d;
						}
					}
				}
			}
		}

////////////////////////////////////////////////////////////////////
// Means (3D)
////////////////////////////////////////////////////////////////////
	for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
		{
		uint L1 = m_Lengths[ShapeIndex1];
		for (uint ShapeIndex2 = 0; ShapeIndex2 <= ShapeIndex1; ++ShapeIndex2)
			{
			uint L2 = m_Lengths[ShapeIndex2];

			double Count = CountMx[ShapeIndex1][ShapeIndex2];
			if (Count == 0)
				continue;

			for (uint Offset1 = 0; Offset1 < L1; ++Offset1)
				{
				for (uint Offset2 = 0; Offset2 < L2; ++Offset2)
					{
					double SumDist =
					  SumDistMx3[ShapeIndex1][ShapeIndex2][Offset1][Offset2];
					m_MeanDistMx3[ShapeIndex1][ShapeIndex2][Offset1][Offset2] =
					  SumDist/Count;
					}
				}
			}
		}

////////////////////////////////////////////////////////////////////
// Sum (x - mean)^2 (3D)
////////////////////////////////////////////////////////////////////
	t_Mx3 SumDiffSquared3;
	InitMx3(SumDiffSquared3);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const vector<string> &Seqs = SeqsVec[ChainIndex];
		vector<uint> PosVec;
		TrainGetPosVec(Chain, Seqs, PosVec);
		asserta(SIZE(PosVec) == ShapeCount);
		for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
			{
			uint Pos1 = PosVec[ShapeIndex1];
			if (Pos1 == UINT_MAX)
				continue;
			uint L1 = m_Lengths[ShapeIndex1];
			for (uint ShapeIndex2 = 0; ShapeIndex2 <= ShapeIndex1; ++ShapeIndex2)
				{
				uint Pos2 = PosVec[ShapeIndex2];
				if (Pos2 == UINT_MAX)
					continue;
				uint L2 = m_Lengths[ShapeIndex2];

				for (uint Offset1 = 0; Offset1 < L1; ++Offset1)
					{
					uint SeqPos1 = Pos1 + Offset1;
					for (uint Offset2 = 0; Offset2 < L2; ++Offset2)
						{
						uint SeqPos2 = Pos2 + Offset2;
						double Dist = Chain.GetDist(SeqPos1, SeqPos2);
						double MeanDist =
						  m_MeanDistMx3[ShapeIndex1][ShapeIndex2][Offset1][Offset2];
						double Diff = Dist - MeanDist;
						SumDiffSquared3[ShapeIndex1][ShapeIndex2][Offset1][Offset2] +=
						  Diff*Diff;
						}
					}
				}
			}
		}

////////////////////////////////////////////////////////////////////
// Standard deviations (3D)
////////////////////////////////////////////////////////////////////
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const vector<string> &Seqs = SeqsVec[ChainIndex];
		vector<uint> PosVec;
		TrainGetPosVec(Chain, Seqs, PosVec);
		asserta(SIZE(PosVec) == ShapeCount);
		for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
			{
			uint L1 = m_Lengths[ShapeIndex1];
			for (uint ShapeIndex2 = 0; ShapeIndex2 <= ShapeIndex1; ++ShapeIndex2)
				{
				uint L2 = m_Lengths[ShapeIndex2];
				double Count = CountMx[ShapeIndex1][ShapeIndex2];
				if (Count == 0)
					continue;

				for (uint Offset1 = 0; Offset1 < L1; ++Offset1)
					{
					for (uint Offset2 = 0; Offset2 < L2; ++Offset2)
						{
						double Sum = SumDiffSquared3[ShapeIndex1][ShapeIndex2][Offset1][Offset2];
						double StdDev = sqrt(Sum/Count);
						m_StdDevMx3[ShapeIndex1][ShapeIndex2][Offset1][Offset2] = StdDev;
						}
					}
				}
			}
		}
	}

void Shapes::Mx2ToFile(FILE *f, const string &MxName,
  const t_Mx2 &Mx) const
	{
	if (f == 0)
		return;

	fprintf(f, "%s\n", MxName.c_str());
	const uint ShapeCount = GetShapeCount();
	for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
		{
		fprintf(f, "%u", ShapeIndex1);
		for (uint ShapeIndex2 = 0; ShapeIndex2 < ShapeIndex1; ++ShapeIndex2)
			fprintf(f, "\t%.6g", Mx[ShapeIndex1][ShapeIndex2]);
		fprintf(f, "\n");
		}
	}

void Shapes::Mx2FromFile(FILE *f, const string &MxName, t_Mx2 &Mx) const
	{
	InitMx2(Mx);
	string Line;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	asserta(Line == MxName);
	const uint ShapeCount = GetShapeCount();
	for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
		{
		vector<string> Fields;
		Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == ShapeIndex1 + 1);
		asserta(StrToUint(Fields[0]) == ShapeIndex1);
		for (uint ShapeIndex2 = 0; ShapeIndex2 < ShapeIndex1; ++ShapeIndex2)
			{
			double x = StrToFloat(Fields[ShapeIndex2+1]);
			Mx[ShapeIndex1][ShapeIndex2] = x;
			}
		}
	}

void Shapes::Mx3ToFile(FILE *f, const string &MxName,
  const t_Mx3 &Mx) const
	{
	if (f == 0)
		return;

	fprintf(f, "%s\n", MxName.c_str());
	const uint ShapeCount = GetShapeCount();
	for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
		{
		uint L1 = m_Lengths[ShapeIndex1];
		for (uint ShapeIndex2 = 0; ShapeIndex2 <= ShapeIndex1; ++ShapeIndex2)
			{
			uint L2 = m_Lengths[ShapeIndex2];
			for (uint Offset1 = 0; Offset1 < L1; ++Offset1)
				{
				fprintf(f, "%u\t%u\t%u", ShapeIndex1, ShapeIndex2, Offset1);
				for (uint Offset2 = 0; Offset2 < L2; ++Offset2)
					fprintf(f, "\t%.6g",
					  Mx[ShapeIndex1][ShapeIndex2][Offset1][Offset2]);
				fprintf(f, "\n");
				}
			}
		}
	}

void Shapes::Mx3FromFile(FILE *f, const string &MxName, t_Mx3 &Mx) const
	{
	InitMx3(Mx);

	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	asserta(Line == MxName);
	const uint ShapeCount = GetShapeCount();
	for (uint ShapeIndex1 = 0; ShapeIndex1 < ShapeCount; ++ShapeIndex1)
		{
		uint L1 = m_Lengths[ShapeIndex1];
		for (uint ShapeIndex2 = 0; ShapeIndex2 <= ShapeIndex1; ++ShapeIndex2)
			{
			uint L2 = m_Lengths[ShapeIndex2];
			for (uint Offset1 = 0; Offset1 < L1; ++Offset1)
				{
				Ok = ReadLineStdioFile(f, Line);
				asserta(Ok);
				Split(Line, Fields, '\t');
				asserta(SIZE(Fields) == L2 + 3);
				asserta(StrToUint(Fields[0]) == ShapeIndex1);
				asserta(StrToUint(Fields[1]) == ShapeIndex2);
				asserta(StrToUint(Fields[2]) == Offset1);
				for (uint Offset2 = 0; Offset2 < L2; ++Offset2)
					{
					double x = StrToFloat(Fields[3+Offset2]);
					Mx[ShapeIndex1][ShapeIndex2][Offset1][Offset2] = x;
					}
				}
			}
		}
	}

void Shapes::ToFile(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	uint ShapeCount = GetShapeCount();
	fprintf(f, "shapes\t%u\n", ShapeCount);
	fprintf(f, "AD_BG_CD\t%u\t%u\t%u\n",
	  m_Offset_D_MotifA,
	  m_Offset_G_MotifB,
	  m_Offset_D_MotifC);
	for (uint ShapeIndex = 0; ShapeIndex < ShapeCount; ++ShapeIndex)
		fprintf(f, "%u\t%s\t%u\t%u\t%u\n",
		  ShapeIndex, m_Names[ShapeIndex].c_str(),
		  m_Lengths[ShapeIndex],
		  m_MinNeighborDists[ShapeIndex],
		  m_MaxNeighborDists[ShapeIndex]);

	Mx3ToFile(f, "MeanDistMx3", m_MeanDistMx3);
	Mx3ToFile(f, "StdDevMx3", m_StdDevMx3);
	Pf(f, "//\n");
	CloseStdioFile(f);
	}

void Shapes::FromFile(const string &FileName)
	{
	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;

	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 2);
	asserta(Fields[0] == "shapes");
	uint ShapeCount = StrToUint(Fields[1]);

	Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 4);
	asserta(Fields[0] == "AD_BG_CD");
	m_Offset_D_MotifA = StrToUint(Fields[1]);
	m_Offset_G_MotifB = StrToUint(Fields[2]);
	m_Offset_D_MotifC = StrToUint(Fields[3]);

	for (uint ShapeIndex = 0; ShapeIndex < ShapeCount; ++ShapeIndex)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 5);
		asserta(StrToUint(Fields[0]) == ShapeIndex);
		const string &Name = Fields[1];
		uint L = StrToUint(Fields[2]);
		uint MinND = StrToUint(Fields[3]);
		uint MaxND = StrToUint(Fields[4]);
		m_Names.push_back(Name);
		m_Lengths.push_back(L);
		m_MinNeighborDists.push_back(MinND);
		m_MaxNeighborDists.push_back(MaxND);
		}

	Mx3FromFile(f, "MeanDistMx3", m_MeanDistMx3);
	Mx3FromFile(f, "StdDevMx3", m_StdDevMx3);

	Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	asserta(Line == "//");

	CloseStdioFile(f);
	}

double Shapes::GetMeanDist3(uint ShapeIndex1, uint ShapeIndex2,
  uint Offset1, uint Offset2) const
	{
	if (ShapeIndex1 <= ShapeIndex2)
		{
		assert(ShapeIndex2 < SIZE(m_MeanDistMx3));
		assert(ShapeIndex1 < SIZE(m_MeanDistMx3[ShapeIndex2]));
		assert(Offset2 < SIZE(m_MeanDistMx3[ShapeIndex2][ShapeIndex1]));
		assert(Offset1 < SIZE(m_MeanDistMx3[ShapeIndex2][ShapeIndex1][Offset2]));
		double d = m_MeanDistMx3[ShapeIndex2][ShapeIndex1][Offset2][Offset1];
		return -d;
		}
	else
		{
		assert(ShapeIndex1 > ShapeIndex2);
		assert(ShapeIndex1 < SIZE(m_MeanDistMx3));
		assert(ShapeIndex2 < SIZE(m_MeanDistMx3[ShapeIndex1]));
		assert(Offset1 < SIZE(m_MeanDistMx3[ShapeIndex1][ShapeIndex2]));
		assert(Offset2 < SIZE(m_MeanDistMx3[ShapeIndex1][ShapeIndex2][Offset1]));
		double d = m_MeanDistMx3[ShapeIndex1][ShapeIndex2][Offset1][Offset2];
		return d;
		}
	}

double Shapes::GetStdDev3(uint ShapeIndex1, uint ShapeIndex2,
  uint Offset1, uint Offset2) const
	{
	if (ShapeIndex1 <= ShapeIndex2)
		{
		assert(ShapeIndex2 < SIZE(m_MeanDistMx3));
		assert(ShapeIndex1 < SIZE(m_MeanDistMx3[ShapeIndex2]));
		assert(Offset2 < SIZE(m_MeanDistMx3[ShapeIndex2][ShapeIndex1]));
		assert(Offset1 < SIZE(m_MeanDistMx3[ShapeIndex2][ShapeIndex1][Offset2]));
		double d = m_StdDevMx3[ShapeIndex2][ShapeIndex1][Offset2][Offset1];
		return -d;
		}
	else
		{
		assert(ShapeIndex1 > ShapeIndex2);
		assert(ShapeIndex1 > ShapeIndex2);
		assert(ShapeIndex1 < SIZE(m_MeanDistMx3));
		assert(ShapeIndex2 < SIZE(m_MeanDistMx3[ShapeIndex1]));
		assert(Offset1 < SIZE(m_MeanDistMx3[ShapeIndex1][ShapeIndex2]));
		assert(Offset2 < SIZE(m_MeanDistMx3[ShapeIndex1][ShapeIndex2][Offset1]));
		double d = m_StdDevMx3[ShapeIndex1][ShapeIndex2][Offset1][Offset2];
		return d;
		}
	}

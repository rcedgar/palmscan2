#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "motifsettings.h"
#include "abcxyz.h"
#include "sort.h"

static uint g_MaxCountA = 0;
static uint g_MaxCountB = 0;
static uint g_MaxCountC = 0;

static uint GetAtomCountNotHydrogen(const PDBChain &Q, uint Pos)
	{
	vector<double> Xs;
	vector<double> Ys;
	vector<double> Zs;
	vector<string> Els;
	vector<string> Names;
	vector<string> Lines;
	Q.GetResidueAtomsInfo(Pos, Xs, Ys, Zs, Els, Names, Lines);
	const uint N = SIZE(Xs);
	uint Count = 0;
	for (uint i = 0; i < N; ++i)
		{
		const string &El = Els[i];
		if (El != "H")
			++Count;
		}
	return Count;
	}


static uint ATOMSA = 70;
static uint ATOMSB = 70;
static uint ATOMSC = 30;

static bool AppendAtomLines(const PDBChain &Q, uint Pos,
  uint AtomCount, vector<string> &Lines)
	{
	const uint QL = Q.GetSeqLength();
	uint Count = 0;
	uint k = 0;
	for (;;)
		{
		vector<string> PosLines;
		if (Pos + k >= QL)
			return false;
		Q.GetATOMLines(Pos + k, PosLines);
		++k;
		const uint n = SIZE(PosLines);
		for (uint i = 0; i < n; ++i)
			{
			const string &Line = PosLines[i];

			string El;
			PDBChain::GetElementNameFromATOMLine(Line, El);
			if (El != "H")
				{
				Lines.push_back(Line);
				++Count;
				if (Count == AtomCount)
					return true;
				}
			}
		}
	return false;
	}

void AItoms(const PDBChain &Q, uint PosA, uint PosB, uint PosC)
	{
	if (PosA == UINT_MAX) return;
	if (PosB == UINT_MAX) return;
	if (PosC == UINT_MAX) return;

	uint LoA = PosA + g_OffAd;
	uint HiA = LoA + 6;

	uint LoB = PosB;
	uint HiB = LoB + 6;

	asserta(g_OffCd > 0);
	uint LoC = PosC + g_OffCd - 1;
	uint HiC = LoC + 2;

	string PepSeqA;
	string PepSeqB;
	string PepSeqC;

	vector<string> Lines;
	bool OkA = AppendAtomLines(Q, LoA, ATOMSA, Lines);
	bool OkB = AppendAtomLines(Q, LoB, ATOMSB, Lines);
	bool OkC = AppendAtomLines(Q, LoC, ATOMSC, Lines);
	if (!OkA || !OkB || !OkC)
		return;

	const uint N = SIZE(Lines);
	asserta(N == ATOMSA + ATOMSB + ATOMSC);

	vector<double> CoordXs;
	vector<double> CoordYs;
	vector<double> CoordZs;
	for (uint i = 0; i < N; ++i)
		{
		const string &Line = Lines[i];

		double CoordX, CoordY, CoordZ;
		PDBChain::GetXYZFromATOMLine(Line, CoordX, CoordY, CoordZ);

		CoordXs.push_back(CoordX);
		CoordYs.push_back(CoordY);
		CoordZs.push_back(CoordZ);
		}

	if (g_fatomdists == 0)
		return;

#pragma omp critical
	{
	fprintf(g_fatomdists, ">%s\n", Q.m_Label.c_str());
	for (uint i = 0; i < N; ++i)
		{
		double Xi = CoordXs[i];
		double Yi = CoordYs[i];
		double Zi = CoordZs[i];
		for (uint j = 0; j < i; ++j)
			{
			double Xj = CoordXs[j];
			double Yj = CoordYs[j];
			double Zj = CoordZs[j];

			double d = GetDist3D(Xi, Yi, Zi, Xj, Yj, Zj);
			fprintf(g_fatomdists, "%u\t%u\t%.1f\n", i, j, d);
			}
		}
	}
	}

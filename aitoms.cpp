#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "motifsettings.h"
#include "abcxyz.h"
#include "sort.h"

static void GetAtomInfos(const PDBChain &Q, uint Lo, uint Hi, char MotifChar,
  vector<vector<double> > &XVec,
  vector<vector<double> > &YVec,
  vector<vector<double> > &ZVec,
  vector<vector<string> > &ElVec,
  vector<vector<string> > &NameVec,
  vector<vector<string> > &LinesVec)
	{
	XVec.clear();
	YVec.clear();
	ZVec.clear();
	ElVec.clear();
	NameVec.clear();
	LinesVec.clear();

	asserta(Lo <= Hi);
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		vector<double> Xs;
		vector<double> Ys;
		vector<double> Zs;
		vector<string> Els;
		vector<string> Names;
		vector<string> Lines;
		Q.GetResidueAtomsInfo(Pos, Xs, Ys, Zs, Els, Names, Lines);

		XVec.push_back(Xs);
		YVec.push_back(Ys);
		ZVec.push_back(Zs);
		ElVec.push_back(Els);
		NameVec.push_back(Names);
		LinesVec.push_back(Lines);
		}

#if	1
	{
	const uint n = SIZE(XVec);
	for (uint i = 0; i < n; ++i)
		{
		const vector<double> &Xs = XVec[i];
		const vector<double> &Ys = YVec[i];
		const vector<double> &Zs = ZVec[i];
		const vector<string> &Els = ElVec[i];
		const vector<string> &Names = NameVec[i];
		const uint m = SIZE(Xs);
		for (uint j = 0; j < m; ++j)
			{
			Log("[%c %4u-%4u i=%2u j=%2u] ", MotifChar, Lo, Hi, i, j);
			Log("  X=%7.1f", Xs[j]);
			Log("  Y=%7.1f", Ys[j]);
			Log("  Z=%7.1f", Zs[j]);
			Log(" %s %s\n", Els[j].c_str(), Names[j].c_str());
			}
		}
	}
#endif
	}

static double GetMinDist(double X, double Y, double Z,
  const vector<vector<double> > &XVec,
  const vector<vector<double> > &YVec,
  const vector<vector<double> > &ZVec)
	{
	const uint N = SIZE(XVec);
	asserta(SIZE(YVec) == N);
	asserta(SIZE(ZVec) == N);
	double MinDist = DBL_MAX;
	for (uint i = 0; i < N; ++i)
		{
		const vector<double> &Xs = XVec[i];
		const vector<double> &Ys = YVec[i];
		const vector<double> &Zs = ZVec[i];
		const uint M = SIZE(Xs);
		asserta(SIZE(Ys) == M);
		asserta(SIZE(Zs) == M);
		for (uint j = 0; j < M; ++j)
			{
			double d = GetDist3D(X, Y, Z, Xs[j], Ys[j], Zs[j]);
			MinDist = min(MinDist, d);
			}
		}
	return MinDist;
	}

static void SelectClosestAtoms(
  const vector<vector<vector<double> > > &XVec,
  const vector<vector<vector<double> > > &YVec,
  const vector<vector<vector<double> > > &ZVec,
  uint MotifIndex, uint ResidueIndex, uint K,
  vector<uint> &Ixs)
	{
	Ixs.clear();

	asserta(MotifIndex < 3);
	uint OtherMotifIndex1 = (MotifIndex + 1)%3;
	uint OtherMotifIndex2 = (OtherMotifIndex1 + 1)%3;

	const uint AtomCount = SIZE(XVec[MotifIndex][ResidueIndex]);
	const uint ResidueCount1 = SIZE(XVec[OtherMotifIndex1]);
	const uint ResidueCount2 = SIZE(XVec[OtherMotifIndex2]);

	vector<double> MinDists;
	for (uint AtomIndex = 0; AtomIndex < AtomCount; ++AtomIndex)
		{
		const double x = XVec[MotifIndex][ResidueIndex][AtomIndex];
		const double y = YVec[MotifIndex][ResidueIndex][AtomIndex];
		const double z = ZVec[MotifIndex][ResidueIndex][AtomIndex];

		double MinDist1 = GetMinDist(x, y, z,
		  XVec[OtherMotifIndex1], YVec[OtherMotifIndex1], ZVec[OtherMotifIndex1]);

		double MinDist2 = GetMinDist(x, y, z,
		  XVec[OtherMotifIndex2], YVec[OtherMotifIndex2], ZVec[OtherMotifIndex2]);

		double MinDist = min(MinDist1, MinDist2);
		MinDists.push_back(MinDist);
		}

	vector<uint> Order(AtomCount);
	QuickSortOrder(MinDists.data(), AtomCount, Order.data());

	for (uint i = 0; i < K; ++i)
		Ixs.push_back(Order[i]);

#if 1
	{
	for (uint i = 0; i < K; ++i)
		{
		uint Ix = Ixs[i];
		Log("Motif %u i=%u Ix=%u mindist=%.2f\n",
		  MotifIndex, i, Ix, MinDists[Ix]);
		}
	}
#endif
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

	vector<vector<vector<double> > > XVec(3);
	vector<vector<vector<double> > > YVec(3);
	vector<vector<vector<double> > > ZVec(3);
	vector<vector<vector<string> > > ElVec(3);
	vector<vector<vector<string> > > NameVec(3);
	vector<vector<vector<string> > > LineVec(3);
	GetAtomInfos(Q, LoA, HiA, 'A',
	  XVec[0], YVec[0], ZVec[0], ElVec[0], NameVec[0], LineVec[0]);
	GetAtomInfos(Q, LoB, HiB, 'B',
	  XVec[1], YVec[1], ZVec[1], ElVec[1], NameVec[1], LineVec[1]);
	GetAtomInfos(Q, LoC, HiC, 'C',
	  XVec[2], YVec[2], ZVec[2], ElVec[2], NameVec[2], LineVec[2]);

	string PepSeqA;
	string PepSeqB;
	string PepSeqC;

	vector<vector<uint> > IxAVec;
	vector<vector<uint> > IxBVec;
	vector<vector<uint> > IxCVec;
	const uint K = 3;
	for (uint PosA = LoA; PosA <= HiA; ++PosA)
		{
		vector<uint> IxAs;
		uint ResidueIndex = PosA - LoA;
		SelectClosestAtoms(XVec, YVec, ZVec, 0, ResidueIndex, K, IxAs);
		asserta(SIZE(IxAs) == K);
		IxAVec.push_back(IxAs);
		PepSeqA += Q.m_Seq[PosA];
		for (uint k = 0; k < K; ++k)
			{
			uint Ix = IxAs[k];
			const string &Line = LineVec[0][ResidueIndex][Ix];
			if (g_fpdb != 0)
				fprintf(g_fpdb, "%s\n", Line.c_str());
			}
		}

	for (uint PosB = LoB; PosB <= HiB; ++PosB)
		{
		vector<uint> IxBs;
		uint ResidueIndex = PosB - LoB;
		SelectClosestAtoms(XVec, YVec, ZVec, 1, ResidueIndex, K, IxBs);
		asserta(SIZE(IxBs) == K);
		IxAVec.push_back(IxBs);
		PepSeqB += Q.m_Seq[PosB];
		for (uint k = 0; k < K; ++k)
			{
			uint Ix = IxBs[k];
			const string &Line = LineVec[1][ResidueIndex][Ix];
			if (g_fpdb != 0)
				fprintf(g_fpdb, "%s\n", Line.c_str());
			}
		}

	for (uint PosC = LoC; PosC <= HiC; ++PosC)
		{
		vector<uint> IxCs;
		uint ResidueIndex = PosC - LoC;
		SelectClosestAtoms(XVec, YVec, ZVec, 2, ResidueIndex, K, IxCs);
		asserta(SIZE(IxCs) == K);
		IxAVec.push_back(IxCs);
		PepSeqC += Q.m_Seq[PosC];
		for (uint k = 0; k < K; ++k)
			{
			uint Ix = IxCs[k];
			const string &Line = LineVec[2][ResidueIndex][Ix];
			if (g_fpdb != 0)
				fprintf(g_fpdb, "%s\n", Line.c_str());
			}
		}

	if (g_fpml != 0)
		{
		fprintf(g_fpml, "");
		fprintf(g_fpml, "cmd.load(\"%s\")\n", opt_pdb);
		fprintf(g_fpml, "hide\n");
		fprintf(g_fpml, "show spheres\n");
		fprintf(g_fpml, "select pepseq %s\n", PepSeqA.c_str());
		fprintf(g_fpml, "color tv_blue, sele\n");
		fprintf(g_fpml, "select pepseq %s\n", PepSeqB.c_str());
		fprintf(g_fpml, "color tv_green, sele\n");
		fprintf(g_fpml, "select pepseq %s\n", PepSeqC.c_str());
		fprintf(g_fpml, "color tv_red, sele\n");
		fprintf(g_fpml, "deselect\n");
		}
	}

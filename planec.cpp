#include "myutils.h"
#include "pdbchain.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "motifsettings.h"

/***
D:\a\doc\notebooks\2023-07-07_fit_points_to_line_or_plane_3d.txt
***/

static int MinX = -45;
static int MaxX = 34;
static int MinY = -27;
static int MaxY = 41;
static int MinZ = -10;
static int MaxZ = 46;
static int RangeX = MaxX - MinX + 1;
static int RangeY = MaxY - MinY + 1;
static int RangeZ = MaxZ - MinZ + 1;
static int TrainCount;
static vector<vector<vector<uint> > > CountMx;
static vector<vector<vector<double> > > ScoreMx;

static void LogCountMx()
	{
	for (int ix = 0; ix < RangeX; ++ix)
		{
		for (int iy = 0; iy < RangeY; ++iy)
			{
			for (int iz = 0; iz < RangeZ; ++iz)
				{
				uint n = CountMx[ix][iy][iz];
				if (n == 0)
					Log(" ");
				else if (n < 4)
					Log(".", n);
				else if (n < 8)
					Log("+");
				else
					Log("*");
				}
			Log("\n");
			}
		}
	}

static void LogScoreMx()
	{
	for (int ix = 0; ix < RangeX; ++ix)
		{
		for (int iy = 0; iy < RangeY; ++iy)
			{
			for (int iz = 0; iz < RangeZ; ++iz)
				{
				double Score = ScoreMx[ix][iy][iz];
				if (Score == 0)
					Log(" ");
				else if (Score <= 1)
					Log(".");
				else if (Score <= 2)
					Log("+");
				else if (Score <= 4)
					Log("*");
				else
					Log("@");
				}
			Log("\n");
			}
		}
	}

static void ScoreMxToFile(FILE *f)
	{
	if (f == 0)
		return;
	fprintf(f, "scoremx\t%d\t%d\t%d\n", RangeX, RangeY, RangeZ);
	for (int ix = 0; ix < RangeX; ++ix)
		{
		for (int iy = 0; iy < RangeY; ++iy)
			{
			for (int iz = 0; iz < RangeZ; ++iz)
				{
				double Score = ScoreMx[ix][iy][iz];
				if (iz > 0)
					fprintf(f, "\t");
				fprintf(f, "%.3g", Score);
				}
			fprintf(f, "\n");
			}
		}
	}

static double CalcScore(int ix, int iy, int iz)
	{
	int d = 2;
	double Sum = 0;
	for (int jx = ix - d; jx <= ix + d; ++jx)
		{
		if (jx < 0 || jx >= RangeX)
			continue;
		for (int jy = iy - d; jy <= iy + d; ++jy)
			{
			if (jy < 0 || jy >= RangeY)
				continue;
			for (int jz = iz - d; jz <= iz + d; ++jz)
				{
				if (jz < 0 || jz >= RangeZ)
					continue;
				int totd = abs(jx - ix) + abs(jy - iy) + abs(jz - iz);
				if (totd > 3)
					continue;
				uint n = CountMx[jx][jy][jz];
				double f = double(n)/(1 + totd);
				Sum += f;
				}
			}
		}
//	double Score = log2(Sum+1);
	double Score = Sum/TrainCount - 0.01;
	return Score;
	}

static void SetScoreMx()
	{
	ScoreMx.clear();
	ScoreMx.resize(RangeX);
	for (int x = 0; x < RangeX; ++x)
		{
		ScoreMx[x].resize(RangeY);
		for (int y = 0; y < RangeY; ++y)
			ScoreMx[x][y].resize(RangeZ, 0);
		}

	for (int ix = 0; ix < RangeX; ++ix)
		for (int iy = 0; iy < RangeY; ++iy)
			for (int iz = 0; iz < RangeZ; ++iz)
				ScoreMx[ix][iy][iz] = CalcScore(ix, iy, iz);
	}

void WriteN(FILE *f, char Chain, double x, double y, double z)
	{
	static uint WaterAtomNr;
	static uint WaterResNr;

	fprintf(f,
"ATOM  %5u  N   MET %c%4u    %8.3f%8.3f%8.3f  0.00  0.00           O\n",
  WaterAtomNr, Chain, WaterResNr, x, y, z);
	}

static void DrawAxis(FILE *f, const vector<double> &Origin,
  char Chain, const vector<double> &Axis)
	{
	asserta(SIZE(Axis) == 3);
	for (uint i = 0; i < 10; ++i)
		{
		vector<double> Pt;
		MulVecScalar(Axis, i*3, Pt);
		vector<double> Dot;
		Add_Vecs(Origin, Pt, Dot);
		WriteN(f, Chain, Dot[0], Dot[1], Dot[2]);
		}
	}

static void DrawAxes(const string &FileName,
  const vector<double> &Origin,
  const vector<vector<double> > &Basis)
	{
	FILE *f = CreateStdioFile(FileName);
	asserta(SIZE(Basis) == 3);
	for (uint i = 0; i < 3; ++i)
		DrawAxis(f, Origin, "XYZ"[i], Basis[i]);
	CloseStdioFile(f);
	}

void WriteCAPt(FILE *f, char Chain, const vector<double> &Pt)
	{
	WriteN(f, Chain, Pt[0], Pt[1], Pt[2]);
	}

void GetClosestPt(vector<double> &Plane, const vector<double> &Pt,
  vector<double> &ClosestPt)
	{
	asserta(SIZE(Plane) == 3);
	asserta(SIZE(Pt) == 3);

// z = ax + by + c
	const double a = Plane[0];
	const double b = Plane[1];
	const double c = Plane[2];

// Ax + By + Cz + D = 0
	const double A = a;
	const double B = b;
	const double C = -1;
	const double D = c;

	double X = Pt[0];
	double Y = Pt[1];
	double Z = Pt[2];

	double t = - (X*A + Y*B + Z*C + D) / (A*A + B*B + C*C);

	ClosestPt.resize(3);
	ClosestPt[0] = X + t*A;
	ClosestPt[1] = Y + t*B;
	ClosestPt[2] = Z + t*C;
	}

void FitPlanePts(const vector<vector<double> > &Pts,
  vector<double> &Plane)
	{
	const uint N = SIZE(Pts);
	asserta(N > 0);
	double Sumx = 0;
	double Sumy = 0;
	double Sumz = 0;

	double Sumxx = 0;
	double Sumyy = 0;
	double Sumxz = 0;
	double Sumyz = 0;

	double Sumxy = 0;

	for (uint i = 0; i < N; ++i)
		{
		const vector<double> &Pt = Pts[i];
		asserta(SIZE(Pt) == 3);
		double x = Pt[0];
		double y = Pt[1];
		double z = Pt[2];

		Sumx += x;
		Sumy += y;
		Sumz += z;

		Sumxx += x*x;
		Sumyy += y*y;

		Sumxy += x*y;
		Sumxz += x*z;
		Sumyz += y*z;
		}

	vector<vector<double> > Mx;
	Resize3x3(Mx);

	Mx[0][0] = Sumxx;
	Mx[1][1] = Sumyy;
	Mx[2][2] = N;

	Mx[1][0] = Mx[0][1] = Sumxy;
	Mx[2][0] = Mx[0][2] = Sumx;

	Mx[1][2] = Mx[2][1] = Sumy;

	vector<double> Vec(3);
	Vec[0] = Sumxz;
	Vec[1] = Sumyz;
	Vec[2] = Sumz;

	vector<vector<double> > InvMx;
	InvertMx(Mx, InvMx);

	MulMxVec(InvMx, Vec, Plane);

#if 0
	{
	double a = Plane[0];
	double b = Plane[1];
	double c = Plane[2];

// z = ax + by + c
	Log("\n");
	Log("a=%.3g, b=%.3g, c=%.3g\n", a, b, c);

	for (uint i = 0; i < N; ++i)
		{
		const vector<double> &Pt = Pts[i];
		double x = Pt[0];
		double y = Pt[1];
		double z = Pt[2];
		double z2 = a*x + b*y + c;

		vector<double> Ptc;
		GetClosestPt(Plane, Pt, Ptc);

		Log("%10.4f  %10.4f  %10.4f  %10.4f\n",
		  x, y, z, z2);
		Log("%10.4f  %10.4f  %10.4f << closest\n\n",
		  Ptc[0], Ptc[1], Ptc[2]);
		}
	}
#endif
	}

// Ad +0 -3
// Bg +5 
// Cd -2 -5 +2 +5

static void AddFitPt(const PDBChain &Chain, uint Pos,
  vector<vector<double> > &Pts)
	{
	vector<double> Pt;
	Chain.GetPt(Pos, Pt);
	Pts.push_back(Pt);
	}

static void GetFitPts(const PDBChain &Chain,
  uint PosA, uint PosB, uint PosC,
  vector<vector<double> > &Pts)
	{
	Pts.clear();
	AddFitPt(Chain, PosC + g_OffCd +2, Pts);
	AddFitPt(Chain, PosC + g_OffCd +5, Pts);
	AddFitPt(Chain, PosC + g_OffCd -2, Pts);
	AddFitPt(Chain, PosC + g_OffCd -5, Pts);
	AddFitPt(Chain, PosA + g_OffAd -3, Pts);
	AddFitPt(Chain, PosB + g_OffBg +5, Pts);
	}

static void UpdateCountMx1(int xn, int yn, int zn)
	{
	if (xn < 0 || xn >= RangeX)
		return;
	if (yn < 0 || yn >= RangeY)
		return;
	if (zn < 0 || zn >= RangeZ)
		return;
	CountMx[xn][yn][zn] += 1;
	}

static void UpdateCountMx(const PDBChain &XF)
	{
	++TrainCount;
	const uint L = XF.GetSeqLength();
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		double x, y, z;
		XF.GetXYZ(Pos, x, y, z);

		int ix = int(x+0.5);
		int iy = int(y+0.5);
		int iz = int(z+0.5);
		int xn = MaxX - ix;
		int yn = MaxY - iy;
		int zn = MaxZ - iz;
		UpdateCountMx1(xn, yn, zn);
		}
	}

static double GetPPScore(const PDBChain &PP)
	{
	const uint L = PP.GetSeqLength();
	double Sum = 0;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		double x, y, z;
		PP.GetXYZ(Pos, x, y, z);

		int ix = int(x+0.5);
		int iy = int(y+0.5);
		int iz = int(z+0.5);
		int xn = MaxX - ix;
		int yn = MaxY - iy;
		int zn = MaxZ - iz;
		if (xn < 0 || xn >= RangeX)
			continue;
		if (yn < 0 || yn >= RangeY)
			continue;
		if (zn < 0 || zn >= RangeZ)
			continue;
		double Score = ScoreMx[xn][yn][zn];
		Sum += Score;
		}
	double PPScore = Sum/L;
	return PPScore;
	}

void cmd_planec()
	{
	const string &FN = opt_planec;

	RdRpModel Model;
	GetRdrpModel(Model);

	RdRpSearcher RS;
	RS.Init(Model);

	vector<PDBChain *> Chains;
	ReadChains(FN, Chains);
	const uint ChainCount = SIZE(Chains);

//  Transformed pp range of X, Y and Z coords
//  X -44 .. 33
//  Y -26 .. 40
//  Z -9 .. 45
	CountMx.clear();
	CountMx.resize(RangeX);
	for (int x = 0; x < RangeX; ++x)
		{
		CountMx[x].resize(RangeY);
		for (int y = 0; y < RangeY; ++y)
			CountMx[x][y].resize(RangeZ, 0);
		}

	vector<PDBChain *> PPs;
	for (uint i = 0; i < ChainCount; ++i)
		{
		const PDBChain &Chain = *Chains[i];
		const string &Label = Chain.m_Label;
		const string &Seq = Chain.m_Seq;
		RS.Search(Label, Seq);

		string PSSM_A, PSSM_B, PSSM_C;
		RS.GetShapesTrainABC(PSSM_A, PSSM_B, PSSM_C);
		
		uint L = Chain.GetSeqLength();

		size_t stPosA = Seq.find(PSSM_A);
		size_t stPosB = Seq.find(PSSM_B);
		size_t stPosC = Seq.find(PSSM_C);
		if (stPosA == string::npos || stPosB == string::npos ||
		  stPosC == string::npos)
			continue;

		asserta(stPosA < L && stPosB < L && stPosC < L);

		uint PosA = uint(stPosA);
		uint PosB = uint(stPosB);
		uint PosC = uint(stPosC);

		if (PosC < PosA)
			continue;

		PDBChain PP;
		Chain.GetRange(PosA, PosC + g_LC, PP);

		PosB -= PosA;
		PosC -= PosA;
		PosA = 0;

		vector<double> PtC1;
		vector<double> PtC2;
		vector<double> PtA1;
		PP.GetPt(PosC + g_OffCd - 5, PtC1);
		PP.GetPt(PosC + g_OffCd - 2, PtC2);
		PP.GetPt(PosA + g_OffAd, PtA1);

		vector<double> Origin = PtC2;

		vector<double> Axis0;
		Sub_Vecs(PtC1, PtC2, Axis0);
		NormalizeVec(Axis0);

		vector<double> CtoA;
		Sub_Vecs(PtC2, PtA1, CtoA);
		NormalizeVec(CtoA);

		vector<double> Axis2;
		CrossProduct(Axis0, CtoA, Axis2);
		NormalizeVec(Axis2);

		vector<double> Axis1;
		CrossProduct(Axis0, Axis2, Axis1);
		NormalizeVec(Axis1);

		vector<vector<double> > Basis;
		Basis.push_back(Axis0);
		Basis.push_back(Axis1);
		Basis.push_back(Axis2);

		AssertUnitBasisA(Basis);

#if 0
		{
		vector<double> Origin2(3);
		string FileName = Label + "_axes.pdb";

		vector<vector<double> > UnitBasis(3);
		UnitBasis[0].resize(3);
		UnitBasis[1].resize(3);
		UnitBasis[2].resize(3);
		UnitBasis[0][0] = 1;
		UnitBasis[1][1] = 1;
		UnitBasis[2][2] = 1;
		DrawAxes(FileName, Origin2, UnitBasis);
		}
#endif
		vector<vector<double> > R;
		GetBasisR(Basis, R);

		vector<double> t;
		t.push_back(-Origin[0]);
		t.push_back(-Origin[1]);
		t.push_back(-Origin[2]);

		PDBChain &XF = *new PDBChain;
		PP.GetXFormChain_tR(t, R, XF);
		PPs.push_back(&XF);

		UpdateCountMx(XF);
#if 0
		{
		string FileName = Label + "_pp.pdb";
		XF.RenumberResidues(1);
		XF.ToPDB(FileName);
		}
#endif
		}

	SetScoreMx();
	LogScoreMx();

	FILE *f = CreateStdioFile("scoremx.data");
	ScoreMxToFile(f);
	CloseStdioFile(f);

	const uint PPCount = SIZE(PPs);
	for (uint i = 0; i < PPCount; ++i)
		{
		const PDBChain &PP = *PPs[i];
		const PDBChain &Chain = *Chains[i];
		double Score = GetPPScore(PP);
		double NegScore = GetPPScore(Chain);
		Log("Train/neg %10.3g  %10.3g\n", Score, NegScore);
		}

	vector<PDBChain *> ICTVs;
	ReadChains("d:/int/ictv_rdrp_structures/out/ictv.cal", ICTVs);
	const uint M = SIZE(ICTVs);
	for (uint i = 0; i < M; ++i)
		{
		const PDBChain &Chain = *ICTVs[i];

		const string &Label = Chain.m_Label;
		const string &Seq = Chain.m_Seq;
		RS.Search(Label, Seq);

		string PSSM_A, PSSM_B, PSSM_C;
		RS.GetShapesTrainABC(PSSM_A, PSSM_B, PSSM_C);
		
		uint L = Chain.GetSeqLength();

		size_t stPosA = Seq.find(PSSM_A);
		size_t stPosB = Seq.find(PSSM_B);
		size_t stPosC = Seq.find(PSSM_C);
		if (stPosA == string::npos || stPosB == string::npos ||
		  stPosC == string::npos)
			continue;

		asserta(stPosA < L && stPosB < L && stPosC < L);

		uint PosA = uint(stPosA);
		uint PosB = uint(stPosB);
		uint PosC = uint(stPosC);

		if (PosC < PosA)
			continue;

		PDBChain PP;
		Chain.GetRange(PosA, PosC + g_LC, PP);

		PosB -= PosA;
		PosC -= PosA;
		PosA = 0;

		vector<double> PtC1;
		vector<double> PtC2;
		vector<double> PtA1;
		PP.GetPt(PosC + g_OffCd - 5, PtC1);
		PP.GetPt(PosC + g_OffCd - 2, PtC2);
		PP.GetPt(PosA + g_OffAd, PtA1);

		vector<double> Origin = PtC2;

		vector<double> Axis0;
		Sub_Vecs(PtC1, PtC2, Axis0);
		NormalizeVec(Axis0);

		vector<double> CtoA;
		Sub_Vecs(PtC2, PtA1, CtoA);
		NormalizeVec(CtoA);

		vector<double> Axis2;
		CrossProduct(Axis0, CtoA, Axis2);
		NormalizeVec(Axis2);

		vector<double> Axis1;
		CrossProduct(Axis0, Axis2, Axis1);
		NormalizeVec(Axis1);

		vector<vector<double> > Basis;
		Basis.push_back(Axis0);
		Basis.push_back(Axis1);
		Basis.push_back(Axis2);

		AssertUnitBasisA(Basis);

#if 0
		{
		vector<double> Origin2(3);
		string FileName = Label + "_axes.pdb";

		vector<vector<double> > UnitBasis(3);
		UnitBasis[0].resize(3);
		UnitBasis[1].resize(3);
		UnitBasis[2].resize(3);
		UnitBasis[0][0] = 1;
		UnitBasis[1][1] = 1;
		UnitBasis[2][2] = 1;
		DrawAxes(FileName, Origin2, UnitBasis);
		}
#endif
		vector<vector<double> > R;
		GetBasisR(Basis, R);

		vector<double> t;
		t.push_back(-Origin[0]);
		t.push_back(-Origin[1]);
		t.push_back(-Origin[2]);

		PDBChain &XF = *new PDBChain;
		PP.GetXFormChain_tR(t, R, XF);
		double Score = GetPPScore(PP);

		Log("ICTV %.3g\n", Score);
		}
	}

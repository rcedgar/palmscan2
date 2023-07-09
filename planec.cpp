#include "myutils.h"
#include "pdbchain.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "motifsettings.h"

void FitPlanePts(const vector<vector<double> > &Pts)
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

	vector<double> Plane(3);
	MulMxVec(InvMx, Vec, Plane);

	double a = Plane[0];
	double b = Plane[1];
	double c = Plane[2];

// z = ax + by + c
	Log("\n");
	Log("a=%.3g, b=%.3g, c=%.3g\n", a, b, c);//@@

	for (uint i = 0; i < N; ++i)
		{
		const vector<double> &Pt = Pts[i];
		double x = Pt[0];
		double y = Pt[1];
		double z = Pt[2];
		double z2 = a*x + b*y + c;

		Log("%10.4g  %10.4g  %10.4g  %10.4g\n",
		  x, y, z, z2);
		}
	}

void FitPlaneC(const PDBChain &Chain, uint PosC, uint CL)
	{
	vector<vector<double> > Pts;

	for (uint i = 0; i < CL; ++i)
		{
		vector<double> Pt;
		Chain.GetPt(PosC+i, Pt);
		Pts.push_back(Pt);
		}

	FitPlanePts(Pts);
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

		uint PosC = uint(stPosC);
		FitPlaneC(*Chains[i], PosC, g_LC);
		}
	}

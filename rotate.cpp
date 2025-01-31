#include "myutils.h"
#include "pdbchain.h"

static void RotatePtAxis(const double InPt[3],
						 uint Axis,
						 double theta,
						 double OutPt[3])
	{
	assert(Axis < 3);
	const double st = sin(theta);
	const double  ct = cos(theta);

	double Mx[3][3];

	Mx[Axis%3][Axis%3] = 1;
	Mx[Axis%3][(Axis+1)%3] = 0;
	Mx[Axis%3][(Axis+2)%3] = 0;

	Mx[(Axis+1)%3][Axis%3] = 0;
	Mx[(Axis+1)%3][(Axis+1)%3] = ct;
	Mx[(Axis+1)%3][(Axis+2)%3] = -st;

	Mx[(Axis+2)%3][Axis%3] = 0;
	Mx[(Axis+2)%3][(Axis+1)%3] = st;
	Mx[(Axis+2)%3][(Axis+2)%3] = ct;

	OutPt[0] = Mx[0][0]*InPt[0] + Mx[0][1]*InPt[1] + Mx[0][2]*InPt[2];
	OutPt[1] = Mx[1][0]*InPt[0] + Mx[1][1]*InPt[1] + Mx[1][2]*InPt[2];
	OutPt[2] = Mx[2][0]*InPt[0] + Mx[2][1]*InPt[1] + Mx[2][2]*InPt[2];
	}

static void RotatePt(const double PtIn[3],
				   double alpha, double beta, double gamma,
				   double PtOut[3])
	{
	double TmpPtX[3];
	double TmpPtY[3];
	RotatePtAxis(PtIn, 0, alpha, TmpPtX);
	RotatePtAxis(TmpPtX, 1, beta, TmpPtY);
	RotatePtAxis(TmpPtY, 2, gamma, PtOut);
	}

void RotateChain(PDBChain &Chain, double alpha, double beta, double gamma)
	{
	uint L = Chain.GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		double Pt[3];
		Pt[0] = Chain.m_Xs[i];
		Pt[1] = Chain.m_Ys[i];
		Pt[2] = Chain.m_Zs[i];
		
		double PtOut[3];
		RotatePt(Pt, alpha, beta, gamma, PtOut);

		Chain.m_Xs[i] = PtOut[0];
		Chain.m_Ys[i] = PtOut[1];
		Chain.m_Zs[i] = PtOut[2];
		}
	const uint n = SIZE(Chain.m_ATOMs);
	if (n == 0)
		return;
	asserta(n == L);
	asserta(SIZE(Chain.m_ATOMs) == L);
	for (uint i = 0; i < L; ++i)
		{
		vector<string> &ResATOMs = Chain.m_ATOMs[i];
		const uint n = SIZE(ResATOMs);
		for (uint j = 0; j < n; ++j)
			{
			string &Line = ResATOMs[j];
			double x, y, z;
			PDBChain::GetXYZFromATOMLine(Line, x, y, z);

			double Pt[3];
			Pt[0] = x;
			Pt[1] = y;
			Pt[2] = z;
		
			double PtOut[3];
			RotatePt(Pt, alpha, beta, gamma, PtOut);

			x = PtOut[0];
			y = PtOut[1];
			z = PtOut[2];

			string UpdatedLine;
			PDBChain::SetXYZInATOMLine(Line, x, y, z, UpdatedLine);
			Chain.m_ATOMs[i][j] = UpdatedLine;
			}
		}
	}

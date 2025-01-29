#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"

void GetThreeFromOne(char aa, string &AAA);

static void WritePDB(const PDBChain &Chain, const string &FileName)
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	const uint L = Chain.GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		char aa = Chain.m_Seq[i];
		string sAAA;
		GetThreeFromOne(aa, sAAA);
		const char *AAA = sAAA.c_str();

		fprintf(f, "ATOM  ");			//  1 -  6        Record name   "ATOM  "
		fprintf(f, "%5u", i+1);			//  7 - 11        Integer       serial       Atom serial number.
		fprintf(f, " ");				// 12
		fprintf(f, " CA ");				// 13 - 16        Atom          name         Atom name.
		fprintf(f, " ");				// 17             Character     altLoc       Alternate location indicator.
		fprintf(f, "%3.3s", AAA);		// 18 - 20        Residue name  resName      Residue name.
		fprintf(f, " ");				// 21
		fprintf(f, "%c", 'A');			// 22             Character     chainID      Chain identifier.
		fprintf(f, "%4u", i+1);			// 23 - 26        Integer       resSeq       Residue sequence number.
		fprintf(f, " ");				// 27             AChar         iCode        Code for insertion of residues.
		fprintf(f, "   ");				// 28 - 30
		fprintf(f, "%8.3f", Chain.m_Xs[i]);	// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		fprintf(f, "%8.3f", Chain.m_Ys[i]);	// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		fprintf(f, "%8.3f", Chain.m_Zs[i]);	// 47 - 57        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		fprintf(f, "%6.2f", 1.0);		// 55 - 60        Real(6.2)     occupancy    Occupancy.
		fprintf(f, "%6.2f", 0.0);		// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		fprintf(f, "          ");		// 67 - 76
		fprintf(f, " C");				// 77 - 78        LString(2)    element      Element symbol, right-justified.
		fprintf(f, "  ");				// 79 - 80        LString(2)    charge       Charge on the atom.

		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}

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

static void RotateChain(PDBChain &Chain, double alpha, double beta, double gamma)
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
	}

static void Join(const PDBChain &a, const PDBChain &b, PDBChain &ab, const string &Label)
	{
	PDBChain A = a;
	PDBChain B = b;
	const uint LA = A.GetSeqLength();
	const uint LB = B.GetSeqLength();

	double xn = A.m_Xs[LA-1];
	double yn = A.m_Ys[LA-1];
	double zn = A.m_Zs[LA-1];
	for (uint i = 0; i < LA; ++i)
		{
		A.m_Xs[i] -= xn;
		A.m_Ys[i] -= yn;
		A.m_Zs[i] -= zn;
		}

	WritePDB(A, "a.pdb");

	double alpha = (randu32()%1000)/200.0;
	double beta = (randu32()%1000)/200.0;
	double gamma = (randu32()%1000)/200.0;
	RotateChain(B, alpha, beta, gamma);

	double x0 = B.m_Xs[0];
	double y0 = B.m_Ys[0];
	double z0 = B.m_Zs[0];
	for (uint i = 0; i < LB; ++i)
		{
		B.m_Xs[i] -= x0;
		B.m_Ys[i] -= y0;
		B.m_Zs[i] -= z0;
		}

	WritePDB(B, "b.pdb");

	//RotateChain(B, 0, 0, 0);
	//WritePDB(B, "b0.pdb");

	//RotateChain(B, 0, -3.14159/4, 0);
	//WritePDB(B, "b2.pdb");

	//RotateChain(B, 0, 3.14159/2, 0);
	//WritePDB(B, "b3.pdb");

	ab.m_Label = Label;
	ab.m_Seq = A.m_Seq + B.m_Seq;
	ab.m_Xs = A.m_Xs;
	ab.m_Ys = A.m_Ys;
	ab.m_Zs = A.m_Zs;

	ab.m_Xs.insert(ab.m_Xs.end(), B.m_Xs.begin(), B.m_Xs.end());
	ab.m_Ys.insert(ab.m_Ys.end(), B.m_Ys.begin(), B.m_Ys.end());
	ab.m_Zs.insert(ab.m_Zs.end(), B.m_Zs.begin(), B.m_Zs.end());
	WritePDB(ab, Label + ".pdb");
	}

static void TestRot()
	{
	double Pt[3] = { 1, 2, 3 };
	double OutPt[3];
	RotatePtAxis(Pt, 2, 3.14159/2, OutPt);
	Log("%.3g %.3g %.3g\n", OutPt[0], OutPt[1], OutPt[2]);
	}

void cmd_fake()
	{
	const string &FN = opt_fake;
	//TestRot();
	//return;

	vector<PDBChain *> Frags;
	ReadChains(FN, Frags);

	PDBChain ab;
	Join(*Frags[0], *Frags[1], ab, "ab");

	PDBChain abc;
	Join(ab, *Frags[2], abc, "abc");

	PDBChain abcd;
	Join(abc, *Frags[3], abcd, "abcd");
	}

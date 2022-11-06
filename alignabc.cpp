#include "myutils.h"
#include "pdb.h"

void cmd_alignabc()
	{
	const string &QueryFN = opt_alignabc;
	const string &RefFN = opt_ref;

	PDB Q;
	PDB R;

	Q.FromFile(QueryFN);
	Q.LogMe();

	R.FromFile(RefFN);

	R.m_MotifPosVec.push_back(0 + 12/2);
	R.m_MotifPosVec.push_back(65 + 12/2);
	R.m_MotifPosVec.push_back(98 + 8/2);

	R.LogMe();

	double LAB, LBC, LAC;
	R.GetMotifTriangle(LAB, LBC, LAC);
	Log("LAB=%.1f, LBC=%.1f, LAC=%.1f\n", LAB, LBC, LAC);

// LAB=14.1, LBC=10.4, LAC=11.4

	TriParams TP;
	TP.LAB = 14.1;
	TP.LBC = 10.4;
	TP.LAC = 11.4;
	TP.Radius = 1;
	TP.NABmin = 10;
	TP.NABmax = 80;
	TP.NBCmin = 10;
	TP.NBCmax = 80;
	TP.NACmin = 80;
	TP.NACmax = 200;

	// LAB=14.1, LBC=10.4, LAC=11.4
	vector<uint> PosAs;
	vector<uint> PosBs;
	vector<uint> PosCs;
	vector<double> RMSDs;
	Q.SearchTriangle(TP, PosAs, PosBs, PosCs, RMSDs);

	const uint N = SIZE(RMSDs);
	Log("%u hits\n", N);
	for (uint i = 0; i < N; ++i)
		{
		uint PosA = PosAs[i];
		uint PosB = PosBs[i];
		uint PosC = PosCs[i];
		double RMSD = RMSDs[i];

		string A, B, C;
		Q.GetMotifSeqFromMidPos(PosA, 12, false, A);
		Q.GetMotifSeqFromMidPos(PosB, 12, false, B);
		Q.GetMotifSeqFromMidPos(PosC, 8, false, C);

		Log("%5u", PosA);
		Log("  %5u", PosB);
		Log("  %5u", PosC);
		Log("  %7.2f", RMSD);
		Log("  %s", A.c_str());
		Log("  %s", B.c_str());
		Log("  %s", C.c_str());
		Log("\n");
		}
	}

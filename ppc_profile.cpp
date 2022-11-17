#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "abcxyz.h"
#include "quarts.h"
#include "outputfiles.h"

void cmd_ppc_profile()
	{
	const string &InputFileName = opt_ppc_profile;

	ChainReader CR;
	CR.Open(InputFileName, false);

	Pf(g_ftsv, "Label");
	Pf(g_ftsv, "\tDGD");
	Pf(g_ftsv, "\tPosA");
	Pf(g_ftsv, "\tPosB");
	Pf(g_ftsv, "\tPosC");
	Pf(g_ftsv, "\tLAB");
	Pf(g_ftsv, "\tLAC");
	Pf(g_ftsv, "\tLBC");
	Pf(g_ftsv, "\tdAB");
	Pf(g_ftsv, "\tdAC");
	Pf(g_ftsv, "\tdBC");
	Pf(g_ftsv, "\n");

	vector<double> VecLAB;
	vector<double> VecLAC;
	vector<double> VecLBC;

	vector<double> VecdAB;
	vector<double> VecdAC;
	vector<double> VecdBC;

	uint ChainCount = 0;
	PDBChain Chain;
	while (CR.GetNext(Chain))
		{
		uint SeqLength = Chain.GetSeqLength();
		if (++ChainCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% scanning profiles\r", sPct.c_str());
			}
		Chain.CheckPPCMotifCoords(true);

		uint PosA = Chain.m_MotifPosVec[A];
		uint PosB = Chain.m_MotifPosVec[B];
		uint PosC = Chain.m_MotifPosVec[C];

		char DGD[4];
		DGD[0] = Chain.m_Seq[PosA];
		DGD[1] = Chain.m_Seq[PosB];
		DGD[2] = Chain.m_Seq[PosC];
		DGD[3] = 0;

		double dAB = Chain.GetDist(PosA, PosB);
		double dAC = Chain.GetDist(PosA, PosC);
		double dBC = Chain.GetDist(PosB, PosC);

		uint LAB = PosB - PosA + 1;
		uint LAC = PosC - PosA + 1;
		uint LBC = PosC - PosB + 1;

		VecLAB.push_back(LAB);
		VecLAC.push_back(LAC);
		VecLBC.push_back(LBC);

		VecdAB.push_back(dAB);
		VecdAC.push_back(dAC);
		VecdBC.push_back(dBC);

		Pf(g_ftsv, "%s", Chain.m_Label.c_str());
		Pf(g_ftsv, "\t%s", DGD);
		Pf(g_ftsv, "\t%u", PosA + 1);
		Pf(g_ftsv, "\t%u", PosB + 1);
		Pf(g_ftsv, "\t%u", PosC + 1);
		Pf(g_ftsv, "\t%u", LAB);
		Pf(g_ftsv, "\t%u", LAC);
		Pf(g_ftsv, "\t%u", LBC);
		Pf(g_ftsv, "\t%.2f", dAB);
		Pf(g_ftsv, "\t%.2f", dAC);
		Pf(g_ftsv, "\t%.2f", dBC);
		Pf(g_ftsv, "\n");
		}

#define X(x)	\
	Log("\n");	\
	QuartsDouble QF_L##x;	\
	QuartsDouble QF_d##x;	\
	GetQuartsDouble(VecL##x, QF_L##x); \
	GetQuartsDouble(Vecd##x, QF_d##x); \
	double MeanL_##x = QF_L##x.Avg; \
	double StdDevL_##x = QF_L##x.StdDev; \
	double Meand_##x = QF_d##x.Avg; \
	double StdDevd_##x = QF_d##x.StdDev; \
	Log("double Mean_L_%s = %.1f;\n", #x, MeanL_##x); \
	Log("double StdDev_L_%s = %.1f;\n", #x, StdDevL_##x); \
	Log("double Mean_d_%s = %.1f;\n", #x, Meand_##x); \
	Log("double StdDev_d_%s = %.1f;\n", #x, StdDevd_##x); \

	X(AB)
	X(AC)
	X(BC)

	ProgressLog("%u chains\n", ChainCount);
	}

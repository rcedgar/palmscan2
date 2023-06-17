#include "myutils.h"
#include "structprof.h"
#include "chainreader.h"
#include "outputfiles.h"
#include "cmpsearcher.h"
#include "gsprof.h"
#include "abcxyz.h"

uint g_APos;
uint g_BPos;
uint g_CPos;
uint g_DPos;
uint g_EPos;
uint g_F1Pos;
uint g_F2Pos;

void WritePML(const PDBChain &Chain, const string &PDBFileName);

static void WriteSubChain(const PDBChain &Chain)
	{
	if (!optset_subchain)
		return;

	if (g_APos == UINT_MAX || g_CPos == UINT_MAX)
		return;

	uint LoPos = 0;
	uint HiPos = min(g_APos, g_CPos);
	if (optset_subchainc)
		{
		LoPos = max(g_APos, g_CPos);
		HiPos = Chain.GetSeqLength() - 1;
		}

	PDBChain SubChain;
	Chain.GetRange(LoPos, HiPos, SubChain);
	SubChain.ToPDB(opt_subchain);
	exit(0);
	}

static void WriteStructProfPos(FILE *f, const StructProf &SP, uint Pos)
	{
	if (f == 0)
		return;

	const PDBChain &Chain = *SP.m_Chain;
	const char *Label = Chain.m_Label.c_str();
	const char *Seq = Chain.m_Seq.c_str();
	char aa = Seq[Pos];
	char ss = Chain.m_SS[Pos];
	uint NU, ND;
	SP.GetHSE(Pos, 12.0, NU, ND);
	uint TSB = SP.GetTSB(Pos, 10.0);
	uint CN = SP.GetCavityNumber(Pos);

	uint Pos_aD = Chain.GetMotifPos(A) + 3;
	uint Pos_bG = Chain.GetMotifPos(B) + 1;
	uint Pos_cD = Chain.GetMotifPos(C) + 3;
	asserta(Seq[Pos_aD] == 'D');
	asserta(Seq[Pos_bG] == 'G');
	asserta(Seq[Pos_cD] == 'D');

	double Dist_aD = Chain.GetDist(Pos, Pos_aD);
	double Dist_bG = Chain.GetDist(Pos, Pos_bG);
	double Dist_cD = Chain.GetDist(Pos, Pos_cD);

	static bool HdrDone = false;
	if (!HdrDone)
		{
		fprintf(f, "Label");
		fprintf(f, "\tPos");
		fprintf(f, "\taa");
		fprintf(f, "\tss");
		fprintf(f, "\tNU");
		fprintf(f, "\tND");
		fprintf(f, "\tTSB");
		fprintf(f, "\tDist_aD");
		fprintf(f, "\tDist_bG");
		fprintf(f, "\tDist_cD");
		fprintf(f, "\tCN");
		fprintf(f, "\tMotif");
		fprintf(f, "\n");
		HdrDone = true;
		}

	int ResNr = SP.m_Chain->GetResidueNr(Pos);
	char Motif = '.';
	if (Pos == g_APos)
		Motif = 'A';
	if (Pos == g_BPos)
		Motif = 'B';
	if (Pos == g_CPos)
		Motif = 'C';
	if (Pos == g_DPos)
		Motif = 'D';
	if (Pos == g_EPos)
		Motif = 'E';
	if (Pos == g_F2Pos)
		Motif = 'F';

	fprintf(f, "%s", Label);
	fprintf(f, "\t%u", Pos+1);
	fprintf(f, "\t%c", aa);
	fprintf(f, "\t%c", ss);
	fprintf(f, "\t%u", NU);
	fprintf(f, "\t%u", ND);
	fprintf(f, "\t%u", TSB);
	fprintf(f, "\t%.1f", Dist_aD);
	fprintf(f, "\t%.1f", Dist_bG);
	fprintf(f, "\t%.1f", Dist_cD);
	fprintf(f, "\t%u", CN);
	fprintf(f, "\t%c", Motif);
	if (ResNr == INT_MAX)
		fprintf(f, "\t.");
	else
		fprintf(f, "\t%d", ResNr);
	fprintf(f, "\n");
	}

static void WriteStructProfPosTsv(const StructProf &SP, uint MinPos, uint MaxPos)
	{
	if (g_ftsv == 0)
		return;
	for (uint Pos = MinPos; Pos <= MaxPos; ++Pos)
		WriteStructProfPos(g_ftsv, SP, Pos);
	}

static bool DoSearch3D(const string &PDBFileName, CMPSearcher &CS,  PDBChain &Chain)
	{
	const uint L = Chain.GetSeqLength();
	uint MinPos = (g_APos < 150 ? 0 : g_APos - 150);
	uint MaxPos = g_CPos + CL + 150 - 1;
	if (MaxPos >= L)
		MaxPos = L - 1;

	StructProf SP;
	SP.SetChain(Chain);
	SP.SetCavityCenterPt();

	g_DPos = SP.FindMofifD_Hueuristics();
	g_EPos = SP.FindMofifE_Hueuristics(g_DPos);
	g_F2Pos = SP.FindMofifF2_Hueuristics(g_APos);

	SP.WriteGSProf(g_fgsprof_out);
	SP.WriteMotifsTsv(g_fmotifs);
	WriteStructProfPosTsv(SP, MinPos, MaxPos);
	WritePML(Chain, PDBFileName);
	WriteSubChain(Chain);

	return true;
	}

static void Search3D(const string &InputFN)
	{
	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;

	CMP Prof;
	Prof.FromFile(ModelFileName);

	CMPSearcher CS;
	CS.SetProf(Prof);

	ChainReader CR;
	CR.Open(InputFN);

	PDBChain Chain;
	PDBChain XChain;
	ProgressStep(0, 1001, "Processing");
	uint LastMil = 0;
	uint DoneCount = 0;
	while (CR.GetNext(Chain))
		{
		uint Mil = CR.GetMilDone();
		if (Mil > 0)
			ProgressStep(Mil, 1001, "Processing %u", DoneCount);

		CS.Search(Chain);

		g_APos = UINT_MAX;
		g_BPos = UINT_MAX;
		g_CPos = UINT_MAX;
		double PalmScore = CS.GetPSSMStarts(g_APos, g_BPos, g_CPos);
		if (PalmScore <= 0)
			continue;

		Chain.SetMotifPosVec(g_APos, g_BPos, g_CPos);
		bool Ok = DoSearch3D(InputFN, CS, Chain);
		++DoneCount;
		if (opt_first_only && Ok)
			break;
		}
	ProgressStep(1000, 1001, "Processing %u", DoneCount);
	}

// for backwards compatibility
void cmd_struct_prof()
	{
	const string &InputFN = opt_struct_prof;
	Search3D(InputFN);
	}

void cmd_search3d()
	{
	const string &InputFN = opt_search3d;
	Search3D(InputFN);
	}

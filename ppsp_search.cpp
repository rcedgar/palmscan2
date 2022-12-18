#include "myutils.h"
#include "pdbchain.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "ppspsearcher.h"

static uint g_DoneCount;
static uint g_HitCount;

void cmd_ppsp_search()
	{
	const string &QueryFN = opt_ppsp_search;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;

	PPSPSearcher PS;
	PS.m_Prof.FromFile(ModelFileName);

	PDBChain Q;

	ChainReader CR;
	CR.m_SaveAtoms = true;
	CR.Open(QueryFN, false);

	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

		if (++g_DoneCount%1000 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}

		const char *Seq = Q.m_Seq.c_str();
		const string &Label = Q.m_Label;
		uint APos = UINT_MAX;
		uint BPos = UINT_MAX;
		uint CPos = UINT_MAX;
		PS.Search(Q);
		double Score = PS.GetPSSMStarts(APos, BPos, CPos);
		if (APos == UINT_MAX)
			continue;
		++g_HitCount;

		if (g_ftsv != 0)
			{
			const char *SeqA = Seq + APos;
			const char *SeqB = Seq + BPos;
			const char *SeqC = Seq + CPos;

			char Gate = SeqA[8];
			const char *GDD = SeqC + 2;

			const char *Cat = "RdRp";
			if (Gate == 'F' || Gate == 'Y')
				Cat = "RT";
			if (!strncmp(GDD, "ADD", 3) ||
			    !strncmp(GDD, "VDD", 3) ||
				!strncmp(GDD, "LDD", 3) ||
				!strncmp(GDD, "MDD", 3))
				Cat = "RT";

			fprintf(g_ftsv, "%.4f", Score);
			fprintf(g_ftsv, "\t%s", Cat);
			fprintf(g_ftsv, "\t%s", Label.c_str());
			fprintf(g_ftsv, "\t%u", APos+1);
			fprintf(g_ftsv, "\t%u", BPos+1);
			fprintf(g_ftsv, "\t%u", CPos+1);
			fprintf(g_ftsv, "\t%*.*s", AL, AL, SeqA);
			fprintf(g_ftsv, "\t%*.*s", BL, BL, SeqB);
			fprintf(g_ftsv, "\t%*.*s", CL, CL, SeqC);
			fprintf(g_ftsv, "\t%c", Gate);
			fprintf(g_ftsv, "\t%3.3s", GDD);
			fprintf(g_ftsv, "\n");
			}
		}

	Progress("100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);
	}

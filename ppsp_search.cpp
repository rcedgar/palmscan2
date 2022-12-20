#include "myutils.h"
#include "pdbchain.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "ppspsearcher.h"
#include <time.h>

static uint g_DoneCount;
static uint g_HitCount;

static void Thread(ChainReader &CR, const PPSP &Prof)
	{
	PPSPSearcher PS;
	PDBChain Q;
	PS.m_Prof = Prof;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

#pragma omp critical
		{
		if (++g_DoneCount%1000 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}
		}

		const char *Seq = Q.m_Seq.c_str();
		const string &Label = Q.m_Label;
		uint APos = UINT_MAX;
		uint BPos = UINT_MAX;
		uint CPos = UINT_MAX;
		PS.Search(Q);
		double PalmScore = PS.GetPSSMStarts(APos, BPos, CPos);
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

			char CmfX = SeqC[2];

			double P_rdrp = PPSP::GetRdRpProb(Gate, CmfX);
			double P_rdrp_gate = PPSP::GetRdRpProb_Gate(Gate);
			double P_rdrp_cmfx = PPSP::GetRdRpProb_CmfX(CmfX);

			double AdjustedPalmScore = 0.75 + PalmScore/4;
			double RdRpScore = AdjustedPalmScore*P_rdrp;

			fprintf(g_ftsv, "%.4f", RdRpScore);
			fprintf(g_ftsv, "\t%.4f", PalmScore);
			fprintf(g_ftsv, "\t%.4f", P_rdrp);
			fprintf(g_ftsv, "\t%s", Label.c_str());
			fprintf(g_ftsv, "\t%u", APos+1);
			fprintf(g_ftsv, "\t%u", BPos+1);
			fprintf(g_ftsv, "\t%u", CPos+1);
			fprintf(g_ftsv, "\t%*.*s", AL, AL, SeqA);
			fprintf(g_ftsv, "\t%*.*s", BL, BL, SeqB);
			fprintf(g_ftsv, "\t%*.*s", CL, CL, SeqC);
			fprintf(g_ftsv, "\t%c(%.4f)", Gate, P_rdrp_gate);
			fprintf(g_ftsv, "\t%3.3s(%.4f)", GDD, P_rdrp_cmfx);
			fprintf(g_ftsv, "\n");
			}
		}
	}

void cmd_ppsp_search()
	{
	const string &QueryFN = opt_ppsp_search;

	time_t tStart = time(0);
	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;

	PPSP Prof;
	Prof.FromFile(ModelFileName);

	ChainReader CR;
	CR.m_SaveAtoms = true;
	CR.Open(QueryFN, false);

	uint ThreadCount = GetRequestedThreadCount();

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, Prof);
	time_t tEnd = time(0);
	uint Secs = uint(tEnd - tStart);
	if (Secs <= 0)
		Secs = 1;
	double Throughput = double(g_DoneCount)/(Secs*ThreadCount);
	ProgressLog("%u done, %u hits, %s secs (%u threads, %.1f/ sec/ thread)\n",
	  g_DoneCount, g_HitCount, IntToStr(Secs), ThreadCount, Throughput);
	}

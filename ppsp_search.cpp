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

		if (++g_DoneCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}

		const string &QSeq = Q.m_Seq;
		const string &QLabel = Q.m_Label;
		PS.Search(Q);
		Log("\n");
		Log(">%s\n", Q.m_Label.c_str());
		const uint HitCount = SIZE(PS.m_Ads);
		Log("%u hits\n", HitCount);
		for (uint i = 0; i < HitCount; ++i)
			{
			Log("%6.3f  %4u  %4u  %4u\n",
			  PS.m_Scores[i], PS.m_Ads[i], PS.m_Bgs[i], PS.m_Cds[i]);
			}
		}

	Progress("100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);
	}

#include "myutils.h"
#include "tma.h"
#include "outputfiles.h"
#include "chainreader.h"

static double g_MinTM = 0.6;
static uint g_QueryCount = 0;
static uint g_HitCount = 0;

void TMAlignPair(TMA &T, const PDBChain &Q, const PDBChain &R);

static void WriteHit(const PDBChain &Q, const PDBChain &R, double TM)
	{
	if (g_ftsv == 0)
		return;
	const char *LabelQ = Q.m_Label.c_str();
	const char *LabelR = R.m_Label.c_str();
#pragma omp critical
	{
	fprintf(g_ftsv, "%s\t%s\t%.4f\n", LabelQ, LabelR, TM);
	}
	}

static void Thread(ChainReader &CR, const vector<PDBChain *> &DBChains)
	{
	const uint DBSize = SIZE(DBChains);
	TMA T;
	uint ThreadIndex = GetThreadIndex();

	PDBChain Q;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;
		string tmps;
		Log("%16s  [%2u] START >%s\n",
		  GetElapsedTimeStr(tmps), ThreadIndex, Q.m_Label.c_str());
		uint t1 = GetElapsedSecs();
		bool AnyHits = false;
		for (uint dbidx = 0; dbidx < DBSize; ++dbidx)
			{
			const PDBChain &R = *DBChains[dbidx];
			double TM = T.AlignChains(Q, R);
			if (TM >= g_MinTM)
				{
				AnyHits = true;
				WriteHit(Q, R, TM);
				}
			}
		uint t2 = GetElapsedSecs();
		Log("%16s  [%2u] END %u secs >%s\n",
		  GetElapsedTimeStr(tmps), ThreadIndex, t2 - t1, Q.m_Label.c_str());
#pragma omp critical
		{
		++g_QueryCount;
		if (AnyHits)
			++g_HitCount;
		if (1 || g_QueryCount%10 == 0)
			{
			double PctDone = CR.GetPctDone();
			uint Pct = uint(PctDone);
			if (Pct == 0)
				Pct = 1;
			if (Pct >= 99)
				Pct = 98;
			Progress("%u / %u hits\r", g_HitCount, g_QueryCount);
			}
		}
		}
	}

void cmd_tm()
	{
	const string &QueryFileName = opt_tm;
	const string &DBFileName = opt_db;

	if (optset_mintm)
		g_MinTM = opt_mintm;
	asserta(g_MinTM >= 0 && g_MinTM <= 1);

	vector<PDBChain *> DBChains;
	ReadChains(DBFileName, DBChains);

	ChainReader CR;
	CR.Open(QueryFileName);

	uint ThreadCount = GetRequestedThreadCount();
	Progress("Aligning\n");
#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, DBChains);
	Progress("%u / %u hits\n", g_HitCount, g_QueryCount);

	CR.Clear();
	}

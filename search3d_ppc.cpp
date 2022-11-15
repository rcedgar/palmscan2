#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "sort.h"
#include "calreader.h"
#include "ppcsearcher.h"
#include "searchparams.h"
#include "sort.h"
#include "omplock.h"
#include "outputfiles.h"

void ReadPpc(const string &FN, vector<PDBChain *> &Chains);

static uint g_DoneCount;
static uint g_HitCount;

static void Thread(CalReader &CR, vector<PDBChain> &Qs,
  vector<PPCSearcher> &PSs)
	{
	uint ThreadIndex = GetThreadIndex();
	PDBChain &Q = Qs[ThreadIndex];
	PPCSearcher &PS = PSs[ThreadIndex];
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			return;
		Lock("Done");
		++g_DoneCount;
		Unlock("Done");
		if (g_DoneCount%100 == 0)
			{
			Lock("Progress");
			string sPct;
			CR.GetPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			Unlock("Progress");
			}

		bool Found = PS.Search(Q);
		PS.WriteOutput();
		if (Found)
			++g_HitCount;
		}
	}

void cmd_search3d_ppc()
	{
	const string &QueryFN = opt_search3d_ppc;
	const string &RefFN = opt_ref;

	CalReader CR;
	CR.Open(QueryFN);

	vector<PDBChain *> RefPDBs;
	ReadPpc(RefFN, RefPDBs);
	for (uint i = 0; i < SIZE(RefPDBs); ++i)
		{
		PDBChain &R = *RefPDBs[i];
		R.CheckPPCMotifCoords();
		R.SetSS();
		}

	const uint RefN = SIZE(RefPDBs);

	uint ThreadCount = GetRequestedThreadCount();
	vector<vector<PDBChain *> > QVecs(ThreadCount);
	vector<PPCSearcher> PSs(ThreadCount);
	for (uint i = 0; i < ThreadCount; ++i)
		{
		PPCSearcher &PS = PSs[i];
		PS.SetRefs(RefPDBs);
		PS.InitClassify();
		}

	uint ThreadFinishedCount = 0;
	uint DoneCount = 0;
	uint HitCount = 0;
	vector<PDBChain> Qs(ThreadCount);
	vector<bool> ThreadDone(ThreadCount);

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, Qs, PSs);

	Progress("100.0% done, %u / %u hits\n", DoneCount, HitCount);
	}

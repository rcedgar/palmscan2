#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "fakechain.h"
#include "omplock.h"
#include <time.h>

void cmd_fake()
	{
	const string &InputFN = g_Arg1;
	FILE *fCal = CreateStdioFile(opt_cal);
	FILE *fTsv = CreateStdioFile(opt_tsv);

	vector<PDBChain *> Frags;
	ReadChains(InputFN, Frags);

	const uint FragCount = SIZE(Frags);
	for (uint i = 0; i < FragCount; ++i)
		{
		ProgressStep(i, FragCount, "Shuffle sequences");
		PDBChain &Chain = *Frags[i];
		random_shuffle(Chain.m_Seq.begin(), Chain.m_Seq.end());
		}

	uint K = 10;
	if (optset_k)
		K = opt_k;

	const double MaxMDL = 0.38;
	const double MaxNENMed = 8.5;
	const uint MinL = 125;
	const uint MaxL = 200;
	const uint Iters = 1000000;

	uint NrLoMDL = 0;
	uint NrLoNENMed = 0;
	uint FailCount = 0;
	uint GoodIdx = 0;
	uint ThreadCount = GetRequestedThreadCount();
	time_t t0 = time(0);
	ProgressStep(GoodIdx, K, "Working %u good, %u failed, %u MDL, %u NEN",
					GoodIdx, FailCount, NrLoMDL, NrLoNENMed);
#pragma omp parallel num_threads(ThreadCount)
	{
	FakeChain FC;
	FC.m_Library = &Frags;
	int ThreadIndex = GetThreadIndex();
	for (;;)
		{
		Lock();
		time_t t = time(0);
		if (t > t0)
			{
			if (GoodIdx < K)
			ProgressStep(GoodIdx, K, "Working %u good, %u failed, %u MDL, %u NEN",
						 GoodIdx, FailCount, NrLoMDL, NrLoNENMed);
			t0 = t;
			}
		Unlock();
		FC.m_Chain.m_Label = "FC";
		uint L = MinL + randu32()%(MaxL - MinL);
		bool Ok = FC.MakeFake(L);
		if (!Ok)
			{
			++FailCount;
			continue;
			}
		FC.Validate();
		double MDL = FC.m_MDL;
		double NENMed = FC.m_NENMed;

		Lock();
		if (MDL > MaxMDL)
			++NrLoMDL;
		if (NENMed > MaxNENMed)
			++NrLoNENMed;
		Unlock();

		if (MDL > MaxMDL || NENMed > MaxNENMed)
			continue;

		Lock();
		++GoodIdx;
		Psa(FC.m_Chain.m_Label, "fake%u", GoodIdx);
		FC.ToTsv(fTsv);
		FC.m_Chain.ToCal(fCal);
		if (optset_outdir)
			{
			asserta(optset_loaddir);
			string FN = opt_outdir;
			if (!EndsWith(FN, "/"))
				FN += "/";
			if (optset_loaddir)
				{
				PDBChain Chain;
				FC.BuildPDB(opt_loaddir, Chain);
				Ps(FN, "%sfake%u.pdb", opt_outdir, GoodIdx);
				Chain.ToPDB(FN);
				}
			}
		Unlock();

		if (GoodIdx >= K)
			break;
		}
	}
	CloseStdioFile(fCal);
	CloseStdioFile(fTsv);
	}

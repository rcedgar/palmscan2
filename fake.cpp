#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "fakechain.h"
#include "omplock.h"
#include <time.h>

static bool IsAccept(double MDL, double NENMed)
	{
	if (MDL < 0.4 && NENMed < 8.5)
		return true;
	if (MDL < 0.5 && NENMed < 9)
		return randu32()%4 == 0;
	if (MDL < 0.6 && NENMed < 10)
		return randu32()%16 == 0;
	//if (MDL < 0.7 && NENMed < 10.5)
	//	return randu32()%32 == 0;
	return randu32()%32 == 0;
	}

void cmd_fake()
	{
	const string &InputFN = g_Arg1;
	FILE *fCal = CreateStdioFile(opt_cal);
	FILE *fTsv = CreateStdioFile(opt_tsv);

	vector<PDBChain *> Frags;
	ReadChains(InputFN, Frags);

	asserta(optset_sample_size);
	asserta(optset_minlen);
	asserta(optset_maxlen);

	const uint SampleSize = opt_sample_size;
	const uint MinL = opt_minlen;
	const uint MaxL = opt_maxlen;

	uint NrReject = 0;
	uint FailCount = 0;
	uint GoodIdx = 0;
	uint ThreadCount = GetRequestedThreadCount();
	time_t t0 = 0;
	time_t ta = 0;

#pragma omp parallel num_threads(ThreadCount)
		{
		FakeChain FC;
		FC.m_Library = &Frags;
		int ThreadIndex = GetThreadIndex();
		for (;;)
			{
			Lock();
			uint MyGoodIdx = GoodIdx;
			time_t t = time(0);
			if (t > t0)
				{
				if (MyGoodIdx+1 < SampleSize)
					ProgressStep(MyGoodIdx, SampleSize,
								 "Working %u good, %u failed, %u rejects",
								 GoodIdx, FailCount, NrReject);
				}
			Unlock();
			if (MyGoodIdx > SampleSize)
				break;

			FC.m_Chain.m_Label = "FC";
			uint L = MinL + randu32()%(MaxL - MinL);
			bool Ok = FC.MakeFake(L);
			FC.Validate();
			if (!Ok)
				{
				++FailCount;
				continue;
				}
			double MDL = FC.m_MDL;
			double NENMed = FC.m_NENMed;

			Lock();
			uint PctDone = GoodIdx*100/SampleSize;
			bool Accept = IsAccept(MDL, NENMed);
			if (Accept)
				ta = time(0);
			else
				++NrReject;
			Unlock();

			if (!Accept)
				continue;

			Lock();
			if (GoodIdx < SampleSize)
				{
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
						PDBChain Chain, ChainShuffledCA;
						FC.BuildPDB(opt_loaddir, Chain, ChainShuffledCA);
						Ps(FN, "%sfake%u.pdb", opt_outdir, GoodIdx);
						Chain.ToPDB(FN);

						if (optset_outdir_shuffledca)
							{
							Ps(FN, "%sfake%u.pdb", opt_outdir_shuffledca, GoodIdx);
							ChainShuffledCA.ToPDB(FN);
							}
						}
					}
				}
			Unlock();

			if (GoodIdx >= SampleSize)
				break;
			}
		}
	ProgressStep(SampleSize-1, SampleSize,
					"Working %u good, %u failed, %u rejects",
					GoodIdx, FailCount, NrReject);
	CloseStdioFile(fCal);
	CloseStdioFile(fTsv);
	}

#include "myutils.h"
#include "tma.h"
#include "outputfiles.h"

void TMAlignPair(TMA &T, const PDBChain &Q, const PDBChain &R);

void cmd_tm_all_vs_all()
	{
	const string &InputFileName = opt_tm_all_vs_all;
	vector<PDBChain *> Chains;
	ReadChains(InputFileName, Chains);
	const uint ChainCount = SIZE(Chains);

	uint ThreadCount = GetRequestedThreadCount();

	vector<TMA *> Ts;
	for (uint i = 0; i < ThreadCount; ++i)
		{
		TMA *T = new TMA;
		Ts.push_back(T);
		}

	uint PairCount = (ChainCount*(ChainCount - 1))/2;

	uint i = 1;
	uint j = 0;
	uint PairIndex = 0;
#pragma omp parallel num_threads(ThreadCount)
	for (;;)
		{
		uint MyPairIndex = UINT_MAX;
#pragma omp critical
		{
		if (PairIndex < PairCount)
			{
			MyPairIndex = PairIndex;
			ProgressStep(PairIndex, (uint) PairCount, "Aligning");
			++PairIndex;
			}
		}
		Log("thread %u, pair %u\n", GetThreadIndex(), MyPairIndex);
		if (MyPairIndex == UINT_MAX)
			break;

		uint Myi = UINT_MAX;
		uint Myj = UINT_MAX;
#pragma omp critical
		{
		Myi = i;
		Myj = j;
		++j;
		if (j == i)
			{
			++i;
			j = 0;
			}
		}

		uint ThreadIndex = GetThreadIndex();
		TMA &T = *Ts[ThreadIndex];
		TMAlignPair(T, *Chains[Myi], *Chains[Myj]);
		}
	}

#include "myutils.h"
#include "tma.h"
#include "outputfiles.h"

void TMAlignPair(TMA &T, const PDBChain &Q, const PDBChain &R);
void Shuffle(vector<unsigned> &v);

void cmd_tm_subsample()
	{
	const string &InputFileName = opt_tm_subsample;
	vector<PDBChain *> Chains;
	ReadChains(InputFileName, Chains);
	const uint ChainCount = SIZE(Chains);
	if (ChainCount < 1)
		{
		Warning("%u chains, no pairs", ChainCount);
		return;
		}

	uint PairCount = (ChainCount*(ChainCount - 1))/2;

	uint ThreadCount = GetRequestedThreadCount();
	uint SampleSize = opt_sample_size;
	if (SampleSize > PairCount)
		{
		Warning("sample_size %u set to nr. pairs %u",
		  SampleSize, PairCount);
		SampleSize = PairCount;
		}

	vector<TMA *> Ts;
	for (uint i = 0; i < ThreadCount; ++i)
		{
		TMA *T = new TMA;
		Ts.push_back(T);
		}

	vector<uint> Indexes;
	for (uint i = 0; i < PairCount; ++i)
		Indexes.push_back(i);
	Shuffle(Indexes);

	vector<pair<uint, uint> > Pairs;
	for (uint i = 0; i < ChainCount; ++i)
		for (uint j = 0; j < i; ++j)
			Pairs.push_back(pair<uint, uint>(i, j));

	uint Counter = 0;
#pragma omp parallel num_threads(ThreadCount)
		{
		for (;;)
			{
			uint MyIndex = UINT_MAX;
			uint MyCounter = UINT_MAX;
#pragma omp critical
			{
			MyCounter = Counter;
			if (Counter < SampleSize)
				{
				MyIndex = Indexes[MyCounter];
				ProgressStep(Counter, SampleSize, "Aligning");
				++Counter;
				}
			} // end critical

			asserta(MyCounter <= SampleSize);
			if (MyCounter == SampleSize)
				break;

			uint i = Pairs[MyIndex].first;
			uint j = Pairs[MyIndex].second;

			uint ThreadIndex = GetThreadIndex();
			TMA &T = *Ts[ThreadIndex];

			TMAlignPair(T, *Chains[i], *Chains[j]);
			}
		} // end parallel
	}

#include "myutils.h"
#include "pdbchain.h"
#include "trisearcher.h"
#include "abcxyz.h"

void InitTS(TriSearcher &TS);

void cmd_alignabc()
	{
	const string &QueryFN = opt_alignabc;
	const string &RefFN = opt_ref;

	TriSearcher TS;
	InitTS(TS);

	vector<PDBChain *> Qs;
	vector<PDBChain *> Rs;

	PDBChain::ReadChainsFromFile(RefFN, Rs);
	PDBChain::ReadChainsFromFile(QueryFN, Qs);

	const uint NQ = SIZE(Qs);
	const uint NR = SIZE(Rs);

	for (uint i = 0; i < NQ; ++i)
		{
		PDBChain &Q = *Qs[i];
//		Q.LogMe();
		for (uint j = 0; j < NR; ++j)
			{
			PDBChain &R = *Rs[j];
//			R.LogMe();

			TS.Search(Q, R);
			TS.SetTriForm();
			TS.LogMe();
			}
		}
	}

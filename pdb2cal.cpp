#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"

void cmd_pdb2cal()
	{
	const string &FN = opt_pdb2cal;

	vector<PDBChain *> Chains;
	ReadChains(FN, Chains);

	const uint N = SIZE(Chains);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Writing output");
		const PDBChain &Chain = *Chains[i];
		Chain.ToCal(g_fcal);
		}
	}

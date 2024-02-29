#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"

void cmd_reverse_chains()
	{
	const string &InputFN = opt_reverse_chains;

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);

	const uint N = SIZE(Chains);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Reversing");
		const PDBChain &Chain = *Chains[i];
		PDBChain RevChain;
		Chain.GetReverse(RevChain);
		RevChain.ToCal(g_fcal);
		if (optset_pdb)
			{
			RevChain.ToPDB(g_fpdb);
			break;
			}
		}
	}

#include "myutils.h"
#include "pdbchain.h"
#include "fakechain.h"

void cmd_cubes()
	{
	const string &InputFN = g_Arg1;

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);
	const uint N = SIZE(Chains);

	FakeChain FC;

	bool Ok = FC.TryAppendFrag(*Chains[0]);
	FC.LogMe();
	asserta(Ok);
	FC.Validate();

	for (uint k = 0; k < 8; ++k)
		{
		for (uint i = 1; i < N; ++i)
			{
			uint j = randu32()%N;
			Ok = FC.TryAppendFrag(*Chains[j]);
			if (Ok)
				{
				FC.Validate();
				Progress("Success %u\n", k);
				break;
				}
			}
		}
	if (!Ok)
		Progress("Fail\n");

	FC.m_Chain.ToPDB("cubes_out.pdb");
	}
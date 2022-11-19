#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"

// Input is CAL which is palmprint-only with coordinates,
//  output is rotated into PPC.
void cmd_rotate()
	{
	const string &InputFN = opt_rotate;

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains, false);

	const uint N = SIZE(Chains);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Rotating");
		const PDBChain &Chain = *Chains[i];
		Chain.CheckMotifCoords();
		uint PosA = Chain.m_MotifPosVec[0];
		uint PosB = Chain.m_MotifPosVec[1];
		uint PosC = Chain.m_MotifPosVec[2];

		PDBChain PPC;
		Chain.GetPPC(PosA, PosB, PosC, PPC);
		asserta(PPC.CheckPPCMotifCoords());
		PPC.ToCal(g_fppc);
		}
	}

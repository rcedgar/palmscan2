#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "fakechain.h"

void cmd_fake()
	{
	const string &InputFN = g_Arg1;
	FILE *fCal = CreateStdioFile(opt_cal);

	FakeChain FC;
	ReadChains(InputFN, FC.m_Library);

	uint K = 10;
	if (optset_k)
		K = opt_k;

	uint FailCount = 0;
	for (uint k = 0; k < K; ++k)
		{
		ProgressStep(k, K, "Working %u failed", FailCount);
		FC.m_Chain.m_Label = "FC";
		bool Ok = FC.MakeFake(150);
		if (!Ok)
			{
			++FailCount;
			continue;
			}
		//FC.LogMe();
		if (optset_outdir)
			{
			string FN = opt_outdir;
			if (!EndsWith(FN, "/"))
				FN += "/";
			Psa(FC.m_Chain.m_Label, "fake%u", k);
			Psa(FN, "fake%u.pdb", k);
			FC.m_Chain.ToPDB(FN);
			FC.m_Chain.ToCal(fCal);
			}
		}
	CloseStdioFile(fCal);
	}

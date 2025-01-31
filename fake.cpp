#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "fakechain.h"

void cmd_fake()
	{
	const string &InputFN = g_Arg1;
	FILE *fCal = CreateStdioFile(opt_cal);
	FILE *fTsv = CreateStdioFile(opt_tsv);

	vector<PDBChain *> Frags;
	ReadChains(InputFN, Frags);

	FakeChain FC;
	FC.m_Library = &Frags;

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
		FC.Validate();
		FC.ToTsv(fTsv);

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

		if (optset_loaddir)
			{
			FC.LogMe();
			FC.m_Chain.ToPDB("fake_caonly.pdb");
			PDBChain Chain;
			FC.BuildPDB(opt_loaddir, Chain);
			Chain.ToPDB("fake.pdb");
			Die("TODO");
			}
		}
	CloseStdioFile(fCal);
	CloseStdioFile(fTsv);
	}

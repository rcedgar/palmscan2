#include "myutils.h"
#include "chainfaker.h"

void cmd_ffake()
	{
	asserta(optset_sample_size);
	const uint SampleSize = opt_sample_size;

	ChainFaker CF;
	CF.m_Trace = true;
	CF.ReadSCOP40(g_Arg1);

	const uint ChainCount = SIZE(CF.m_SCOP40);

	string PDBOutDir = opt_outdir;
	if (!EndsWith(PDBOutDir, "/"))
		PDBOutDir += "/";

	FILE *fCal = CreateStdioFile(opt_cal);

	uint n = 0;
	const uint TRIES = 3*SampleSize;
	for (uint Try = 0; Try < TRIES; ++Try)
		{
		ProgressStep(n, SampleSize, "Faking");

		uint ChainIdx = randu32()%ChainCount;
		PDBChain FakeChain;
		bool Ok = CF.MakeFake(ChainIdx, FakeChain);
		if (!Ok)
			continue;
		++n;
		string Fold;
		CF.GetFoldStr(Fold);
		string Label;
		Ps(Label, "Fake_%u_%s", n+1, Fold.c_str());
		assert(Fold[1] == '.');
		Fold[1] = '_';
		CF.m_FakeChain->m_Label = Label;

		CF.m_FakeChain->ToCal(fCal);
		if (optset_outdir)
			CF.ToPDB(PDBOutDir + Label + ".pdb");
		if (n >= SampleSize)
			break;
		}
	if (n < SampleSize)
		Warning("Generated %u / %u", n, SampleSize);

	CloseStdioFile(fCal);
	}

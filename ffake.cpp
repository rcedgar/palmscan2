#include "myutils.h"
#include "chainfaker.h"

void cmd_ffake()
	{
	asserta(optset_sample_size);
	const uint SampleSize = opt_sample_size;

	ChainFaker CF;
	CF.m_Trace = true;

	asserta(optset_minlen);
	asserta(optset_maxlen);
	asserta(optset_minreppct);
	asserta(optset_shufflepct);
	CF.m_MinTakeoutLen = opt_minlen;
	CF.m_MaxTakeoutLen = opt_maxlen;
	CF.m_MinReplacedPct = opt_minreppct;
	CF.m_ShufflePct = opt_shufflepct;

	ProgressLog("N=%u rep [%u-%u] %u%%, shuffle %u%%\n",
				SampleSize,
				CF.m_MinTakeoutLen,
				CF.m_MaxTakeoutLen,
				CF.m_MinReplacedPct,
				CF.m_ShufflePct);

	CF.ReadSCOP40(g_Arg1);

	const uint ChainCount = SIZE(CF.m_SCOP40);

	string PDBOutDir = opt_outdir;
	if (!EndsWith(PDBOutDir, "/"))
		PDBOutDir += "/";

	FILE *fCal = CreateStdioFile(opt_cal);
	FILE *fTsv = CreateStdioFile(opt_tsv);

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
		string SF;
		CF.GetSFStr(SF);
		string Label;
		Ps(Label, "Fake_%u_%s", n, SF.c_str());
		CF.m_FakeChain->m_Label = Label;

		CF.m_FakeChain->ToCal(fCal);
		CF.ToTsv(fTsv);
		if (optset_outdir)
			CF.ToPDB(PDBOutDir + Label + ".pdb");
		if (n >= SampleSize)
			break;
		}
	if (n < SampleSize)
		Warning("Generated %u / %u", n, SampleSize);

	CloseStdioFile(fCal);
	CloseStdioFile(fTsv);
	}

#include "myutils.h"
#include "filler.h"

void cmd_filler_train()
	{
	const string &FN = opt_filler_train;

	vector<PDBChain *> Chains;
	ReadChains(FN, Chains);
	const uint ChainCount = SIZE(Chains);

	Filler F;
	F.Train(Chains);
	FILE *fMx = CreateStdioFile("filler.mx");
	F.ScoreMxToFile(fMx);
	CloseStdioFile(fMx);

	vector<string> Lines;
	ReadLinesFromFile("filler.mx", Lines);

	Filler F2;
	F2.FromLines(Lines);
	}

#include "myutils.h"
#include "scop40bit.h"
#include "sort.h"

void cmd_scop40bit_fam()
	{
	asserta(optset_output);
	const string &FN = g_Arg1;
	SCOP40Bit SB;
	SB.ReadDomInfo();
	SB.ReadHits_Bin(FN);
	}
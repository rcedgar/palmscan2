#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "outputfiles.h"

void cmd_pdb2cal()
	{
	const string &FN = opt_pdb2cal;

	ChainReader CR;
	CR.Open(FN);

	PDBChain Chain;
	uint Count = 0;
	for (;;)
		{
		bool Ok = CR.GetNext(Chain);
		if (!Ok)
			break;

		if (++Count%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u converted >%s\r",
			  sPct.c_str(), Count, Chain.m_Label.c_str());
			}
		Chain.ToCal(g_fcal);
		}
	Progress("100.0%% done, %u converted\n", Count);
	}

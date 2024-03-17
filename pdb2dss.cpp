#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "outputfiles.h"
#include "dss.h"

void cmd_pdb2dss()
	{
	const string &FN = opt_pdb2dss;
	asserta(optset_output);
	FILE *fout = CreateStdioFile(opt_output);

	ChainReader CR;
	CR.Open(FN);

	DSS D;
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
		D.Init(Chain);
		string s;
		D.ToStr(s);
		SeqToFasta(fout, Chain.m_Label, s); 
		}
	Progress("100.0%% done, %u converted\n", Count);
	}

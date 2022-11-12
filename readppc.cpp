#include "myutils.h"
#include "pdbchain.h"
#include "calreader.h"

void ReadPpc(const string &FN, vector<PDBChain *> &Chains)
	{
	Chains.clear();
	CalReader CR;
	CR.Open(FN);
	uint n = 0;
	for (;;)
		{
		++n;
		if (n%100 == 0)
			{
			string sPct;
			CR.GetPctDone(sPct);
			Progress("Reading %s (%s%%)\r", FN.c_str(), sPct.c_str());
			}
		PDBChain *Chain = new PDBChain;
		bool Ok = CR.GetNext(*Chain);
		if (!Ok)
			return;
		Chains.push_back(Chain);
		}
	Progress("Reading %s (100.0%%)\n");
	}

#include "myutils.h"
#include "pdbchain.h"

void ReadLinesFromFile(const string &FileName, vector<string> &Lines);

void ReadPDBs(const vector<string> &FileNames, vector<PDBChain *> &Structures)
	{
	Structures.clear();
	const uint N = SIZE(FileNames);
	for (uint i = 0; i < N; ++i)
		{
		const string &FileName = FileNames[i];
		PDBChain *Structure = new PDBChain;
		vector<string> Lines;
		ReadLinesFromFile(FileName, Lines);
		vector<PDBChain *> Chains;
		PDBChain::ReadChainsFromLines(FileName, Lines, Chains);
		const uint NC = SIZE(Chains);
		for (uint j = 0; j < NC; ++j)
			Structures.push_back(Chains[j]);
		}
	}

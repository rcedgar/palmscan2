#include "myutils.h"
#include "pdbchain.h"

void GetFileNames(const string &SpecFileName, vector<string> &FileNames);
void ReadLinesFromFile(const string &FileName, vector<string> &Lines);

void ReadPDBs(const vector<string> &FileNames, vector<PDBChain *> &Structures)
	{
	Structures.clear();
	const uint N = SIZE(FileNames);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Reading PDBs");
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

void ReadPDBs(const string &FileName, vector<PDBChain *> &Structures)
	{
	vector<string> FileNames;
	GetFileNames(FileName, FileNames);
	ReadPDBs(FileNames, Structures);
	}

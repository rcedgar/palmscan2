#include "myutils.h"
#include "pdb.h"

void ReadPDBs(const vector<string> &FileNames, vector<PDB *> &Structures)
	{
	Structures.clear();
	const uint N = SIZE(FileNames);
	for (uint i = 0; i < N; ++i)
		{
		const string &FileName = FileNames[i];
		PDB *Structure = new PDB;
		Structure->FromFile(FileName);
		Structures.push_back(Structure);
		}
	}

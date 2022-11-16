#include "myutils.h"
#include "pdbchain.h"
#include "calreader.h"

void ReadLinesFromFile(const string &FileName, vector<string> &Lines);

void GetLabelFromFileName(const string &FileName, string &Label)
	{
	vector<string> Fields;
	Split(FileName, Fields, '/');
	uint FieldCount = SIZE(Fields);
	asserta(FieldCount > 0);
	Label = Fields[FieldCount-1];
	if (EndsWith(Label, ".pdb"))
		{
		uint n = SIZE(Label);
		if (n > 4)
			Label = Label.substr(0, n-4);
		}

	if (Label.size() == 6 && Label[4] == '_')
		Label = Label.substr(0, 4);
	}

void ReadChains(const vector<string> &FileNames,
  vector<PDBChain *> &Structures, bool SaveAtoms)
	{
	Structures.clear();
	const uint N = SIZE(FileNames);
	for (uint i = 0; i < N; ++i)
		{
		const string &FileName = FileNames[i];
		string Label;
		GetLabelFromFileName(FileName, Label);

		vector<string> Lines;
		ReadLinesFromFile(FileName, Lines);

		ProgressStep(i, N, "Reading chains %s", Label.c_str());

		vector<PDBChain *> Chains;
		PDBChain::ChainsFromLines(Label, Lines, Chains, SaveAtoms);
		const uint NC = SIZE(Chains);
		for (uint j = 0; j < NC; ++j)
			Structures.push_back(Chains[j]);
		}
	}

void ReadChainsCal(const string &FileName, vector<PDBChain *> &Structures)
	{
	CalReader CR;
	CR.Open(FileName);
	uint N = 0;
	for (;;)
		{
		PDBChain &Q = *new PDBChain;
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;
		Structures.push_back(&Q);
		++N;
		if (N%100 == 0)
			{
			string sPct;
			CR.GetPctDone(sPct);
			Progress("Reading %s %s%% done\r",
			  FileName.c_str(), sPct.c_str());
			}
		}
	Progress("Reading %s 100.0%% done\n", FileName.c_str());
	}

void ReadChains(const string &FileName, vector<PDBChain *> &Structures,
  bool SaveAtoms)
	{
	if (FileName.empty())
		Die("Missing chains filename");

	if (EndsWith(FileName, ".cal") || EndsWith(FileName, ".ppc"))
		{
		ReadChainsCal(FileName, Structures);
		return;
		}
	else if (EndsWith(FileName, ".files"))
		{
		vector<string> FileNames;
		ReadLinesFromFile(FileName, FileNames);
		ReadChains(FileNames, Structures, SaveAtoms);
		}

	string Label;
	GetLabelFromFileName(FileName, Label);

	vector<string> Lines;
	ReadLinesFromFile(FileName, Lines);

	vector<PDBChain *> Chains;
	PDBChain::ChainsFromLines(Label, Lines, Chains, SaveAtoms);
	const uint NC = SIZE(Chains);
	for (uint j = 0; j < NC; ++j)
		Structures.push_back(Chains[j]);
	}

#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include <set>

void ReadLinesFromFile(const string &FileName, vector<string> &Lines);
void ReadChains(const string &FileName, vector<PDBChain *> &Structures);

void cmd_getchains()
	{
	const string &InputFileName = opt_getchains;

	asserta(opt_labels != "");

	vector<string> Labels;
	ReadLinesFromFile(opt_labels, Labels);
	set<string> LabelSet;
	const uint LabelCount = SIZE(Labels);
	for (uint i = 0; i < LabelCount; ++i)
		{
		const string &Label = Labels[i];
		LabelSet.insert(Label);
		}

	vector<PDBChain *> Chains;
	ReadChains(InputFileName, Chains);
	const uint N = SIZE(Chains);

	uint n = 0;
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N+1, "Searching %u / %u found", n, LabelCount);

		PDBChain &Q = *Chains[i];
		const string &Label = Q.m_Label;
		if (LabelSet.find(Label) != LabelSet.end())
			{
			Q.ToCal(g_fcal);
			++n;
			}
		}
	ProgressStep(N, N+1, "Searching %u / %u found", n, LabelCount);
	}

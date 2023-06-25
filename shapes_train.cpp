#include "myutils.h"
#include "shapesearcher.h"
#include "rdrpsearcher.h"
#include <map>

/***
shapes_train
	Input: (a) tsv file with shape sequences or (b) PSSMs
	Output: calculate and save contact map profile for
	  all-vs-all shapes.

	Tests ShapeSearcher ABC detection against training set.
***/

void cmd_shapes_train()
	{
	const string &InputFileName1 = opt_shapes_train;
	asserta(optset_input);
	const string &InputFileName2 = opt_input;
	
	vector<string> ChainLabels;
	vector<string> MotifNames;
	vector<uint> MotifLengths;
	vector<vector<string> > MotifSeqsVec;
	vector<PDBChain *> Chains;
	ReadChains(InputFileName1, Chains);
	const uint ChainCount = SIZE(Chains);

	GetTrainingMotifs(InputFileName2, Chains,
	  ChainLabels, MotifNames, MotifLengths, MotifSeqsVec, opt_extend_abc);
	const uint N = SIZE(ChainLabels);
	asserta(SIZE(MotifSeqsVec) == N);

	map<string, uint> ChainLabelToIndex;
	for (uint i = 0; i < N; ++i)
		{
		const string &ChainLabel = ChainLabels[i];
		asserta(ChainLabelToIndex.find(ChainLabel) == ChainLabelToIndex.end());
		ChainLabelToIndex[ChainLabel] = i;
		}

	vector<PDBChain *> Chains2;
	vector<vector<string> > MotifSeqsVec2;
	uint LabelNotFoundCount = 0;
	for (uint i = 0; i < ChainCount; ++i)
		{
		PDBChain &Chain = *Chains[i];
		const string &Label = Chain.m_Label;
		map<string, uint>::const_iterator p = ChainLabelToIndex.find(Label);
		if (p == ChainLabelToIndex.end())
			{
			++LabelNotFoundCount;
			continue;
			}
		uint Index = p->second;
		Chains2.push_back(&Chain);
		MotifSeqsVec2.push_back(MotifSeqsVec[Index]);
		}

	if (LabelNotFoundCount > 0)
		Warning("%u chain labels not found", LabelNotFoundCount);

	Shapes S;
	S.Init(MotifNames, MotifLengths);
	S.Train(Chains2, MotifSeqsVec2, opt_extend_abc);
	S.ToFile(opt_output);

	ShapeSearcher::TestABC(S, Chains2, MotifSeqsVec2);
	}

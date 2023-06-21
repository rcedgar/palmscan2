#include "myutils.h"
#include "shapes.h"
#include <map>

static void ReadMotifsFile(const string &FileName,
  vector<string> &ChainLabels,
  vector<string> &MotifNames, vector<uint> &MotifLengths,
  vector<vector<string> > &MotifSeqsVec)
	{
	ChainLabels.clear();
	MotifNames.clear();
	MotifLengths.clear();
	MotifSeqsVec.clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	vector<string> Fields;
	Split(Line, Fields, '\t');
	asserta(Fields[0] == "Label");
	asserta(SIZE(Fields) > 1);
	const uint MotifCount = SIZE(Fields) - 1;
	for (uint i = 0; i < MotifCount; ++i)
		{
		MotifLengths.push_back(UINT_MAX);
		MotifNames.push_back(Fields[i+1]);
		}

	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == MotifCount + 1);

		const string &ChainLabel = Fields[0];
		ChainLabels.push_back(ChainLabel);

		vector<string> MotifSeqs;
		for (uint i = 0; i < MotifCount; ++i)
			MotifSeqs.push_back(Fields[i+1]);
		MotifSeqsVec.push_back(MotifSeqs);

		for (uint i = 0; i < MotifCount; ++i)
			{
			const string &MotifSeq = MotifSeqs[i];
			if (MotifSeq == "" || MotifSeq == ".")
				continue;
			uint L = SIZE(MotifSeq);
			if (MotifLengths[i] == UINT_MAX)
				MotifLengths[i] = L;
			else
				asserta(L == MotifLengths[i]);
			}
		}
	CloseStdioFile(f);
	}

void cmd_shapes_train()
	{
	const string &InputFileName1 = opt_shapes_train;
	asserta(optset_input);
	const string &InputFileName2 = opt_input;
	
	vector<string> ChainLabels;
	vector<string> MotifNames;
	vector<uint> MotifLengths;
	vector<vector<string> > MotifSeqsVec;
	ReadMotifsFile(InputFileName2,
	  ChainLabels, MotifNames, MotifLengths, MotifSeqsVec);
	const uint N = SIZE(ChainLabels);
	asserta(SIZE(MotifSeqsVec) == N);

	map<string, uint> ChainLabelToIndex;
	for (uint i = 0; i < N; ++i)
		{
		const string &ChainLabel = ChainLabels[i];
		asserta(ChainLabelToIndex.find(ChainLabel) == ChainLabelToIndex.end());
		ChainLabelToIndex[ChainLabel] = i;
		}

	vector<PDBChain *> Chains;
	ReadChains(InputFileName1, Chains);
	const uint ChainCount = SIZE(Chains);

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
	S.Train(Chains2, MotifSeqsVec2);
	ProgressLog("Done.\n");
	}

#include "myutils.h"
#include "shapesearcher.h"
#include "rdrpsearcher.h"

static void GetTrainingMotifs_PSSMs(const vector<PDBChain *> &Chains,
  vector<string> &ChainLabels, vector<string> &MotifNames,
  vector<uint> &MotifLengths, vector<vector<string> > &MotifSeqsVec,
  bool ExtendABC)
	{
	const double MIN_RS_SCORE = 10;

	ChainLabels.clear();
	MotifNames.clear();
	MotifLengths.clear();
	MotifSeqsVec.clear();

	MotifNames.push_back("A");
	MotifNames.push_back("B");
	MotifNames.push_back("C");

// A -2 +2			12 => 16
// B 0 8			14 => 22
// C -2 +2			8 => 12
	if (ExtendABC)
		{
		MotifLengths.push_back(12+2+2);
		MotifLengths.push_back(14+0+8);
		MotifLengths.push_back(8+2+2);
		}
	else
		{
		MotifLengths.push_back(12);
		MotifLengths.push_back(14);
		MotifLengths.push_back(8);
		}

	RdRpModel Model;
	GetRdrpModel(Model);

	RdRpSearcher RS;
	RS.Init(Model);

	const uint N = SIZE(Chains);
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Chain = *Chains[i];
		const string &Label = Chain.m_Label;
		const string &Seq = Chain.m_Seq;
		RS.Search(Label, Seq);
		if (RS.m_TopPalmHit.m_Score < MIN_RS_SCORE)
			continue;

		uint GroupIndex = RS.m_TopPalmHit.m_GroupIndex;
		bool Permuted = RS.m_TopPalmHit.m_Permuted;

		string GroupName;
		RS.GetGroupName(GroupIndex, GroupName);

		uint PosA = RS.GetMotifPos(0);
		uint PosB = RS.GetMotifPos(1);
		uint PosC = RS.GetMotifPos(2);
		if (PosA == UINT_MAX || PosB == UINT_MAX || PosC == UINT_MAX)
			continue;
		if (Permuted)
			continue;
		if (GroupName == "Birna")
			{
			Warning("Birna");
			continue;
			}

		vector<string> MotifSeqs(3);
		if (ExtendABC)
			{
			Chain.GetSubSeq(PosA-2, 12+2+2, MotifSeqs[0]);
			Chain.GetSubSeq(PosB, 14+0+8, MotifSeqs[1]);
			Chain.GetSubSeq(PosC-2, 8+2+2, MotifSeqs[2]);
			}
		else
			{
			Chain.GetSubSeq(PosA, 12, MotifSeqs[0]);
			Chain.GetSubSeq(PosB, 14, MotifSeqs[1]);
			Chain.GetSubSeq(PosC, 8, MotifSeqs[2]);
			}

		ChainLabels.push_back(Label);
		MotifSeqsVec.push_back(MotifSeqs);
		}
	}

// If FileName is tsv file with motifs or "PSSMs"
void GetTrainingMotifs(const string &FileName,
  const vector<PDBChain *> &Chains, vector<string> &ChainLabels,
  vector<string> &MotifNames, vector<uint> &MotifLengths,
  vector<vector<string> > &MotifSeqsVec, bool ExtendABC)
	{
	if (FileName == "PSSMs")
		{
		GetTrainingMotifs_PSSMs(Chains, ChainLabels,
		  MotifNames, MotifLengths, MotifSeqsVec, ExtendABC);
		return;
		}

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

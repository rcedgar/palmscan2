#include "myutils.h"
#include "shapesearcher.h"
#include "rdrpsearcher.h"
#include "outputfiles.h"
#include <map>

/***
make_trunc_bench
	Inputs:
		structures (PDB, cal, .files)
		tsv file with ABC sequences

	Output:
		cal truncated around pp boundaries 
		  to check for FP motifs due to edge effects
***/

static void Do1(const PDBChain &Chain, const string &A, const string &C, uint CL);

void cmd_make_trunc_bench()
	{
	const string &InputFileName1 = opt_make_trunc_bench;
	asserta(optset_input);
	const string &InputFileName2 = opt_input;
	
	vector<string> ChainLabels;
	vector<string> MotifNames;
	vector<uint> MotifLengths;
	vector<vector<string> > MotifSeqsVec;
	vector<PDBChain *> Chains;
	Progress("Reading chains\n");
	ReadChains(InputFileName1, Chains);
	const uint ChainCount = SIZE(Chains);

	GetTrainingMotifs(InputFileName2, Chains,
	  ChainLabels, MotifNames, MotifLengths, MotifSeqsVec);
	asserta(SIZE(MotifNames) == 3);
	asserta(MotifNames[0] == "A");
	asserta(MotifNames[1] == "B");
	asserta(MotifNames[2] == "C");
	uint CL = MotifLengths[2];

	const uint N = SIZE(ChainLabels);
	asserta(SIZE(MotifSeqsVec) == N);

	map<string, uint> ChainLabelToIndex;
	for (uint i = 0; i < N; ++i)
		{
		const string &ChainLabel = ChainLabels[i];
		asserta(ChainLabelToIndex.find(ChainLabel) == ChainLabelToIndex.end());
		ChainLabelToIndex[ChainLabel] = i;
		}

	for (uint i = 0; i < ChainCount; ++i)
		{
		PDBChain &Chain = *Chains[i];
		const string &Label = Chain.m_Label;
		map<string, uint>::const_iterator p = ChainLabelToIndex.find(Label);
		if (p == ChainLabelToIndex.end())
			continue;
		uint Index = p->second;
		const vector<string> &MotifSeqs = MotifSeqsVec[Index];
		Do1(Chain, MotifSeqs[0], MotifSeqs[2], CL);
		}
	}

static void Do1(const PDBChain &Chain, const string &A, const string &C, uint CL)
	{
	if (A == "" || C == "")
		return;
	if (A == "." || C == ".")
		return;

	const string &Seq = Chain.m_Seq;
	size_t stPosA = Seq.find(A);
	size_t stPosC = Seq.find(C);
	if (stPosA == string::npos || stPosC == string::npos)
		return;
	if (stPosC < stPosA)
		return;
	const char *Label = Chain.m_Label.c_str();

	uint L = SIZE(Seq);
	uint PosA = uint(stPosA);
	uint PosC = uint(stPosC) + CL;
	PDBChain TruncChain;
	Chain.GetRange(PosA, PosC, TruncChain);
	string Annot;
	Ps(Annot, ".00 A:%s C:%s", A.c_str(), C.c_str());
	TruncChain.m_Label += Annot;
	TruncChain.ToCal(g_fcal);
	return;//@@@@@@@@@@@
	for (uint i = 0; i < 4; ++i)
		{
		uint Left = randu32()%32 + 4;
		uint Right = randu32()%32 + 4;

		int LoShrink = int(PosA) + int(Left);
		int HiShrink = int(PosC) - int(Right);

		int LoGrow = int(PosA) - int(Left);
		int HiGrow = int(PosC) + int(Right);
		if (LoGrow < 0)
			LoGrow = 0;
		if (HiGrow >= int(L))
			HiGrow = int(L) - 1;

		if (HiShrink - LoShrink > 40)
			{
			Chain.GetRange(uint(LoShrink), uint(HiShrink), TruncChain);
			TruncChain.m_Label += ".ss";
			TruncChain.ToCal(g_fcal);
			}

		if (HiGrow - LoShrink > 40)
			{
			Chain.GetRange(uint(LoShrink), uint(HiGrow), TruncChain);
			TruncChain.m_Label += ".sg";
			TruncChain.ToCal(g_fcal);
			}

		if (HiShrink - LoGrow > 40)
			{
			Chain.GetRange(uint(LoGrow), uint(HiShrink), TruncChain);
			TruncChain.m_Label += ".gs";
			TruncChain.ToCal(g_fcal);
			}

		if (HiGrow - LoGrow > 40)
			{
			Chain.GetRange(uint(LoGrow), uint(HiGrow), TruncChain);
			TruncChain.m_Label += ".gg";
			TruncChain.ToCal(g_fcal);
			}
		}
	}

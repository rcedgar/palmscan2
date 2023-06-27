#include "myutils.h"
#include "shapesearcher.h"
#include "quarts.h"
#include <map>

/***
shapes_train
	Input: tsv file with shape sequences
	Input: contact map profile for all-vs-all shapes.
	Output: ShapeSearcher scores for training shapes
***/

static void ShapesTrainScore(ShapeSearcher &SS,
  const PDBChain &Chain, const vector<string> &ShapeSeqs,
  vector<vector<double> > &ScoresVec)
	{
	uint ShapeCount = SS.GetShapeCount();
	asserta(SIZE(ScoresVec) == ShapeCount);
	asserta(SIZE(ShapeSeqs) == ShapeCount);
	const string &ChainSeq = Chain.m_Seq;
	vector<uint> ShapePosVec;
	for (uint ShapeIndex = 0; ShapeIndex < ShapeCount; ++ShapeIndex)
		{
		const char *ShapeName = SS.GetShapeName(ShapeIndex);
		const string &ShapeSeq = ShapeSeqs[ShapeIndex];
		uint Pos = UINT_MAX;
		if (ShapeSeq != "" && ShapeSeq != ".")
			{
			size_t n = ChainSeq.find(ShapeSeq);
			asserta(n != string::npos);
			Pos = uint(n);
			}
		ShapePosVec.push_back(Pos);
		}

	SS.SetQuery(Chain);
	uint QL = Chain.GetSeqLength();
#if 0
	Log("\n");
	Log(">%s\n", Chain.m_Label.c_str());
#endif
	for (uint ShapeIndex = 0; ShapeIndex < ShapeCount; ++ShapeIndex)
		{
		uint Pos = ShapePosVec[ShapeIndex];
		if (Pos == UINT_MAX)
			continue;
		double Score = SS.GetSelfScore(ShapeIndex, Pos);
		ScoresVec[ShapeIndex].push_back(Score);
#if 0
		vector<uint> HitPosVec;
		const char *ShapeName = SS.GetShapeName(ShapeIndex);
		const string &ShapeSeq = ShapeSeqs[ShapeIndex];
		uint L = SS.GetShapeLength(ShapeIndex);
		vector<double> HitScores;
		char Letter = 0;
		uint LetterOffset = UINT_MAX;
		uint Lo = 0;
		uint Hi = QL - L;
		//Lo = Pos;
		//Hi = Pos;
		double MinScore = 0.6;
		SS.SearchOneShapeSelf(ShapeIndex, MinScore, Lo, Hi,
		  Letter, LetterOffset, HitPosVec, HitScores);
		Log("%3.3s", ShapeName);
		Log("  %5u", Pos);
		Log("  %16.16s", ShapeSeq.c_str());
		Log("  %.4g", Score);
		Log("\n");
#endif
		}
	}

void cmd_shapes_train_score()
	{
	const string &InputFileName1 = opt_shapes_train_score;
	asserta(optset_input);
	const string &InputFileName2 = opt_input;
	asserta(optset_shapes);
	const string &ShapesFileName = opt_shapes;
	
	vector<string> ChainLabels;
	vector<string> MotifNames;
	vector<uint> MotifLengths;
	vector<vector<string> > MotifSeqsVec;
	vector<PDBChain *> Chains;
	ReadChains(InputFileName1, Chains);
	const uint ChainCount = SIZE(Chains);

	GetTrainingMotifs(InputFileName2, Chains, ChainLabels,
	  MotifNames, MotifLengths, MotifSeqsVec, opt_extend_abc);
	const uint N = SIZE(ChainLabels);
	asserta(SIZE(MotifSeqsVec) == N);

	map<string, uint> ChainLabelToIndex;
	for (uint i = 0; i < N; ++i)
		{
		const string &ChainLabel = ChainLabels[i];
		asserta(ChainLabelToIndex.find(ChainLabel) == ChainLabelToIndex.end());
		ChainLabelToIndex[ChainLabel] = i;
		}


	Shapes S;
	S.FromFile(ShapesFileName);

	ShapeSearcher SS;
	SS.Init(S);

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

	const uint ChainCount2 = SIZE(Chains2);
	const uint ShapeCount = SS.GetShapeCount();
	vector<vector<double> > ScoresVec(ShapeCount);
	for (uint i = 0; i < ChainCount2; ++i)
		{
		ProgressStep(i, ChainCount2, "Training");

		const PDBChain &Chain = *Chains2[i];
		const vector<string> &ShapeSeqs = MotifSeqsVec2[i];
		ShapesTrainScore(SS, Chain, ShapeSeqs, ScoresVec);
		}

	for (uint ShapeIndex = 0; ShapeIndex < ShapeCount; ++ShapeIndex)
		{
		const vector<double> &Scores = ScoresVec[ShapeIndex];
		QuartsDouble Q;
		GetQuartsDouble(Scores, Q);
		const char *ShapeName = SS.GetShapeName(ShapeIndex);

		Log("\n");
		Log("%u", ShapeIndex);
		Log("  %3.3s: ", ShapeName);
		Q.LogMe();
		}
	}

#include "myutils.h"
#include "shapesearcher.h"
#include "outputfiles.h"
#include "motifsettings.h"

/***
shapes_train1
	Add one motif to ABC

	Inputs: 
		tsv1 file with label+ABC sequences
		tsv2 file with label+X sequences for one new motif X

	Output:
		Shapes for ABCX.
***/

static const double MINSCOREX = 0.4;

void ShapeSearcher::LogTrainStats(double MinScoreX,
  const vector<PDBChain *> &Chains,
  vector<string> &X1s,
  vector<string> &X2s,
  vector<double> &Score1s,
  vector<double> &Score2s)
	{
	const uint ChainCount = SIZE(Chains);
	asserta(SIZE(X1s) == ChainCount);
	asserta(SIZE(X2s) == ChainCount);
	asserta(SIZE(Score1s) == ChainCount);
	asserta(SIZE(Score2s) == ChainCount);

	uint SameSeq = 0;
	uint ScoreHigher = 0;
	uint ScoreLower = 0;
	uint Found = 0;
	uint Lost = 0;
	uint Neither = 0;
	uint Both = 0;
	uint ShiftGt0 = 0;
	uint ShiftGt4 = 0;
	uint Hits = 0;
	for (uint i = 0; i < ChainCount; ++i)
		{
		const string &Seq = Chains[i]->m_Seq;
		const string &X1 = X1s[i];
		const string &X2 = X2s[i];
		if (X1 == X2)
			++SameSeq;
		double Score1 = Score1s[i];
		double Score2 = Score2s[i];
		size_t Pos1 = Seq.find(X1);
		size_t Pos2 = Seq.find(X2);
		if (Score2 >= MinScoreX)
			++Hits;
		if (Score1 < MinScoreX && Score2 < MinScoreX)
			++Neither;
		else if (Score1 < MinScoreX && Score2 >= MinScoreX)
			++Found;
		else if (Score1 >= MinScoreX && Score2 < MinScoreX)
			++Lost;
		else if (Score1 >= MinScoreX && Score2 >= MinScoreX)
			++Both;
		else
			asserta(false);

		double Diff = fabs(Score2 - Score1);
		if (Diff > 0.05)
			{
			if (Score2 > Score1)
				++ScoreHigher;
			else
				++ScoreLower;
			}

		if (Pos1 != string::npos && Pos2 != string::npos)
			{
			int Shift = abs(int(Pos1) - int(Pos2));
			if (Shift > 0 && Shift < 4)
				++ShiftGt0;
			else if (Shift > 4)
				++ShiftGt4;
			}
		}
	ProgressLog("Min score %.4f\n", MinScoreX);
#define X(x)	ProgressLog("%16.16s  %u (%.1f%%)\n", #x, x, GetPct(x, ChainCount));
	X(ChainCount);
	X(Hits);
	X(SameSeq);
	X(ScoreHigher);
	X(ScoreLower);
	X(Found);
	X(Lost);
	X(Neither);
	X(Both);
	X(ShiftGt0);
	X(ShiftGt4);
#undef X
	}

void cmd_shapes_train1()
	{
	const string &ChainsFN = opt_shapes_train1;

	asserta(optset_input);
	asserta(optset_input1);
	const string &TsvFN_ABC = opt_input;
	const string &TsvFN_X = opt_input1;

	vector<PDBChain *> Chains;
	Progress("Reading chains\n");
	ReadChains(ChainsFN, Chains);
	const uint ChainCount = SIZE(Chains);
	vector<string> ChainLabels;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;
		ChainLabels.push_back(Label);
		}

	map<string, vector<string> > LabelToABC;
	map<string, string > LabelToX;

	FILE *f = OpenStdioFile(TsvFN_ABC);
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(Fields[0] == "Label");
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 4);
		vector<string> ABC;
		const string &Label = Fields[0];
		ABC.push_back(Fields[1]);
		ABC.push_back(Fields[2]);
		ABC.push_back(Fields[3]);
		LabelToABC[Label] = ABC;
		CheckABC(ABC);
		}
	CloseStdioFile(f);

	f = OpenStdioFile(TsvFN_X);
	Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 2);
	asserta(Fields[0] == "Label");
	const string MotifNameX = Fields[1];
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		const string &Label = Fields[0];
		const string &X = Fields[1];
		LabelToX[Label] = X;
		}
	CloseStdioFile(f);

	vector<string> X1s;
	vector<vector<string> > ABC1s;
	vector<vector<string> > ABCX1s;
	vector<string> MotifNames1;
	vector<uint> MotifLengths1;
	bool BeforeABC = 
	  ShapeSearcher::JoinABCX1(Chains, ChainLabels, LabelToABC, LabelToX,
	  X1s, ABC1s, ABCX1s, MotifNames1, MotifLengths1);
	ProgressLog("%s  %s\n", 
	  MotifNameX.c_str(),
	  BeforeABC ? "before ABC" : "after ABC");

	Shapes S1;
	Shapes S2;
	vector<string> X2s;
	vector<string> X3s;
	vector<double> ScoreX2s;
	vector<double> ScoreX3s;
	ShapeSearcher::TrainABCX(Chains, ABC1s, X1s, BeforeABC,
	  S1, S2, X2s, X3s, ScoreX2s, ScoreX3s);

	ShapeSearcher::LogTrainStats(0.4, Chains, X2s, X3s, ScoreX2s, ScoreX3s);

	Shapes S3;
	Shapes S4;
	vector<string> X4s;
	vector<string> X5s;
	vector<double> ScoreX4s;
	vector<double> ScoreX5s;
	ShapeSearcher::TrainABCX(Chains, ABC1s, X3s, BeforeABC,
	  S3, S4, X4s, X5s, ScoreX4s, ScoreX5s);

	ShapeSearcher::LogTrainStats(0.4, Chains, X4s, X5s, ScoreX4s, ScoreX5s);

	if (g_ftsv != 0)
		{
		fprintf(g_ftsv, "Label\t%s\n", MotifNameX.c_str());
		asserta(SIZE(X5s) == ChainCount);
		for (uint i = 0; i < ChainCount; ++i)
			{
			const char *Label = Chains[i]->m_Label.c_str();
			const string &X = X5s[i];
			fprintf(g_ftsv, "%s\t%s\n", Label, X.c_str());
			}
		}
	}

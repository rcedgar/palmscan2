#include "myutils.h"
#include "shapesearcher.h"

static const double MINSCOREX = 0.4;
void ShapeSearcher::TrainABCX(const vector<PDBChain *> &Chains,
  const vector<vector<string> > &ABCs,
  const vector<string> &SeqX1s,
  bool BeforeABC,
  Shapes &S1,
  Shapes &S2,
  vector<string> &SeqX2s,
  vector<string> &SeqX3s,
  vector<double> &ScoreX2s,
  vector<double> &ScoreX3s)
	{
	const uint ChainCount = SIZE(Chains);
	asserta(SIZE(ABCs) == ChainCount);
	asserta(SIZE(SeqX1s) == ChainCount);

	vector<vector<string> > ABCX1s;
	vector<uint> ShapeLengths1;
	vector<string> ShapeNames1;
	JoinABCX2(ABCs, SeqX1s, BeforeABC, ABCX1s, ShapeNames1, ShapeLengths1);

	S1.Init(ShapeNames1, ShapeLengths1);
	S1.Train(Chains, ABCX1s);

	ShapeSearcher SS1;
	SS1.Init(S1);
	ShapeSearcher::SearchABCX(S1, Chains, MINSCOREX, SeqX2s, ScoreX2s);

	vector<vector<string> > ABCX2s;
	vector<uint> ShapeLengths2;
	vector<string> ShapeNames2;
	JoinABCX2(ABCs, SeqX2s, BeforeABC, ABCX2s, ShapeNames2, ShapeLengths2);
	asserta(ShapeNames2 == ShapeNames1);
	asserta(ShapeLengths2 == ShapeLengths1);

	S2.Init(ShapeNames2, ShapeLengths2);
	S2.Train(Chains, ABCX2s);

	ShapeSearcher SS2;
	SS2.Init(S2);
	ShapeSearcher::SearchABCX(S2, Chains, MINSCOREX, SeqX3s, ScoreX3s);
	}

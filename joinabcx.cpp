#include "myutils.h"
#include "shapesearcher.h"

void ShapeSearcher::JoinABCX2(
  const vector<vector<string> > &ABCs,
  const vector<string> &SeqXs,
  bool BeforeABC,
  vector<vector<string> > &ABCXs,
  vector<string> &Names,
  vector<uint> &Lengths)
	{
	ABCXs.clear();
	Lengths.clear();
	Names.clear();

	const uint N = SIZE(ABCs);
	asserta(SIZE(ABCs) == N);
	asserta(SIZE(SeqXs) == N);

	uint LA = UINT_MAX;
	uint LB = UINT_MAX;
	uint LC = UINT_MAX;
	uint LX = UINT_MAX;
	for (uint i = 0; i < N; ++i)
		{
		asserta(SIZE(ABCs[i]) == 3);
		const string &A = ABCs[i][0];
		const string &B = ABCs[i][1];
		const string &C = ABCs[i][2];
		const string &X = SeqXs[i];

		vector<string> ABCX;
		if (BeforeABC)
			ABCX.push_back(X);
		ABCX.push_back(A);
		ABCX.push_back(B);
		ABCX.push_back(C);
		if (!BeforeABC)
			ABCX.push_back(X);
		asserta(SIZE(ABCX) == 4);
		ABCXs.push_back(ABCX);

		if (A != "" && A != ".")
			{
			if (LA == UINT_MAX)
				LA = SIZE(A);
			else
				asserta(LA == SIZE(A));
			}
		if (B != "" && B != ".")
			{
			if (LB == UINT_MAX)
				LB = SIZE(B);
			else
				asserta(LB == SIZE(B));
			}
		if (C != "" && C != ".")
			{
			if (LC == UINT_MAX)
				LC = SIZE(C);
			else
				asserta(LC == SIZE(C));
			}
		if (X != "" && X != ".")
			{
			if (LX == UINT_MAX)
				LX = SIZE(X);
			else
				asserta(LX == SIZE(X));
			}
		}

	const uint ShapeIndexX = BeforeABC ? 0 : 3;

	if (BeforeABC)
		Lengths.push_back(LX);
	Lengths.push_back(LA);
	Lengths.push_back(LB);
	Lengths.push_back(LC);
	if (!BeforeABC)
		Lengths.push_back(LX);

	if (BeforeABC)
		Names.push_back("X");
	Names.push_back("A");
	Names.push_back("B");
	Names.push_back("C");
	if (!BeforeABC)
		Names.push_back("X");

	asserta(SIZE(SeqXs) == N);
	asserta(SIZE(ABCs) == N);
	asserta(SIZE(ABCXs) == N);
	asserta(SIZE(Names) == 4);
	asserta(SIZE(Lengths) == 4);
	}

void ShapeSearcher::JoinABCX1(const vector<string> &ChainLabels,
  const map<string, vector<string> > &LabelToABC,
  const map<string, string> &LabelToX,
  bool BeforeABC,
  vector<string> &Xs,
  vector<vector<string> > &ABCs,
  vector<vector<string> > &ABCXs,
  vector<string> &Names,
  vector<uint> &Lengths)
	{
	ABCs.clear();
	Xs.clear();

	const uint N = SIZE(ChainLabels);

	for (uint i = 0; i < N; ++i)
		{
		const string &Label = ChainLabels[i];
		map<string, vector<string> >::const_iterator pABC =
		  LabelToABC.find(Label);
		string A = ".";
		string B = ".";
		string C = ".";
		if (pABC != LabelToABC.end())
			{
			const vector<string> &ABC = pABC->second;
			asserta(SIZE(ABC) == 3);
			A = ABC[0];
			B = ABC[1];
			C = ABC[2];
			}

		string X = ".";
		map<string, string>::const_iterator pX = LabelToX.find(Label);
		if (pX != LabelToX.end())
			X = pX->second;

		Xs.push_back(X);

		vector<string> ABC;
		ABC.push_back(A);
		ABC.push_back(B);
		ABC.push_back(C);
		ABCs.push_back(ABC);
		}

	JoinABCX2(ABCs, Xs, BeforeABC, ABCXs, Names, Lengths);
	}

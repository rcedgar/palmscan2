#include "myutils.h"
#include "shapesearcher.h"

void ShapeSearcher::SearchABCX(const Shapes &S,
  const vector<PDBChain *> &Chains, double MinScoreX,
  vector<string> &SeqXs, vector<double> &Scores)
	{
	const uint ChainCount = SIZE(Chains);
	SeqXs.clear();
	Scores.clear();

	ShapeSearcher SS;
	SS.Init(S);

	bool BeforeABC = false;
	uint ShapeIndexX = UINT_MAX;
	if (SS.m_ShapeIndexA == 0 && SS.m_ShapeIndexB == 1 && SS.m_ShapeIndexC == 2)
		{
		BeforeABC = false;
		ShapeIndexX = 3;
		}
	else if (SS.m_ShapeIndexA == 1 && SS.m_ShapeIndexB == 2 && SS.m_ShapeIndexC == 3)
		{
		BeforeABC = true;
		ShapeIndexX = 0;
		}
	else
		asserta(false);

	uint LX = S.GetShapeLength(ShapeIndexX);

	uint MinDist = UINT_MAX;
	uint MaxDist = UINT_MAX;
	if (BeforeABC)
		SS.GetDistRange(ShapeIndexX, SS.m_ShapeIndexA, MinDist, MaxDist);
	else
		SS.GetDistRange(SS.m_ShapeIndexC, ShapeIndexX, MinDist, MaxDist);

	vector<uint> ShapeIndexes;
	ShapeIndexes.push_back(SS.m_ShapeIndexA);
	ShapeIndexes.push_back(SS.m_ShapeIndexB);
	ShapeIndexes.push_back(SS.m_ShapeIndexC);

	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const char *Label = Chain.m_Label.c_str();
		SS.SetQuery(Chain);

		vector<string> MotifSeqs;
		MotifSeqs.push_back(".");
		MotifSeqs.push_back(".");
		MotifSeqs.push_back(".");
		MotifSeqs.push_back(".");

		double ScoreABC = SS.SearchABC();
		if (ScoreABC < SS.m_MinScoreABC)
			{
			SeqXs.push_back(".");
			Scores.push_back(0);
			continue;
			}

		vector<uint> PosVec;
		if (BeforeABC)
			PosVec.push_back(UINT_MAX);
		PosVec.push_back(SS.m_PosA);
		PosVec.push_back(SS.m_PosB);
		PosVec.push_back(SS.m_PosC);
		if (!BeforeABC)
			PosVec.push_back(UINT_MAX);

		uint Lo = UINT_MAX;
		uint Hi = UINT_MAX;

		uint QL = SS.GetQL();
		if (BeforeABC)
			{
			if (SS.m_PosA != UINT_MAX)
				{
				int iLo = int(SS.m_PosA) - int(MaxDist);
				int iHi = int(SS.m_PosA) - int(MinDist);
				if (iHi >= 0)
					{
					if (iLo < 0)
						iLo = 0;
					Hi = uint(iHi);
					Lo = uint(iLo);
					asserta(Lo <= Hi);
					}
				}
			}
		else
			{
			if (SS.m_PosC != UINT_MAX)
				{
				Lo = SS.m_PosC + MinDist;
				Hi = SS.m_PosC + MaxDist;
				}
			}

		uint PosX;
		double ScoreX;
		SS.SearchShapeTopHit(ShapeIndexX, PosVec, MinScoreX, Lo, Hi,
		  0, UINT_MAX, PosX, ScoreX);

		string SeqX = ".";
		if (PosX != UINT_MAX)
			SS.GetSubSeq(PosX, LX, SeqX);
		SeqXs.push_back(SeqX);
		Scores.push_back(ScoreX);
		}
	asserta(SIZE(SeqXs) == ChainCount);
	asserta(SIZE(Scores) == ChainCount);
	}

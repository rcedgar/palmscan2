#include "myutils.h"
#include "shapesearcher.h"

void ShapeSearcher::SearchPalm(const PDBChain &Q)
	{
	m_ScoreABC = 0;
	SetQuery(Q);
	uint QL = GetQL();
	const uint ShapeCount = GetShapeCount();

	m_ScoreABC = SearchABC();
	if (m_ScoreABC == 0)
		return;

	const Shapes &S = *m_Shapes;
	asserta(m_ShapeIndexB == m_ShapeIndexA + 1);
	asserta(m_ShapeIndexC == m_ShapeIndexB + 1);

	vector<uint> ShapeIndexOrder;
	for (int ii = int(m_ShapeIndexA) - 1; ii >= 0; --ii)
		ShapeIndexOrder.push_back(uint(ii));
	for (uint i = m_ShapeIndexC + 1; i < ShapeCount; ++i)
		ShapeIndexOrder.push_back(i);

	for (uint k = 0; k < SIZE(ShapeIndexOrder); ++k)
		{
		uint ShapeIndex = ShapeIndexOrder[k];
		uint ShapeLength = GetShapeLength(ShapeIndex);

		uint MinDist, MaxDist;
		uint Lo, Hi;
		if (ShapeIndex < m_ShapeIndexA)
			{
			uint Pos2 = m_ShapePosVec[ShapeIndex+1];
			if (Pos2 == UINT_MAX)
				continue;
			GetDistRange(ShapeIndex, ShapeIndex+1, MinDist, MaxDist);
			if (MaxDist > Pos2)
				Lo = 0;
			else
				Lo = Pos2 - MaxDist;
			if (MinDist > Pos2)
				continue;
			else
				Hi = Pos2 - MinDist;
			}
		else if (ShapeIndex > m_ShapeIndexC)
			{
			uint Pos2 = m_ShapePosVec[ShapeIndex-1];
			if (Pos2 == UINT_MAX)
				continue;
			GetDistRange(ShapeIndex-1, ShapeIndex, MinDist, MaxDist);
			Lo = Pos2 + MinDist;
			if (Lo + ShapeLength >= QL)
				continue;
			Hi = Pos2 + MaxDist;
			if (Hi + ShapeLength >= QL)
				Hi = QL - ShapeLength - 1;
			}
		else
			asserta(false);

		uint HitPos = UINT_MAX;
		double HitScore = 0;
		SearchShapeTopHit(ShapeIndex, m_ShapePosVec, 0.5, Lo, Hi,
		  0, UINT_MAX, HitPos, HitScore);

		m_ShapeScores[ShapeIndex] = HitScore;
		m_ShapePosVec[ShapeIndex] = HitPos;
		}

	m_PalmScore = 0;
	if (m_ScoreABC > 0)
		{
		vector<uint> FoundShapeIndexes;
		vector<uint> FoundPosVec;
		for (uint i = 0; i < ShapeCount; ++i)
			{
			uint Pos = m_ShapePosVec[i];
			if (Pos != UINT_MAX)
				{
				FoundShapeIndexes.push_back(i);
				FoundPosVec.push_back(Pos);
				}
			}
		asserta(SIZE(FoundPosVec) >= 3);
		m_PalmScore = GetScoreShapes(FoundShapeIndexes, FoundPosVec);
		}
	}

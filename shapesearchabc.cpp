#include "myutils.h"
#include "shapesearcher.h"
#include "sort.h"

double ShapeSearcher::SearchABC()
	{
	m_ScoreABC = 0;
	m_PosA = UINT_MAX;
	m_PosB = UINT_MAX;
	m_PosC = UINT_MAX;

	m_ShapeIndexABCs.resize(3);
	m_CharABCs.resize(3);
	m_OffsetABCs.resize(3);

	m_TopPosABCs.resize(0);
	m_TopScoreABCs.resize(0);

	m_ShapeIndexABCs[0] = m_ShapeIndexA;
	m_ShapeIndexABCs[1] = m_ShapeIndexB;
	m_ShapeIndexABCs[2] = m_ShapeIndexC;

	m_CharABCs[0] = 'D';
	m_CharABCs[1] = 'G';
	m_CharABCs[2] = 'D';

	asserta(optset_abcoffsets);
	vector<string> Fields;
	Split(opt_abcoffsets, Fields, '/');
	asserta(SIZE(Fields) == 3);
	m_OffsetABCs[0] = StrToUint(Fields[0]);
	m_OffsetABCs[1] = StrToUint(Fields[1]);
	m_OffsetABCs[2] = StrToUint(Fields[2]);

	m_TopPosABCs.resize(3);
	m_TopScoreABCs.resize(3);
	for (uint k = 0; k < 3; ++k)
		{
		SearchABC1(k);
		if (SIZE(m_TopPosABCs[k]) == 0)
			return 0;
		}

	uint MinDistAB, MaxDistAB;
	uint MinDistBC, MaxDistBC;
	GetDistRange(m_ShapeIndexA, m_ShapeIndexB, MinDistAB, MaxDistAB);
	GetDistRange(m_ShapeIndexB, m_ShapeIndexC, MinDistBC, MaxDistBC);

	MinDistAB *= 0.8;
	MaxDistAB *= 1.2;

	MinDistBC *= 0.8;
	MaxDistBC *= 1.2;

	uint NA = SIZE(m_TopPosABCs[0]);
	uint NB = SIZE(m_TopPosABCs[1]);
	uint NC = SIZE(m_TopPosABCs[2]);
	vector<uint> PosVec(3);
	for (uint ia = 0; ia < NA; ++ia)
		{
		uint PosA = m_TopPosABCs[0][ia];
		for (uint ib = 0; ib < NB; ++ib)
			{
			uint PosB = m_TopPosABCs[1][ib];
			if (PosB < PosA + MinDistAB || PosB > PosA + MaxDistAB)
				continue;
			for (uint ic = 0; ic < NC; ++ic)
				{
				uint PosC = m_TopPosABCs[2][ic];
				if (PosC < PosB + MinDistBC || PosC > PosB + MaxDistBC)
					continue;

				PosVec[0] = PosA;
				PosVec[1] = PosB;
				PosVec[2] = PosC;
				double Score = GetScoreShapes(m_ShapeIndexABCs, PosVec);
				if (Score > m_ScoreABC)
					{
					m_ScoreABC = Score;
					m_PosA = PosA;
					m_PosB = PosB;
					m_PosC = PosC;
#if 0
					{
					double PredScoreSelfA = GetSelfScore(m_ShapeIndexA, PosA);
					double PredScoreSelfB = GetSelfScore(m_ShapeIndexB, PosB);
					double PredScoreSelfC = GetSelfScore(m_ShapeIndexC, PosC);
					Log("PredA pos %u self %.3g\n", PosA, PredScoreSelfA);
					Log("PredB pos %u self %.3g\n", PosB, PredScoreSelfB);
					Log("PredC pos %u self %.3g\n", PosC, PredScoreSelfC);
					Log("all %.3g\n", Score);
					}
#endif
					}
				}
			}
		}

	m_ShapePosVec[0] = m_PosA;
	m_ShapePosVec[1] = m_PosB;
	m_ShapePosVec[2] = m_PosC;
	return m_ScoreABC;
	}

void ShapeSearcher::SearchABC1(uint k)
	{
	uint ShapeIndex = m_ShapeIndexABCs[k];
	char c = m_CharABCs[k];
	uint Offset = m_OffsetABCs[k];
	uint ShapeLength = GetShapeLength(ShapeIndex);
	m_TopPosABCs[k].resize(0);
	m_TopScoreABCs[k].resize(0);

	uint QL = GetQL();

	if (QL < 50)
		return;

	uint MinStart = 0;
	uint MaxStart = QL - ShapeLength;

	vector<uint> Hits;
	vector<double> Scores;
	SearchShapeSelf(ShapeIndex, m_MinScoreABC, MinStart, MaxStart,
	  c, Offset, Hits, Scores);
	const uint NX = SIZE(Hits);
	if (NX == 0)
		return;

	vector<uint> Order(NX);
	QuickSortOrderDesc(Scores.data(), NX, Order.data());
	for (uint i = 0; i < min(NX, m_MaxTopHitCountABC); ++i)
		{
		uint j = Order[i];
		uint Pos = Hits[j];
		double Score = Scores[j];
		m_TopPosABCs[k].push_back(Pos);
		m_TopScoreABCs[k].push_back(Score);
		}
	}

#include "myutils.h"
#include "shapesearcher.h"
#include "motifsettings.h"
#include "sort.h"

void ShapeSearcher::SearchCAB(bool DoTrace)
	{
	m_PosA = UINT_MAX;
	m_PosB = UINT_MAX;
	m_PosC = UINT_MAX;

	asserta(m_ShapeIndexA != UINT_MAX);
	asserta(m_ShapeIndexB != UINT_MAX);
	asserta(m_ShapeIndexC != UINT_MAX);

	uint QL = GetQL();
	uint AL = GetShapeLength(m_ShapeIndexA);
	uint BL = GetShapeLength(m_ShapeIndexB);
	uint CL = GetShapeLength(m_ShapeIndexC);

	if (DoTrace)
		Log("ShapeSearcher::SearchCAB >%s QL=%u\n",
		  m_Query->m_Label.c_str(), QL);

	m_ScoreABC = 0;

	m_ScoreABC = 0;
	uint L2 = 2*AL + 2*BL + 2*CL;
	if (QL < L2)
		return;

	uint MinStartC = 0;
	uint MaxStartC = QL - L2;

	vector<uint> HitsC;
	vector<double> ScoresC;
	SearchShapeSelf(m_ShapeIndexC, m_MinScoreABC, MinStartC, MaxStartC,
	  'D', g_OffCd, HitsC, ScoresC);
	const uint NC = SIZE(HitsC);

	if (DoTrace)
		Log("SearchShapeSelf(C, MinScore=%.2f, MinStartC=%u, MaxStartC=%u) %u hits\n",
		  m_MinScoreABC, MinStartC, MaxStartC, NC);

	vector<uint> ABCIndexes;
	ABCIndexes.push_back(m_ShapeIndexA);
	ABCIndexes.push_back(m_ShapeIndexB);
	ABCIndexes.push_back(m_ShapeIndexC);

	for (uint ic = 0; ic < NC; ++ic)
		{
		uint PosC = HitsC[ic];

		uint MinStartA = PosC + CL - 10;
		uint MaxStartA = PosC + 150;

		vector<uint> HitsA;
		vector<double> ScoresA;
		SearchShapeSelf(m_ShapeIndexA, m_MinScoreABC, MinStartA, MaxStartA,
		  'D', g_OffAd, HitsA, ScoresA);
		const uint NA = SIZE(HitsA);
		if (DoTrace)
			Log(" [%u] PosA=%u SearchShapeSelf(B, MinScore=%.2f, MinStartB=%u, MaxStartB=%u) %u hits\n",
			  ic, PosC, m_MinScoreABC, MinStartA, MaxStartA, NA);

		for (uint ia = 0; ia < NA; ++ia)
			{
			uint PosA = HitsA[ia];

			uint MinStartB = PosA + AL;
			uint MaxStartB = PosA + 150;

			vector<uint> HitsB;
			vector<double> ScoresB;
			SearchShapeSelf(m_ShapeIndexB, m_MinScoreABC, MinStartB, MaxStartB,
			  'G', g_OffBg, HitsB, ScoresB);
			const uint NB = SIZE(HitsB);
			if (DoTrace)
				Log("  [%u,%u] PosC,A=%u,%u SearchShapeSelf(B, MinScore=%.2f, MinStartB=%u, MaxStartB=%u) %u hits\n",
				  ic, ia, PosC, PosA, m_MinScoreABC, MinStartB, MaxStartB, NB);
			for (uint ib = 0; ib < NB; ++ib)
				{
				uint PosB = HitsB[ib];

				vector<uint> PosVec;
				PosVec.push_back(PosA);
				PosVec.push_back(PosB);
				PosVec.push_back(PosC);
				double ScoreABC = GetScoreShapes(ABCIndexes, PosVec);
				if (DoTrace)
					Log("   [%u,%u,%u] PosB %u ScoreABC %.4f\n",
					  ic, ia, ib, PosB, ScoreABC);

				if (ScoreABC > m_ScoreABC)
					{
					m_ScoreABC = ScoreABC;
					m_PosA = PosA;
					m_PosB = PosB;
					m_PosC = PosC;

					m_ShapePosVec[m_ShapeIndexA] = PosA;
					m_ShapePosVec[m_ShapeIndexB] = PosB;
					m_ShapePosVec[m_ShapeIndexC] = PosC;

					double ScoreA = ScoresA[ia];
					double ScoreB = ScoresB[ib];
					double ScoreC = ScoresC[ic];

					m_ShapeScores[m_ShapeIndexA] = ScoreA;
					m_ShapeScores[m_ShapeIndexB] = ScoreB;
					m_ShapeScores[m_ShapeIndexC] = ScoreC;
					}
				}
			}
		}
	}

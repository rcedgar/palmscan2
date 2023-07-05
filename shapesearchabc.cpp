#include "myutils.h"
#include "shapesearcher.h"
#include "motifsettings.h"
#include "sort.h"

double ShapeSearcher::SearchABC(bool DoTrace)
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
		Log("ShapeSearcher::SearchABC >%s QL=%u\n",
		  m_Query->m_Label.c_str(), QL);

	double TopScore = 0;

	uint MaxStartC = QL - CL;

	uint MinDistAC, MaxDistAC;
	GetDistRange(m_ShapeIndexA, m_ShapeIndexC, MinDistAC, MaxDistAC);

	uint MinStartA = 0;
	if (MinDistAC > MaxStartC)
		return 0;
	uint MaxStartA = MaxStartC - MinDistAC;

	vector<uint> HitsA;
	vector<double> ScoresA;
	SearchShapeSelf(m_ShapeIndexA, m_MinScoreABC, MinStartA, MaxStartA,
	  'D', g_OffAd, HitsA, ScoresA);
	const uint NA = SIZE(HitsA);

	if (DoTrace)
		Log("SearchShapeSelf(A, MinScore=%.2f, MinStartA=%u, MaxStartA=%u) %u hits\n",
		  m_MinScoreABC, MinStartA, MaxStartA, NA);

	vector<uint> ABCIndexes;
	ABCIndexes.push_back(m_ShapeIndexA);
	ABCIndexes.push_back(m_ShapeIndexB);
	ABCIndexes.push_back(m_ShapeIndexC);

	for (uint ia = 0; ia < NA; ++ia)
		{
		uint PosA = HitsA[ia];
		uint MinDistAB, MaxDistAB;
		GetDistRange(m_ShapeIndexA, m_ShapeIndexB, MinDistAB, MaxDistAB);

		uint MinStartB = PosA + MinDistAB;
		if (MinStartB >= QL - BL - CL)
			{
			if (DoTrace)
				Log("PosA %u, MinStartB > QL - BL - CL\n", PosA);
			continue;
			}
		uint MaxStartB = PosA + MaxDistAB;
		if (MaxStartB > QL - BL - CL)
			MaxStartB = QL - BL - CL;

		vector<uint> HitsB;
		vector<double> ScoresB;
		SearchShapeSelf(m_ShapeIndexB, m_MinScoreABC, MinStartB, MaxStartB,
		  'G', g_OffBg, HitsB, ScoresB);
		const uint NB = SIZE(HitsB);
		if (DoTrace)
			Log(" [%u] PosA=%u SearchShapeSelf(B, MinScore=%.2f, MinStartB=%u, MaxStartB=%u) %u hits\n",
			  ia, PosA, m_MinScoreABC, MinStartB, MaxStartB, NB);

		for (uint ib = 0; ib < NB; ++ib)
			{
			uint PosB = HitsB[ib];
			uint MinDistBC, MaxDistBC;
			GetDistRange(m_ShapeIndexB, m_ShapeIndexC, MinDistBC, MaxDistBC);

			uint MinStartC = PosB + MinDistBC;
			if (MinStartC >= QL - CL)
				{
				if (DoTrace)
					Log("PosB %u, MinStartC > QL - CL\n", PosB);
				continue;
				}
			uint MaxStartC = PosB + MaxDistBC;
			if (MaxStartC > QL - CL)
				MaxStartC = QL - CL;

			vector<uint> HitsC;
			vector<double> ScoresC;
			SearchShapeSelf(m_ShapeIndexC, m_MinScoreABC, MinStartC, MaxStartC,
			  'D', g_OffCd, HitsC, ScoresC);
			const uint NC = SIZE(HitsC);
			if (DoTrace)
				Log("  [%u,%u] PosA,B=%u,%u SearchShapeSelf(C, MinScore=%.2f, MinStartC=%u, MaxStartC=%u) %u hits\n",
				  ia, ib, PosA, PosB, m_MinScoreABC, MinStartC, MaxStartC, NC);
			for (uint ic = 0; ic < NC; ++ic)
				{
				uint PosC = HitsC[ic];

				vector<uint> PosVec;
				PosVec.push_back(PosA);
				PosVec.push_back(PosB);
				PosVec.push_back(PosC);
				double ScoreABC = GetScoreShapes(ABCIndexes, PosVec);
				if (DoTrace)
					Log("   [%u,%u,%u] PosC %u ScoreABC %.4f\n",
					  ia, ib, ic, PosC, ScoreABC);

				if (ScoreABC > TopScore)
					{
					TopScore = ScoreABC;
					m_PosA = PosA;
					m_PosB = PosB;
					m_PosC = PosC;

					m_ShapePosVec[m_ShapeIndexA] = PosA;
					m_ShapePosVec[m_ShapeIndexB] = PosB;
					m_ShapePosVec[m_ShapeIndexC] = PosC;

					m_ShapeScores[m_ShapeIndexA] = TopScore;
					m_ShapeScores[m_ShapeIndexB] = TopScore;
					m_ShapeScores[m_ShapeIndexC] = TopScore;
					}
				}
			}
		}

	return TopScore;
	}

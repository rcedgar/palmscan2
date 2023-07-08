#include "myutils.h"
#include "shapesearcher.h"
#include "motifsettings.h"
#include "sort.h"

void ShapeSearcher::SearchABC(bool DoTrace)
	{
	m_PosA = UINT_MAX;
	m_PosB = UINT_MAX;
	m_PosC = UINT_MAX;
	m_Permuted = false;

	asserta(m_ShapeIndexA != UINT_MAX);
	asserta(m_ShapeIndexB != UINT_MAX);
	asserta(m_ShapeIndexC != UINT_MAX);

	uint QL = GetQL();
	uint AL = GetShapeLength(m_ShapeIndexA);
	uint BL = GetShapeLength(m_ShapeIndexB);
	uint CL = GetShapeLength(m_ShapeIndexC);

	m_ABCScore = 0;
	if (QL < 2*AL + 2*BL + 2*CL)
		return;

	if (DoTrace)
		Log("ShapeSearcher::SearchABC >%s QL=%u\n",
		  m_Query->m_Label.c_str(), QL);


	uint MaxStartC = QL - CL;

	uint MinDistAC, MaxDistAC;
	GetDistRange(m_ShapeIndexA, m_ShapeIndexC, MinDistAC, MaxDistAC);

	uint MinStartA = 0;
	if (MinDistAC > MaxStartC)
		return;
	uint MaxStartA = MaxStartC - MinDistAC;

	vector<uint> HitsA;
	vector<double> ScoresA;
	SearchShapeSelf(m_ShapeIndexA, m_MinSelfScoreABC, MinStartA, MaxStartA,
	  'D', g_OffAd, HitsA, ScoresA);
	const uint NA = SIZE(HitsA);

	if (DoTrace)
		Log("SearchShapeSelf(A, MinScore=%.2f, MinStartA=%u, MaxStartA=%u) %u hits\n",
		  m_MinSelfScoreABC, MinStartA, MaxStartA, NA);

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
		SearchShapeSelf(m_ShapeIndexB, m_MinSelfScoreABC, MinStartB, MaxStartB,
		  'G', g_OffBg, HitsB, ScoresB);
		const uint NB = SIZE(HitsB);
		if (DoTrace)
			Log(" [%u] PosA=%u scoreA %.4f SearchShapeSelf(B, MinScore=%.2f, MinStartB=%u, MaxStartB=%u) %u hits\n",
			  ia, PosA, ScoresA[ia], m_MinSelfScoreABC, MinStartB, MaxStartB, NB);

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

		// Non-permuted
			vector<uint> HitsC;
			vector<double> ScoresC;
			double TopC1 = SearchShapeSelf(m_ShapeIndexC, m_MinSelfScoreABC,
			  MinStartC, MaxStartC, 'D', g_OffCd, HitsC, ScoresC);

		// Permuted
			if (PosA > CL)
				{
				int iMinStartC2 = int(PosA) - 40;
				if (iMinStartC2 < 0)
					iMinStartC2 = 0;
				uint MinStartC2 = uint(iMinStartC2);
				uint MaxStartC2 = PosA - CL + 8;
				vector<uint> HitsC2;
				vector<double> ScoresC2;
				double TopC2 = SearchShapeSelf(m_ShapeIndexC, m_MinSelfScoreABC,
				  MinStartC2, MaxStartC2, 'D', g_OffCd, HitsC2, ScoresC2);
				if (TopC2 > TopC1)
					{
					for (uint k = 0; k < SIZE(HitsC2); ++k)
						{
						HitsC.push_back(HitsC2[k]);
						ScoresC.push_back(ScoresC2[k]);
						}
					}
				}

			const uint NC = SIZE(HitsC);
			if (DoTrace)
				Log("  [%u,%u] PosA,B=%u,%u scoreB %.4f SearchShapeSelf(C, MinScore=%.2f, MinStartC=%u, MaxStartC=%u) %u hits\n",
				  ia, ib, PosA, PosB, ScoresB[ib], m_MinSelfScoreABC,
				  MinStartC, MaxStartC, NC);
			for (uint ic = 0; ic < NC; ++ic)
				{
				uint PosC = HitsC[ic];

				vector<uint> PosVec(m_ShapeCount, UINT_MAX);
				PosVec[m_ShapeIndexA] = PosA;
				PosVec[m_ShapeIndexB] = PosB;
				PosVec[m_ShapeIndexC] = PosC;
				double ScoreABC = GetScoreShapes(PosVec);
				if (DoTrace)
					{
					string SeqC;
					GetSubSeq(PosC, CL, SeqC);
					Log("   [%u,%u,%u] PosC %u %s scoreC %.4f ScoreABC %.4f\n",
					  ia, ib, ic, PosC, SeqC.c_str(), ScoresC[ic], ScoreABC);
					}

				if (ScoreABC > m_ABCScore)
					{
					m_ABCScore = ScoreABC;
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
	m_Permuted = (m_PosC != UINT_MAX && m_PosA != UINT_MAX && m_PosC < m_PosA);
	}

#include "myutils.h"
#include "shapesearcher.h"

#define TRACE	0

void ShapeSearcher::Init(const Shapes &S)
	{
	m_Shapes = &S;
	m_ShapeCount = S.GetShapeCount();
	m_ShapeIndexA = UINT_MAX;
	m_ShapeIndexB = UINT_MAX;
	m_ShapeIndexC = UINT_MAX;
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		if (S.m_Names[i] == "A")
			m_ShapeIndexA = i;
		else if (S.m_Names[i] == "B")
			m_ShapeIndexB = i;
		else if (S.m_Names[i] == "C")
			m_ShapeIndexC = i;
		}
	}

void ShapeSearcher::SetQuery(const PDBChain &Query)
	{
	ClearSearch();
	m_Query = &Query;

	const uint ShapeCount = GetShapeCount();
	m_ShapePosVec.resize(ShapeCount);
	m_ShapeScores.resize(ShapeCount);
	for (uint i = 0; i < ShapeCount; ++i)
		{
		m_ShapePosVec[i] = UINT_MAX;
		m_ShapeScores[i] = 0;
		}
	}

void ShapeSearcher::GetShapeSeq(uint ShapeIndex, string &Seq) const
	{
	uint ShapeLength = GetShapeLength(ShapeIndex);
	asserta(ShapeIndex < SIZE(m_ShapePosVec));
	uint Pos = m_ShapePosVec[ShapeIndex];
	m_Query->GetSubSeq(Pos, ShapeLength, Seq);
	}

// dist = start_j - start_i
void ShapeSearcher::GetDistRange(uint ShapeIndexi, uint ShapeIndexj,
  uint &MinDist, uint &MaxDist) const
	{
	asserta(ShapeIndexj > ShapeIndexi);

	uint QL = GetQL();
	MinDist = 0;
	MaxDist = 0;
	for (uint ShapeIndex = ShapeIndexi; ShapeIndex < ShapeIndexj;
	  ++ShapeIndex)
		{
		MinDist += m_Shapes->m_MinNeighborDists[ShapeIndex];
		MaxDist += m_Shapes->m_MaxNeighborDists[ShapeIndex];
		}
	asserta(MinDist < QL);
	asserta(MaxDist < QL);
	}

double ShapeSearcher::GetScoreShapes(const vector<uint> &ShapeIndexes,
  const vector<uint> &PosVec) const
	{
	double Sum = 0;
	const uint N = SIZE(ShapeIndexes);
	if (N == 0)
		return 0;
	uint Pairs = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint ShapeIndexi = ShapeIndexes[i];
		uint Posi = PosVec[i];
		for (uint j = 0; j <= i; ++j)
			{
			uint ShapeIndexj = ShapeIndexes[j];
			uint Posj = PosVec[j];
			Sum += GetScoreShapePair(ShapeIndexi, ShapeIndexj, Posi, Posj);
			++Pairs;
			}
		}
	return Sum/Pairs;
	}

double ShapeSearcher::GetScoreShapePair(uint ShapeIndex1, uint ShapeIndex2,
  uint Pos1, uint Pos2) const
	{
	double Sum = 0;
	uint L1 = GetShapeLength(ShapeIndex1);
	uint L2 = GetShapeLength(ShapeIndex2);
	uint QL = GetQL();
	asserta(Pos1 + L1 <= QL);
	asserta(Pos2 + L2 <= QL);
	uint n = 0;
	for (uint Offset1 = 0; Offset1 < L1; ++Offset1)
		{
		uint StartOffset2 = (ShapeIndex1 == ShapeIndex2 ? Offset1 + 1 : 0);
		for (uint Offset2 = StartOffset2; Offset2 < L2; ++Offset2)
			{
			Sum += GetScoreResiduePair(ShapeIndex1, ShapeIndex2,
			  Pos1, Pos2, Offset1, Offset2);
			++n;
			}
		}
	assert(n > 0);
	double Score = Sum/n;
	return Score;
	}

double ShapeSearcher::GetScoreResiduePair(uint ShapeIndex1, uint ShapeIndex2,
  uint Pos1, uint Pos2, uint Offset1, uint Offset2) const
	{
	bool Diag = (ShapeIndex1 == ShapeIndex2);
	double XS = (Diag ? 1.5 : 2);
	const PDBChain &Query = *m_Query;
	double Observed_d = Query.GetDist(Pos1+Offset1, Pos2+Offset2);
	Observed_d = fabs(Observed_d);
	double Mu = GetMeanDist3(ShapeIndex1, ShapeIndex2, Offset1, Offset2);
	Mu = fabs(Mu);
	double Sigma = GetStdDev3(ShapeIndex1, ShapeIndex2, Offset1, Offset2);
	double y = GetNormal(Mu, XS*Sigma, Observed_d);
	double Max = GetNormal(Mu, XS*Sigma, Mu);
	double Score = y/Max;
	return Score;
	}

double ShapeSearcher::GetSelfScore(uint ShapeIndex, uint Pos) const
	{
	double Score = GetScoreShapePair(ShapeIndex, ShapeIndex, Pos, Pos);
	return Score;
	}

double ShapeSearcher::GetScore(uint ShapeIndex, uint Pos,
  const vector<uint> &PosVec) const
	{
	const uint ShapeCount = GetShapeCount();
	asserta(SIZE(PosVec) == ShapeCount);
	double Sum = 0;
	uint n = 0;
	for (uint ShapeIndex2 = 0; ShapeIndex2 < ShapeCount; ++ShapeIndex2)
		{
		uint Pos2 = PosVec[ShapeIndex2];
		if (Pos2 == UINT_MAX)
			continue;
		++n;
		Sum += GetScoreShapePair(ShapeIndex, ShapeIndex2, Pos, Pos2);
		}
	asserta(n > 0);
	double Score = Sum/n;
	return Score;
	}

void ShapeSearcher::SearchShapeSelf(uint ShapeIndex, double MinScore,
  uint Lo, uint Hi, char Letter, uint LetterOffset,
  vector<uint> &HitPosVec, vector<double> &HitScores) const
	{
	HitPosVec.clear();
	HitScores.clear();
	uint L = GetShapeLength(ShapeIndex);
	uint QL = GetQL();
	asserta(Hi + L <= QL);
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		if (Letter != 0 && LetterOffset != UINT_MAX &&
		  m_Query->m_Seq[Pos+LetterOffset] != Letter)
			continue;
		double Score = GetSelfScore(ShapeIndex, Pos);
		//Log("Pos=%u score=%.3g\n", Pos, Score);
		if (Score < MinScore)
			continue;
		uint HitCount = SIZE(HitPosVec);
		if (HitCount > 0)
			{
			uint LastHit = HitPosVec[HitCount-1];
			if (Pos - LastHit < L)
				{
				double LastScore = HitScores[HitCount-1];
				if (Score > LastScore)
					{
					HitPosVec[HitCount-1] = Pos;
					HitScores[HitCount-1] = Score;
					}
				continue;
				}
			}
		HitPosVec.push_back(Pos);
		HitScores.push_back(Score);
		}
	}

void ShapeSearcher::SearchShapeTopHit(uint ShapeIndex,
  const vector<uint> &PosVec, double MinScore, uint Lo, uint Hi,
  char Letter, uint LetterOffset, uint &Pos, double &Score) const
	{
	if (Lo == UINT_MAX)
		return;
	asserta(Hi != UINT_MAX);

	Pos = UINT_MAX;
	Score = 0;

	vector<uint> HitPosVec;
	vector<double> Scores;
	SearchShape(ShapeIndex, PosVec, MinScore, Lo, Hi,
	  Letter, LetterOffset, HitPosVec, Scores);

	const uint N = SIZE(HitPosVec);
	for (uint i = 0; i < N; ++i)
		{
		if (Scores[i] > Score)
			{
			Pos = HitPosVec[i];
			Score = Scores[i];
			}
		}
	}

void ShapeSearcher::SearchShape(uint ShapeIndex, const vector<uint> &PosVec,
  double MinScore, uint Lo, uint Hi, char Letter, uint LetterOffset,
  vector<uint> &HitPosVec, vector<double> &HitScores) const
	{
	HitPosVec.clear();
	HitScores.clear();
	if (Lo == UINT_MAX)
		return;
	asserta(Hi != UINT_MAX);
	uint L = GetShapeLength(ShapeIndex);
	uint QL = GetQL();
	asserta(Hi + L <= QL);
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		if (Letter != 0 && LetterOffset != UINT_MAX &&
		  m_Query->m_Seq[Pos+LetterOffset] != Letter)
			continue;
		double Score = GetScore(ShapeIndex, Pos, PosVec);
		if (Score < MinScore)
			continue;
		uint HitCount = SIZE(HitPosVec);
		if (HitCount > 0)
			{
			uint LastHit = HitPosVec[HitCount-1];
			if (Pos - LastHit < L)
				{
				double LastScore = HitScores[HitCount-1];
				if (Score > LastScore)
					{
					HitPosVec[HitCount-1] = Pos;
					HitScores[HitCount-1] = Score;
					}
				continue;
				}
			}
		HitPosVec.push_back(Pos);
		HitScores.push_back(Score);
		}
	}

double ShapeSearcher::SearchABC()
	{
	asserta(m_ShapeIndexA != UINT_MAX);
	asserta(m_ShapeIndexB != UINT_MAX);
	asserta(m_ShapeIndexC != UINT_MAX);

	uint QL = GetQL();
	uint AL = GetShapeLength(m_ShapeIndexA);
	uint BL = GetShapeLength(m_ShapeIndexB);
	uint CL = GetShapeLength(m_ShapeIndexC);

	uint Offset_D_MotifA = m_Shapes->m_Offset_D_MotifA;
	uint Offset_G_MotifB = m_Shapes->m_Offset_G_MotifB;
	uint Offset_D_MotifC = m_Shapes->m_Offset_D_MotifC;

#if TRACE
	Log("ShapeSearcher::SearchABC >%s QL=%u\n", m_Query->m_Label.c_str(), QL);
#endif

	double TopScore = 0;

	uint MaxStartC = QL - CL;

	uint MinDistAC, MaxDistAC;
	GetDistRange(m_ShapeIndexA, m_ShapeIndexC, MinDistAC, MaxDistAC);

	uint MinStartA = 0;
	uint MaxStartA = MaxStartC - MinDistAC;

	vector<uint> HitsA;
	vector<double> ScoresA;
	SearchShapeSelf(m_ShapeIndexA, m_MinScoreABC, MinStartA, MaxStartA,
	  'D', Offset_D_MotifA, HitsA, ScoresA);
	const uint NA = SIZE(HitsA);
#if TRACE
	Log("SearchShapeSelf(A, MinScore=%.2f, MinStartA=%u, MaxStartA=%u) %u hits\n",
	  m_MinScoreABC, MinStartA, MaxStartA, NA);
#endif

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
#if TRACE
			Log("PosA %u, MinStartB > QL - BL - CL\n", PosA);
#endif
			continue;
			}
		uint MaxStartB = PosA + MaxDistAB;
		if (MaxStartB > QL - BL - CL)
			MaxStartB = QL - BL - CL;

		vector<uint> HitsB;
		vector<double> ScoresB;
		SearchShapeSelf(m_ShapeIndexB, m_MinScoreABC, MinStartB, MaxStartB,
		  'G', Offset_G_MotifB, HitsB, ScoresB);
		const uint NB = SIZE(HitsB);
#if TRACE
		Log(" [%u] PosA=%u SearchShapeSelf(B, MinScore=%.2f, MinStartB=%u, MaxStartB=%u) %u hits\n",
		  ia, PosA, m_MinScoreABC, MinStartB, MaxStartB, NB);
#endif

		for (uint ib = 0; ib < NB; ++ib)
			{
			uint PosB = HitsB[ib];
			uint MinDistBC, MaxDistBC;
			GetDistRange(m_ShapeIndexB, m_ShapeIndexC, MinDistBC, MaxDistBC);

			uint MinStartC = PosB + MinDistBC;
			if (MinStartC >= QL - CL)
				{
#if TRACE
				Log("PosB %u, MinStartC > QL - CL\n", PosB);
#endif
				continue;
				}
			uint MaxStartC = PosB + MaxDistBC;
			if (MaxStartC > QL - CL)
				MaxStartC = QL - CL;

			vector<uint> HitsC;
			vector<double> ScoresC;
			SearchShapeSelf(m_ShapeIndexC, m_MinScoreABC, MinStartC, MaxStartC,
			  'D', Offset_D_MotifC, HitsC, ScoresC);
			const uint NC = SIZE(HitsC);
#if TRACE
			Log("  [%u,%u] PosA,B=%u,%u SearchShapeSelf(C, MinScore=%.2f, MinStartC=%u, MaxStartC=%u) %u hits\n",
			  ia, ib, PosA, PosB, m_MinScoreABC, MinStartC, MaxStartC, NC);
#endif
			for (uint ic = 0; ic < NC; ++ic)
				{
				uint PosC = HitsC[ic];

				vector<uint> PosVec;
				PosVec.push_back(PosA);
				PosVec.push_back(PosB);
				PosVec.push_back(PosC);
				double ScoreABC = GetScoreShapes(ABCIndexes, PosVec);
#if TRACE
				Log("   [%u,%u,%u] PosC %u ScoreABC %.4f\n",
				  ia, ib, ic, PosC, ScoreABC);
#endif

				if (ScoreABC > TopScore)
					{
					TopScore = ScoreABC;
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

#if TRACE
	Log("TopScore %.4f PosA,B,C %u,%u,%u\n",
	  TopScore, TopPosA, TopPosB, TopPosC);
#endif
	return TopScore;
	}

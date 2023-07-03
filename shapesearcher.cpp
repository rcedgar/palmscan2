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

void ShapeSearcher::GetSubSeq(uint Pos, uint n, string &Seq) const
	{
	m_Query->GetSubSeq(Pos, n, Seq);
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

	MinDist = 0;
	MaxDist = 0;
	for (uint ShapeIndex = ShapeIndexi; ShapeIndex < ShapeIndexj;
	  ++ShapeIndex)
		{
		MinDist += m_Shapes->m_MinNeighborDists[ShapeIndex];
		MaxDist += m_Shapes->m_MaxNeighborDists[ShapeIndex];
		}
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
	if (Pos1 + L1 > QL)
		return 0;
	if (Pos2 + L2 > QL)
		return 0;
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
	if (Lo == UINT_MAX || Hi == UINT_MAX)
		return;
	HitPosVec.clear();
	HitScores.clear();
	uint L = GetShapeLength(ShapeIndex);
	uint QL = GetQL();
	if (Hi + L > QL)
		return;
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
	Pos = UINT_MAX;
	Score = 0;
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
	//asserta(Hi + L <= QL);
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

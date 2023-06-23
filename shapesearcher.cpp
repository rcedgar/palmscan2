#include "myutils.h"
#include "shapesearcher.h"

void ShapeSearcher::SetQuery(const PDBChain &Query,
  uint PosA, uint PosB, uint PosC)
	{
	ClearSearch();
	if (PosA == UINT_MAX || PosB == UINT_MAX || PosC == UINT_MAX)
		return;
	m_Query = &Query;
	m_PosA = PosA;
	m_PosB = PosB;
	m_PosC = PosC;
	}

void ShapeSearcher::GetLoHi(uint ShapeIndex, uint ShapeIndex2,
  uint Pos2, uint &Lo, uint &Hi) const
	{
	if (ShapeIndex == ShapeIndex2 && Pos2 != UINT_MAX)
		{
		Lo = Pos2;
		Hi = Pos2;
		return;
		}

	Lo = UINT_MAX;
	Hi = UINT_MAX;
	if (Pos2 == UINT_MAX)
		return;

	const Shapes &S = *m_Shapes;

	uint QL = GetQL();
	double MeanDist = S.GetMeanDist2(ShapeIndex2, ShapeIndex);
	double StdDev = S.GetStdDev2(ShapeIndex2, ShapeIndex);

	double Var = m_Sigmas*StdDev;
	double dLo = double(Pos2) + MeanDist - Var - 0.5;
	double dHi = double(Pos2) + MeanDist + Var + 0.5;
	Lo = (dLo < 0 ? 0 : uint(dLo));
	if (Lo >= QL)
		Lo = QL - 1;
	Hi = uint(dHi);
	if (Hi >= QL)
		Hi = QL - 1;
	if (Lo > Hi)
		{
		Lo = UINT_MAX;
		Hi = UINT_MAX;
		}
#if 0
	{
	Log("%s  %s %u + { %+.3g +/- %.3g } = %u .. %u\n",
	  GetShapeName(ShapeIndex),
	  GetShapeName(ShapeIndex2),
	  Pos2,
	  MeanDist, StdDev,
	  Lo, Hi);
	}
#endif
	}

void ShapeSearcher::GetMinLoMaxHi(uint ShapeIndex,
  const vector<uint> &PosVec, uint &MinLo, uint &MaxHi) const
	{
	uint ShapeCount = GetShapeCount();
	asserta(SIZE(PosVec) == ShapeCount);
	MinLo = UINT_MAX;
	MaxHi = UINT_MAX;
	uint Pos = PosVec[ShapeIndex];
	if (Pos != UINT_MAX)
		{
		MinLo = Pos;
		MaxHi = Pos;
		return;
		}

	uint FirstShapeIndex = UINT_MAX;
	uint SecondShapeIndex = UINT_MAX;
	for (uint ShapeIndex2 = 0; ShapeIndex2 < ShapeCount; ++ShapeIndex2)
		{
		uint Pos = PosVec[ShapeIndex2];
		if (Pos == UINT_MAX)
			continue;
		if (ShapeIndex2 < ShapeIndex)
			FirstShapeIndex = ShapeIndex2;
		if (SecondShapeIndex == UINT_MAX && ShapeIndex2 > ShapeIndex)
			{
			SecondShapeIndex = ShapeIndex2;
			break;
			}
		}
	if (FirstShapeIndex == UINT_MAX && SecondShapeIndex != UINT_MAX)
		FirstShapeIndex = SecondShapeIndex;
	if (FirstShapeIndex != UINT_MAX && SecondShapeIndex == UINT_MAX)
		SecondShapeIndex = FirstShapeIndex;

	if (FirstShapeIndex == UINT_MAX || SecondShapeIndex == UINT_MAX)
		return;

	uint FirstPos = PosVec[FirstShapeIndex];
	asserta(FirstPos != UINT_MAX);
	uint Hi;
	GetLoHi(ShapeIndex, FirstShapeIndex, FirstPos, MinLo, Hi);

	uint SecondPos = PosVec[FirstShapeIndex];
	asserta(FirstPos != UINT_MAX);
	uint Lo;
	GetLoHi(ShapeIndex, SecondShapeIndex, FirstPos, Lo, MaxHi);
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

void ShapeSearcher::SearchOneShapeSelf(uint ShapeIndex, double MinScore,
  uint Lo, uint Hi, char Letter, uint LetterOffset,
  vector<uint> &HitPosVec, vector<double> &HitScores) const
	{
	uint L = GetShapeLength(ShapeIndex);
	uint QL = GetQL();
	asserta(Hi + L <= QL);
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		if (Letter != 0 && LetterOffset != UINT_MAX &&
		  m_Query->m_Seq[Pos] != Letter)
			continue;
		double Score = GetSelfScore(ShapeIndex, Pos);
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

#pragma once

#include "shapes.h"

class ShapeSearcher
	{
public:
	const PDBChain *m_Query = 0;
	uint m_PosA = UINT_MAX;
	uint m_PosB = UINT_MAX;
	uint m_PosC = UINT_MAX;
	const Shapes *m_Shapes = 0;
	uint m_ShapeCount = UINT_MAX;
	double m_Sigmas = 2.5;
	uint m_ShapeIndexA = UINT_MAX;
	uint m_ShapeIndexB = UINT_MAX;
	uint m_ShapeIndexC = UINT_MAX;
	double m_MinScoreABC = 0.4;

public:
	void ClearSearch()
		{
		m_Query = 0;
		m_PosA = UINT_MAX;
		m_PosB = UINT_MAX;
		m_PosC = UINT_MAX;
		}

	void Init(const Shapes &S);

	uint GetQL() const
		{
		return m_Query->GetSeqLength();\
		}

	const char *GetShapeName(uint ShapeIndex) const
		{
		return m_Shapes->m_Names[ShapeIndex].c_str();
		}

	uint GetShapeLength(uint ShapeIndex) const
		{
		return m_Shapes->m_Lengths[ShapeIndex];
		}

	double GetMeanDist3(uint ShapeIndex1, uint ShapeIndex2,
	  uint Offset1, uint Offset2) const
		{
		return m_Shapes->GetMeanDist3(ShapeIndex1, ShapeIndex2, Offset1, Offset2);
		}

	double GetStdDev3(uint ShapeIndex1, uint ShapeIndex2,
	  uint Offset1, uint Offset2) const
		{
		return m_Shapes->GetStdDev3(ShapeIndex1, ShapeIndex2, Offset1, Offset2);
		}

	uint GetShapeCount() const { return m_ShapeCount; }

	void SetQuery(const PDBChain &Query, uint PosA, uint PosB, uint PosC);
	void SetQuery(const PDBChain &Query);

	void GetDistRange(uint ShapeIndex, uint ShapeIndex2, 
	  uint &MinDist, uint &MaxDist) const;

	void SearchShape(uint ShapeIndex, const vector<uint> &PosVec,
	  double MinScore, uint Lo, uint Hi, char Letter, uint LetterOffset,
	  vector<uint> &HitPosVec, vector<double> &HitScores) const;

	void SearchShapeSelf(uint ShapeIndex, double MinScore,
	  uint Lo, uint Hi, char Letter, uint LetterOffset,
	  vector<uint> &HitPosVec, vector<double> &HitScores) const;

	double GetSelfScore(uint ShapeIndex, uint Pos) const;

	double GetScore(uint ShapeIndex, uint Pos,
	  const vector<uint> &PosVec) const;

	double GetScoreShapePair(uint ShapeIndex1, uint ShapeIndex2,
	  uint Pos1, uint Pos2) const;

	double GetScoreShapes(const vector<uint> &ShapeIndexes,
	  const vector<uint> &PosVec) const;

	double GetScoreResiduePair(uint ShapeIndex1, uint ShapeIndex2,
	  uint Pos1, uint Pos2, uint Offset1, uint Offset2) const;

	double SearchABC(uint &PosA, uint &PosB, uint &PosC);

	void TestABC1(const PDBChain &Chain,
	  const vector<string> &MotifSeqs);

public:
	static void TestABC(const Shapes &S,
	  const vector<PDBChain *> &Chains,
	  vector<vector<string> > &MotifSeqsVec);
	};

double GetNormal(double Mu, double Sigma, double x);

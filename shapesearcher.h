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

public:
	void Init(const Shapes &S)
		{
		m_Shapes = &S;
		m_ShapeCount = S.GetShapeCount();
		}

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

	double GetMeanDist2(uint ShapeIndex1, uint ShapeIndex2) const
		{
		return m_Shapes->GetMeanDist2(ShapeIndex1, ShapeIndex2);
		}

	double GetStdDev2(uint ShapeIndex1, uint ShapeIndex2) const
		{
		return m_Shapes->GetStdDev2(ShapeIndex1, ShapeIndex2);
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

	void ClearSearch() {  }
	void Search(const PDBChain &Query, uint PosA, uint PosB, uint PosC);
	uint GetShapeCount() const { return m_ShapeCount; }
	void GetMinLoMaxHi(uint ShapeIndex, const vector<uint> &PosVec,
	   uint &MinLo, uint &MaxHi) const;
	void GetLoHi(uint ShapeIndex, uint ShapeIndex2, 
	  uint Pos2, uint &Lo, uint &Hi) const;

	double GetScore(uint ShapeIndex, uint Pos,
	  const vector<uint> &OtherShapeIndexes) const;

	double GetScoreShapePair(uint ShapeIndex1, uint ShapeIndex2,
	  uint Pos1, uint Pos2) const;

	double GetScoreResiduePair(uint ShapeIndex1, uint ShapeIndex2,
	  uint Pos1, uint Pos2, uint Offset1, uint Offset2) const;
	};

double GetNormal(double Mu, double Sigma, double x);

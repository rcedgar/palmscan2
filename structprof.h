#pragma once

#include "pdbchain.h"

class StructProf
	{
public:
	const PDBChain *m_Chain = 0;
	uint m_MinPos = UINT_MAX;
	uint m_MaxPos = UINT_MAX;

public:
	void Clear()
		{
		m_Chain = 0;
		m_MinPos = UINT_MAX;
		m_MaxPos = UINT_MAX;
		}

	void SetChain(const PDBChain &Chain);
	void SetMinMaxPos(uint MinPos, uint MaxPos);
	void GetHSE(uint Pos, double Radius,
	  uint &NU, uint &ND) const;
	uint GetTSB(uint Pos, double Radius) const;
	void GetSphere(const vector<double> &CenterPt, double Radius,
	  vector<uint> &PosVec) const;
	uint SearchDist(uint Pos, uint Lo, uint Hi,
	  bool Maximize, double X, double &BestDist) const;
	uint FindMofifD_Hueuristics() const;
	uint FindMofifE_Hueuristics(uint Pos_MotifD) const;
	};

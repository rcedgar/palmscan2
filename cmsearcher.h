#pragma once

#include "cmpsearcher.h"
#include "pdbchain.h"
#include "cmp.h"

class DSHit;

class CMSearcher
	{
public:
	static vector<const CMP *> m_Profs;

public:
	CMPSearcher m_PS;

	PDBChain *m_Query = 0;

	uint m_RefIndex = UINT_MAX;
	uint m_PosA = UINT_MAX;
	uint m_PosB = UINT_MAX;
	uint m_PosC = UINT_MAX;
	double m_Score = 0;

public:
	void ClearSearch()
		{
		m_Query = 0;
		m_RefIndex = UINT_MAX;
		uint m_PosA = UINT_MAX;
		uint m_PosB = UINT_MAX;
		uint m_PosC = UINT_MAX;
		double m_Score = 0;
		}

	void Search(PDBChain &Query);
	double GetPSSMStarts(uint &PosA, uint &PosB, uint &PosC) const;

public:
	static void ProfsFromFile(const string &FileName);
	static bool MeansFromFile(FILE *f, string &Label,
	  vector<vector<double> > &Means);
	};

#pragma once

#include "ppspsearcher.h"
#include "pdbchain.h"
#include "ppsp.h"

class DSHit;

class PPSPDBSearcher
	{
public:
	vector<PPSP *> m_Profs;
	PPSPSearcher m_PS;

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
	void ProfsFromFile(const string &FileName);
	};

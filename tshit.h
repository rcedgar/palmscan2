#pragma once

#include "trisearcher.h"
#include "pdbchain.h"

class TSHit
	{
public:
	const PDBChain *m_Query = 0;
	const PDBChain *m_Ref = 0;
	uint m_QPosA = UINT_MAX;
	uint m_QPosB = UINT_MAX;
	uint m_QPosC = UINT_MAX;

	uint m_RPosA = UINT_MAX;
	uint m_RPosB = UINT_MAX;
	uint m_RPosC = UINT_MAX;

	double m_TriRMSD2 = DBL_MAX;
	double m_MotifRMSD2 = DBL_MAX;

public:
	void Clear()
		{
		m_Query = 0;
		m_Ref = 0;
		m_QPosA = UINT_MAX;
		m_QPosB = UINT_MAX;
		m_QPosC = UINT_MAX;

		m_RPosA = UINT_MAX;
		m_RPosB = UINT_MAX;
		m_RPosC = UINT_MAX;

		m_TriRMSD2 = DBL_MAX;
		m_MotifRMSD2 = DBL_MAX;
		}

	void WriteAln(FILE *f) const;
	void WriteTsv(FILE *f) const;
	void WritePalmprintFasta(FILE *f) const;
	void WritePalmprintPDB(const string &FileNamePrefix) const;
	void WriteSketch(FILE* f) const;
	};

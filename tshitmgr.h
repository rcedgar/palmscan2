#pragma once

#include "trisearcher.h"

class TSHitMgr
	{
public:
	const PDBChain *m_Query = 0;
	vector<TSHit> m_Hits;
	TSHit *m_TopHit = 0;

public:
	void SetQuery(const PDBChain &Query);
	void Add(TSHit &TH);
	void WriteOutput();
	void SetTopHit();
	void WriteReport(FILE *f) const;
	void WritePPC(FILE* f) const;
	};

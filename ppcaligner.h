#pragma once

#include "pdbchain.h"
#include "tshit.h"

class PpcAligner
	{
	const PDBChain *m_Q = 0;
	const PDBChain *m_R = 0;

	vector<double> m_RPt;
	vector<double> m_QPt;

public:
	PpcAligner()
		{
		Clear();
		}

	void Clear()
		{
		m_Q = 0;
		m_R = 0;
		m_RPt.clear();
		m_QPt.clear();
		m_RPt.resize(3);
		m_QPt.resize(3);
		}

	void SetQuery(const PDBChain &Q);
	void SetRef(const PDBChain &R);
	double GetRMSD2Segment(uint QPos, uint RPos, uint n);
	double GetMotifRMSD();
	};

#pragma once

#include "cdinfo.h"

class CDData
	{
public:
	const CDInfo *m_Info;
	vector<vector<double> > m_Data;

public:
	void Clear()
		{
		m_Info = 0;
		m_Data.clear();
		}

	void Init(const CDInfo &Info);
	double Get(uint MotifIndex, uint i) const;
	void Set(uint MotifIndex, uint i, double Value);
	};

#pragma once

#include "myutils.h"
#include "pdbchain.h"

class DSS
	{
public:
	const PDBChain *m_Chain = 0;
	uint m_L = 0;

	vector<double> m_NUDX_ScaledValues;
	string m_SS;

	int m_NUDX_W = 50;
	double m_NUDX_Radius = 20.0;
	int m_SSD2_W = 100;
	int m_SSD2_w = 12;
	double m_SSD2_Step1 = 5;
	double m_SSD2_Step2 = 7;
	double m_SSD2_Step3 = 9;

public:
	static vector<vector<double> > g_BinLos;
	static vector<vector<double> > g_Scores;

public:
	void Init(const PDBChain &Chain)
		{
		m_Chain = &Chain;
		m_NUDX_ScaledValues.clear();
		m_SS.clear();
		}
	uint GetSeqLength() const { return m_Chain->GetSeqLength(); }

	uint GetFeature(uint FeatureIndex, uint Pos);
	double Get_NUDX(uint Pos);
	char Get_SSX(uint Pos);
	char Get_SSX2(uint Pos);
	double Get_SSD2(uint Pos);
	void Get_NUDX_Lo(uint Pos, double &NU, double &ND) const;
	void Set_NUDXVec();

public:
	static const char *GetFeatureName(uint FeatureIndex);
	};

#pragma once

#include "myutils.h"
#include "pdbchain.h"

const uint DSSFeatureIndex = 99;

// Discrete Structure States
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
	static uint m_NUDX_Bins;
	static string m_AlphaStr;
	static byte m_AminoLetterToCompressedLetter[20];
	static uint m_CompressedAminoAlphaSize;

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
	static uint GetAlphaSize(uint FeatureIndex);
	static void SetCharToLetter();
	static double GetScore(
		uint i0, uint i1, uint i2, uint i3, uint i4,
		uint j0, uint j1, uint j2, uint j3, uint j4);
	static uint GetIdx(uint i0, uint i1, uint i2, uint i3, uint i4);
	static double GetExpectedScore(const vector<vector<double> > &ObsFreqs);
	};

#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "alpha.h"
#include "dss.h"

uint DSS::m_NUDX_Bins = 8;
string DSS::m_AlphaStr = "AST,C,DN,EHKQR,FWY,G,ILMV,P";
byte DSS::m_AminoLetterToCompressedLetter[20];
uint DSS::m_CompressedAminoAlphaSize;

static bool StaticInit()
	{
	DSS::SetCharToLetter();
	return true;
	}
static bool g_InitDone = StaticInit();

void DSS::SetCharToLetter()
	{
	vector<string> Fields;
	Split(m_AlphaStr, Fields, ',');
	m_CompressedAminoAlphaSize = SIZE(Fields);
	for (uint i = 0; i < 20; ++i)
		m_AminoLetterToCompressedLetter[i] = INVALID_LETTER;

	for (uint Group = 0; Group < m_CompressedAminoAlphaSize; ++Group)
		{
		const string &GroupStr = Fields[Group];
		for (uint i = 0; i < SIZE(GroupStr); ++i)
			{
			byte c = GroupStr[i];
			uint Letter = g_CharToLetterAmino[c];
			m_AminoLetterToCompressedLetter[Letter] = Group;
			}
		}
	for (uint i = 0; i < 20; ++i)
		if (m_AminoLetterToCompressedLetter[i] == INVALID_LETTER)
			Die("Letter %c missing in DSS::m_AlphaStr='%s'",
			  g_LetterToCharAmino[i], m_AlphaStr.c_str());
	}

char DSS::Get_SSX(uint Pos)
	{
	if (m_SS.empty())
		m_Chain->GetSS(m_SS);
	asserta(Pos < SIZE(m_SS));
	return m_SS[Pos];
	}

double DSS::Get_NUDX(uint Pos)
	{
	if (m_NUDX_ScaledValues.empty())
		Set_NUDXVec();
	asserta(Pos < SIZE(m_NUDX_ScaledValues));
	return m_NUDX_ScaledValues[Pos];
	}

void DSS::Set_NUDXVec()
	{
	const uint L = GetSeqLength();
	vector<double> Values;
	m_NUDX_ScaledValues.clear();
	double MinValue = 999;
	double MaxValue = 0;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		double NU, ND;
		Get_NUDX_Lo(Pos, NU, ND);
		double Value = (NU == DBL_MAX || ND == DBL_MAX) ? DBL_MAX : NU + ND;
		Values.push_back(Value);
		if (Value != DBL_MAX)
			{
			MinValue = min(MinValue, Value);
			MaxValue = max(MaxValue, Value);
			}
		}

	double Range = (MaxValue - MinValue);
	if (Range < 1)
		Range = 1;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		double Value = Values[Pos];
		if (Value == DBL_MAX)
			{
			m_NUDX_ScaledValues.push_back(DBL_MAX);
			continue;
			}
		double ScaledValue = (Value - MinValue)/Range;
		asserta(ScaledValue >= 0 && ScaledValue <= 1);
		m_NUDX_ScaledValues.push_back(ScaledValue);
		}
	}

void DSS::Get_NUDX_Lo(uint Pos, double &NU, double &ND) const
	{
	const double Radius = 12.0;
	NU = 0;
	ND = 0;
	const PDBChain &Chain = *m_Chain;
	const uint L = SIZE(Chain.m_Seq);
	if (Pos == 0 || Pos+1 >= L)
		{
		NU = DBL_MAX;
		ND = DBL_MAX;
		return;
		}

	vector<double> PtPrevCA;
	vector<double> PtNextCA;
	vector<double> PtCA;
	vector<double> PtCB;
	Chain.GetPt(Pos-1, PtPrevCA);
	Chain.GetPt(Pos, PtCA);
	Chain.GetPt(Pos+1, PtNextCA);

	vector<double> d1;
	vector<double> d2;
	Sub_Vecs(PtCA, PtPrevCA, d1);
	Sub_Vecs(PtCA, PtNextCA, d2);

	vector<double> VecPAB;
	Add_Vecs(d1, d2, VecPAB);
	NormalizeVec(VecPAB);

	vector<double> Pt2;
	vector<double> Vec12;
	//const int W = 50;
	//const double RADIUS = 20.0;
	int iLo = int(Pos) - m_NUDX_W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + m_NUDX_W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + 3 >= Pos && Pos2 <= Pos + 3)
			continue;
		double Dist = Chain.GetDist(Pos, Pos2);
		double DistFactor = exp(-Dist/m_NUDX_Radius);
		Chain.GetPt(Pos2, Pt2);
		Sub_Vecs(Pt2, PtCA, Vec12);
		double Theta = GetTheta_Vecs(VecPAB, Vec12);
		double Deg = degrees(Theta);
		if (Deg < 90)
			NU += DistFactor;
		else
			ND += DistFactor;
		}
	}

double DSS::Get_SSD2(uint Pos)
	{
	if (m_SS.empty())
		m_Chain->GetSS(m_SS);
	asserta(Pos < SIZE(m_SS));
	const uint L = GetSeqLength();
	//const int W = 100;
	//const uint w = 12;
	int iLo = int(Pos) - m_SSD2_W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + m_SSD2_W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	//for (uint Pos2 = 0; Pos2 < L; ++Pos2)
	double MinDist = 999;
	uint MinPos = UINT_MAX;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + m_SSD2_w >= Pos && Pos2 <= Pos + m_SSD2_w)
			continue;
		double Dist = m_Chain->GetDist(Pos, Pos2);
		if (Dist < MinDist)
			{
			MinDist = Dist;
			MinPos = Pos2;
			}
		}
	return MinDist;
	}

char DSS::Get_SSX2(uint Pos)
	{
	if (m_SS.empty())
		m_Chain->GetSS(m_SS);
	asserta(Pos < SIZE(m_SS));
	const uint L = GetSeqLength();
	//const int W = 100;
	//const uint w = 12;
	int iLo = int(Pos) - m_SSD2_W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + m_SSD2_W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	//for (uint Pos2 = 0; Pos2 < L; ++Pos2)
	double MinDist = 999;
	uint MinPos = UINT_MAX;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + m_SSD2_w >= Pos && Pos2 <= Pos + m_SSD2_w)
			continue;
		double Dist = m_Chain->GetDist(Pos, Pos2);
		if (Dist < MinDist)
			{
			MinDist = Dist;
			MinPos = Pos2;
			}
		}
	if (MinPos == UINT_MAX)
		return '~';
	char c = m_SS[MinPos];
	return c;
	}

const char *DSS::GetFeatureName(uint FeatureIndex)
	{
	switch (FeatureIndex)
		{
	case 0: return "NUDX";
	case 1: return "SSX";
	case 2: return "SSX2";
	case 3: return "SSD2";
	case 4: return "CMAA";
	case DSSFeatureIndex: return "DSS";
		}
	asserta(false);
	return "??";
	}

uint DSS::GetAlphaSize(uint FeatureIndex)
	{
	switch (FeatureIndex)
		{
	case 0: return m_NUDX_Bins;
	case 1: return 4;
	case 2: return 4;
	case 3: return 4;
	case 4: return 8;
	case DSSFeatureIndex:
		{
		uint n0 = GetAlphaSize(0);
		uint n1 = GetAlphaSize(1);
		uint n2 = GetAlphaSize(2);
		uint n3 = GetAlphaSize(3);
		uint n4 = GetAlphaSize(4);
		return n0*n1*n2*n3*n4;
		}
		}
	asserta(false);
	return UINT_MAX;
	}

void DSS::ToStr(string &s)
	{
	s.clear();
	const uint L = GetSeqLength();
	const string &Seq = m_Chain->m_Seq;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		char aa = toupper(Seq[Pos]);
		s += aa;
		uint i0 = GetFeature(0, Pos);
		uint i1 = GetFeature(1, Pos);
		uint i2 = GetFeature(2, Pos);
		uint i3 = GetFeature(3, Pos);
		uint Idx = GetIdx_NoCMAA(i0, i1, i2, i3);
		char Tmp[4];
		sprintf(Tmp, "%03x", Idx);
		s += Tmp[0];
		s += Tmp[1];
		s += Tmp[2];
		s += ' ';
		}
	}

uint DSS::GetFeature(uint FeatureIndex, uint Pos)
	{
	switch (FeatureIndex)
		{
	case 0:
		{
		double Value = Get_NUDX(Pos);
		if (m_NUDX_Bins == 8)
			{
			if (Value < 0.263661) return 0;
			if (Value < 0.404948) return 1;
			if (Value < 0.511572) return 2;
			if (Value < 0.600825) return 3;
			if (Value < 0.681155) return 4;
			if (Value < 0.759807) return 5;
			if (Value < 0.848306) return 6;
			return 7;
			}
		else if (m_NUDX_Bins == 4)
			{
			if (Value < 0.404948) return 0;
			if (Value < 0.600825) return 1;
			if (Value < 0.759807) return 2;
			return 3;
			}
		else
			asserta(false);
		}

	case 1:
		{
		char c = Get_SSX(Pos);
		switch (c)
			{
		case 'h': return 0;
		case 's': return 1;
		case 't': return 2;
		case '~': return 3;
			}
		asserta(false);
		return UINT_MAX;
		}

	case 2:
		{
		char c = Get_SSX2(Pos);
		switch (c)
			{
		case 'h': return 0;
		case 's': return 1;
		case 't': return 2;
		case '~': return 3;
			}
		asserta(false);
		return UINT_MAX;
		}

	case 3:
		{
		double d = Get_SSD2(Pos);
		if (d < m_SSD2_Step1)
			return 0;
		else if (d < m_SSD2_Step2)
			return 1;
		else if (d < m_SSD2_Step3)
			return 2;
		else
			return 3;
		}

	case 4:
		{
		byte c = m_Chain->m_Seq[Pos];
		uint Letter = g_CharToLetterAmino[c];
		if (Letter < 20)
			return m_AminoLetterToCompressedLetter[Letter];
		return UINT_MAX;
		}

	case DSSFeatureIndex:
		{
		uint i0 = GetFeature(0, Pos);
		uint i1 = GetFeature(1, Pos);
		uint i2 = GetFeature(2, Pos);
		uint i3 = GetFeature(3, Pos);
		uint i4 = GetFeature(4, Pos);

		uint n0 = GetAlphaSize(0);
		uint n1 = GetAlphaSize(1);
		uint n2 = GetAlphaSize(2);
		uint n3 = GetAlphaSize(3);

		uint i = i0 + i1*n0 + i2*n0*n1 + i3*n0*n1*n2 + i4*n0*n1*n2*n3;
		return i;
		}

	default:
		break;
		}
	asserta(false);
	return UINT_MAX;
	}

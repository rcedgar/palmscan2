#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "dss.h"

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
		}
	asserta(false);
	return "??";
	}

uint DSS::GetFeature(uint FeatureIndex, uint Pos)
	{
	switch (FeatureIndex)
		{
	case 0:
		{
		const uint BINS = 8;
		double Value = Get_NUDX(Pos);
		if (Value == DBL_MAX)
			return BINS/2;
		uint Bin = min(uint(Value*(BINS+2)), BINS-1);
		return Bin;
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

	default:
		break;
		}
	asserta(false);
	return UINT_MAX;
	}

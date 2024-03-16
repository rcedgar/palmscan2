#include "myutils.h"
#include "chainreader.h"
#include "outputfiles.h"
#include "abcxyz.h"
#include "pdbchain.h"
#include "xprof.h"
#include "quarts.h"

void GetHSE(const PDBChain &Chain, uint Pos, double Radius,
  uint MinPos, uint MaxPos, uint &NU, uint &ND);

void cmd_xprof()
	{
	const string &InputFN = opt_xprof;

	ChainReader CR;
	CR.Open(InputFN);

	XProf XP;
	PDBChain Chain;
	uint Counter = 0;
	while (CR.GetNext(Chain))
		{
		if (++Counter%100 == 0)
			Progress("%s\r", Chain.m_Label.c_str());
		XP.Init(Chain);
		XP.ToCfv(g_fcfv);
		}
	Progress("\n");
	}

void XProf::Init(const PDBChain &Chain)
	{
	m_Chain = &Chain;
	m_L = m_Chain->GetSeqLength();
	m_NUDX_ScaledValues.clear();
	m_SS.clear();
	}

void XProf::PosToCfv(FILE *f, uint Pos) const
	{
	if (f == 0)
		return;

	fprintf(f, "%u\t%c", Pos+1, m_Chain->m_Seq[Pos]);

	double x = m_Chain->m_Xs[Pos];
	double y = m_Chain->m_Ys[Pos];
	double z = m_Chain->m_Zs[Pos];

	fprintf(f, "\t%.1f\t%.1f\t%.1f", x, y, z);

	const uint FeatureCount = GetFeatureCount();
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount;
	  ++FeatureIndex)
		{
		double Value;
		uint iValue;
		GetFeature(FeatureIndex, Pos, Value, iValue);
		if (Value == DBL_MAX)
			fprintf(f, "\t.");
		else
			fprintf(f, "\t%.3g", Value);
		if (iValue == UINT_MAX)
			fprintf(f, "\t.");
		else
			fprintf(f, "\t%u", iValue);
		}
	fprintf(f, "\n");
	}

void XProf::ToCfv(FILE *f) const
	{
	if (f == 0)
		return;

	const uint L = GetSeqLength();
	const char *Label = m_Chain->m_Label.c_str();
	const uint FeatureCount = GetFeatureCount();

	fprintf(f, ">\t%s", Label);
	fprintf(f, "\t%u", L);
	fprintf(f, "\t%u", FeatureCount);
	for (uint i = 0; i < FeatureCount; ++i)
		{
		const char *Name = GetFeatureName(i);
		fprintf(f, "\t%s", Name);
		}
	fprintf(f, "\n");

	for (uint Pos = 0; Pos < L; ++Pos)
		PosToCfv(f, Pos);
	}

uint XProf::GetSphereNr(uint Pos, double Radius) const
	{
	vector<uint> PosVec;
	m_Chain->GetSphere(Pos, Radius, 0, m_L-1, PosVec);
	uint n = SIZE(PosVec);
	return n;
	}

void XProf::GetFeatures(uint Pos, vector<double> &Values) const
	{
	Values.clear();
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		{
		double Value = GetFeatureValue(FeatureIndex, Pos);
		Values.push_back(Value);
		}
	}

double XProf::GetFeatureValue(uint FeatureIndex, uint Pos) const
	{
	uint notused_iValue;
	double Value;
	GetFeature(FeatureIndex, Pos, Value, notused_iValue);
	return Value;
	}

void XProf::GetFeature(uint FeatureIndex, uint Pos,
  double &Value, uint &iValue) const
	{
	asserta(FeatureIndex < XFEATS);
	Value = DBL_MAX;
	iValue = UINT_MAX;
	switch (FeatureIndex)
		{
	case 0: Get_Ang_m2_p2(Pos, Value, iValue); return;
	case 1: Get_Ang_m3_p3(Pos, Value, iValue); return;
	case 2: Get_ED_p4(Pos, Value, iValue); return;
	case 3: Get_ED_m4(Pos, Value, iValue); return;
	case 4: Get_NU(Pos, Value, iValue); return;
	case 5: Get_ND(Pos, Value, iValue); return;
		}
	asserta(false);
	}

uint XProf::GetFeatureCount()
	{
	return 6;
	}

uint XProf::GetFeatureIndex(const string &FeatureName)
	{
	if (FeatureName == "Ang_m2_p2") return 0;
	if (FeatureName == "Ang_m3_p3") return 1;
	if (FeatureName == "ED_p4")		return 2;
	if (FeatureName == "ED_m4")		return 3;
	if (FeatureName == "NU")		return 4;
	if (FeatureName == "ND")		return 5;
	asserta(false);
	return UINT_MAX;
	}

const char *XProf::GetFeatureName(uint FeatureIndex)
	{
	switch (FeatureIndex)
		{
	case 0: return "Ang_m2_p2";
	case 1: return "Ang_m3_p3";
	case 2: return "ED_p4";
	case 3: return "ED_m4";
	case 4: return "NU";
	case 5: return "ND";
		}
	asserta(false);
	return "???";
	}

double XProf::GetAngle(uint PosA1, uint PosA2, 
  uint PosB1, uint PosB2) const
	{
	vector<double> PtA1;
	vector<double> PtA2;
	vector<double> PtB1;
	vector<double> PtB2;

	m_Chain->GetPt(PosA1, PtA1);
	m_Chain->GetPt(PosA2, PtA2);

	m_Chain->GetPt(PosB1, PtB1);
	m_Chain->GetPt(PosB2, PtB2);

	vector<double> vA;
	vector<double> vB;
	Sub_Vecs(PtA2, PtA1, vA);
	Sub_Vecs(PtB2, PtB1, vB);

	double Radians = GetTheta_Vecs(vA, vB);
	double Deg = degrees(Radians);
	return Deg;
	}

uint XProf::GetFeatureX(uint FeatureIndex, uint Pos)
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
		if (d < 5)
			return 0;
		else if (d < 7)
			return 1;
		else if (d < 9)
			return 2;
		else
			return 3;
		}

	case 4:
		{
		double Angle = Get_SSAngle2(Pos);
		uint i = uint((Angle*90)/360);
		i = i%4;
		return i;
		}

	default:
		break;
		}
	asserta(false);
	return UINT_MAX;
	}

char XProf::Get_SSX2(uint Pos)
	{
	if (m_SS.empty())
		m_Chain->GetSS(m_SS);
	asserta(Pos < SIZE(m_SS));
	const uint L = GetSeqLength();
	const int W = 100;
	const uint w = 12;
	int iLo = int(Pos) - W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	//for (uint Pos2 = 0; Pos2 < L; ++Pos2)
	double MinDist = 999;
	uint MinPos = UINT_MAX;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + w >= Pos && Pos2 <= Pos + w)
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

double XProf::Get_SSAngle2(uint Pos)
	{
	const uint L = GetSeqLength();
	const int W = 100;
	const uint w = 12;
	int iLo = int(Pos) - W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	//for (uint Pos2 = 0; Pos2 < L; ++Pos2)
	double MinDist = 999;
	uint MinPos = UINT_MAX;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + w >= Pos && Pos2 <= Pos + w)
			continue;
		double Dist = m_Chain->GetDist(Pos, Pos2);
		if (Dist < MinDist)
			{
			MinDist = Dist;
			MinPos = Pos2;
			}
		}
	if (MinPos == UINT_MAX)
		return 0;

	int PosA1 = int(Pos) - 1;
	int PosA2 = int(Pos) + 1;
	int PosB1 = int(MinPos) - 1;
	int PosB2 = int(MinPos) + 1;
	if (PosA1 < 0 || PosA2 >= int(L))
		return 0;
	if (PosB1 < 0 || PosB2 >= int(L))
		return 0;
	double Angle = GetAngle(PosA1, PosA2, PosB1, PosB2);
	return Angle;
	}

double XProf::Get_SSD2(uint Pos)
	{
	if (m_SS.empty())
		m_Chain->GetSS(m_SS);
	asserta(Pos < SIZE(m_SS));
	const uint L = GetSeqLength();
	const int W = 100;
	const uint w = 12;
	int iLo = int(Pos) - W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	//for (uint Pos2 = 0; Pos2 < L; ++Pos2)
	double MinDist = 999;
	uint MinPos = UINT_MAX;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + w >= Pos && Pos2 <= Pos + w)
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

char XProf::Get_SSX(uint Pos)
	{
	if (m_SS.empty())
		m_Chain->GetSS(m_SS);
	asserta(Pos < SIZE(m_SS));
	return m_SS[Pos];
	}

double XProf::Get_NUDX(uint Pos)
	{
	if (m_NUDX_ScaledValues.empty())
		Set_NUDXVec();
	asserta(Pos < SIZE(m_NUDX_ScaledValues));
	return m_NUDX_ScaledValues[Pos];
	}

void XProf::Set_NUDXVec()
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

void XProf::Get_NUDX_Lo(uint Pos, double &NU, double &ND) const
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
	const int W = 50;
	//const double RADIUS = 12.0;
	const double RADIUS = 20.0;
	int iLo = int(Pos) - W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	//for (uint Pos2 = 0; Pos2 < L; ++Pos2)
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + 3 >= Pos && Pos2 <= Pos + 3)
			continue;
		double Dist = Chain.GetDist(Pos, Pos2);
		double DistFactor = exp(-Dist/RADIUS);
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

void XProf::Get_NU(uint Pos, double &Value, uint &iValue) const
	{
	uint NU, ND;
	const uint L = GetSeqLength();
	GetHSE(*m_Chain, Pos, 12.0, 0, L-1, NU, ND);
	Value = double(NU);
	iValue = NU/3;
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_ND(uint Pos, double &Value, uint &iValue) const
	{
	uint NU, ND;
	const uint L = GetSeqLength();
	GetHSE(*m_Chain, Pos, 12.0, 0, L-1, NU, ND);
	Value = double(ND);
	iValue = ND;
	if (iValue > 10)
		iValue -= 10;
	iValue /= 2;
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_Ang_m1_p1(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos == 0 || Pos + 1 >= m_L)
		{
		Value = 0;
		iValue = 0;
		return;
		}
	Value = GetAngle(Pos-1, Pos, Pos, Pos+1);
	iValue = uint(Value/20.0);
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_Ang_m2_p2(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos < 2 || Pos + 2 >= m_L)
		{
		Value = 0;
		iValue = 0;
		return;
		}
	Value = GetAngle(Pos-2, Pos, Pos, Pos+2);
	iValue = IntScale(Value, 1.0, 120.0);
	}

void XProf::Get_Ang_m3_p3(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos < 3 || Pos + 3 >= m_L)
		{
		Value = 0;
		iValue = 0;
		return;
		}
	Value = GetAngle(Pos-3, Pos, Pos, Pos+3);
	iValue = IntScale(Value, 1.0, 90.0);
	}

void XProf::Get_Ang01_23(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos + 3 >= m_L)
		{
		Value = 0;
		iValue = 0;
		return;
		}
	Value = GetAngle(Pos, Pos+1, Pos+2, Pos+3);
	iValue = uint(Value*4.0/80.0);
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_Ang01_34(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos + 4 >= m_L)
		{
		Value = 0;
		iValue = 0;
		return;
		}
	Value = GetAngle(Pos, Pos+1, Pos+3, Pos+4);
	iValue = uint(Value*4.0/80.0);
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_Ang01_45(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos + 5 >= m_L)
		{
		Value = 0;
		iValue = 0;
		return;
		}
	Value = GetAngle(Pos, Pos+1, Pos+4, Pos+5);
	iValue = uint(Value*4.0/80.0);
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_ED_p4(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos + 4 >= m_L)
		{
		Value = 0;
		iValue = 0;
		return;
		}
	Value = m_Chain->GetDist(Pos, Pos+4);
	iValue = IntScale(Value, 4.5, 10.0);
	}

void XProf::Get_ED_m4(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos < 4)
		{
		Value = 0;
		iValue = 0;
		return;
		}
	Value = m_Chain->GetDist(Pos-4, Pos);
	iValue = IntScale(Value, 4.5, 10.0);
	}

uint XProf::IntScale(double Value, double MinVal, double HiQ) const
	{
	double x = Value;
	if (x < MinVal)
		x = MinVal;
	x -= MinVal;
	uint i = (uint) ((x*XBINS)/(HiQ));
	if (i >= XBINS)
		i = XBINS - 1;
	return i;
	}

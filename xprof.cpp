#include "myutils.h"
#include "chainreader.h"
#include "outputfiles.h"
#include "abcxyz.h"
#include "xprof.h"

void cmd_xprof()
	{
	const string &InputFN = opt_xprof;

	ChainReader CR;
	CR.Open(InputFN);

	XProf XP;
	PDBChain Chain;
	while (CR.GetNext(Chain))
		{
		XP.Init(Chain);
		XP.ToTsv(g_ftsv);
		}
	}

void XProf::Init(const PDBChain &Chain)
	{
	m_Chain = &Chain;
	m_L = m_Chain->GetSeqLength();
	}

void XProf::PosToTsv(FILE *f, uint Pos) const
	{
	fprintf(f, "%s\t%u", m_Chain->m_Label.c_str(), Pos+1);
	const uint FeatureCount = GetFeatureCount();
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount;
	  ++FeatureIndex)
		{
		double Value = GetFeature(Pos, FeatureIndex);
		if (Value == DBL_MAX)
			fprintf(f, "\t0");
		else
			fprintf(f, "\t%.3g", Value);
		}
	fprintf(f, "\n");
	}

void XProf::ToTsv(FILE *f) const
	{
	const uint L = GetSeqLength();
	for (uint Pos = 0; Pos < L; ++Pos)
		PosToTsv(f, Pos);
	}

double XProf::GetFeature0(uint Pos) const
	{
	if (Pos + 3 >= m_L)
		return DBL_MAX;
	double Deg = GetAngle(Pos, Pos+1, Pos+2, Pos+3);
	return Deg;
	}

double XProf::GetFeature1(uint Pos) const
	{
	if (Pos + 4 >= m_L)
		return DBL_MAX;
	double Deg = GetAngle(Pos, Pos+1, Pos+3, Pos+4);
	return Deg;
	}

double XProf::GetFeature2(uint Pos) const
	{
	if (Pos + 5 >= m_L)
		return DBL_MAX;
	double Deg = GetAngle(Pos, Pos+1, Pos+4, Pos+5);
	return Deg;
	}

double XProf::GetFeature3(uint Pos) const
	{
	uint n = GetSphereNr(Pos, 6.0);
	return double(n);
	}

double XProf::GetFeature4(uint Pos) const
	{
	uint n = GetSphereNr(Pos, 12.0);
	return double(n);
	}

uint XProf::GetSphereNr(uint Pos, double Radius) const
	{
	vector<uint> PosVec;
	m_Chain->GetSphere(Pos, Radius, 0, m_L-1, PosVec);
	uint n = SIZE(PosVec);
	return n;
	}

uint XProf::GetFeatureCount() const
	{
	return 5;
	}

double XProf::GetFeature(uint Pos, uint FeatureIndex) const
	{
	switch (FeatureIndex)
		{
	case 0: return GetFeature0(Pos);
	case 1: return GetFeature1(Pos);
	case 2: return GetFeature2(Pos);
	case 3: return GetFeature3(Pos);
	case 4: return GetFeature4(Pos);
		}
	asserta(false);
	return DBL_MAX;
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

#include "myutils.h"
#include "structprof.h"
#include "chainreader.h"
#include "outputfiles.h"
#include "cmpsearcher.h"
#include "abcxyz.h"

void StructProf::SetChain(const PDBChain &Chain)
	{
	Clear();
	m_Chain = &Chain;
	m_MinPos = 0;
	uint L = Chain.GetSeqLength();
	asserta(L > 0);
	m_MaxPos = L - 1;
	}

void StructProf::SetMinMaxPos(uint MinPos, uint MaxPos)
	{
	asserta(MinPos <= MaxPos && MaxPos < m_Chain->GetSeqLength());
	m_MinPos = MinPos;
	m_MaxPos = MaxPos;
	}

uint StructProf::SearchDist(uint Pos, uint Lo, uint Hi,
  bool Maximize, double X, double &BestDist) const
	{
	const PDBChain &Chain = *m_Chain;
	const uint L = Chain.GetSeqLength();
	asserta(Lo < Hi && Hi < L);
	BestDist = (Maximize ? 0 : DBL_MAX);
	uint BestPos = UINT_MAX;
	for (uint Pos2 = Lo; Pos2 <= Hi; ++Pos2)
		{
		double Dist = Chain.GetDist(Pos2, Pos);
		if (Maximize)
			{
			if (Dist > BestDist)
				{
				BestDist = Dist;
				BestPos = Pos2;
				}
			else if (Dist < BestDist - X)
				break;
			}
		else
			{
			if (Dist < BestDist)
				{
				BestDist = Dist;
				BestPos = Pos2;
				}
			else if (Dist > BestDist + X)
				break;
			}
		}
	return BestPos;
	}

void StructProf::GetSphere(const vector<double> &CenterPt,
  double Radius, vector<uint> &PosVec) const
	{
	PosVec.clear();
	const uint L = m_Chain->GetSeqLength();
	asserta(m_MinPos <= m_MaxPos && m_MaxPos < L);
	vector<double> Pt;
	for (uint Pos2 = m_MinPos; Pos2 <= m_MaxPos; ++Pos2)
		{
		m_Chain->GetPt(Pos2, Pt);
		double d = GetDist(CenterPt, Pt);
		if (d <= Radius)
			PosVec.push_back(Pos2);
		}
	}

void StructProf::GetHSE(uint Pos, double Radius,
  uint &NU, uint &ND) const
	{
	NU = 0;
	ND = 0;
	const uint L = m_Chain->GetSeqLength();
	if (Pos == 0 || Pos+1 >= L)
		return;

	vector<double> PtPrevCA;
	vector<double> PtNextCA;
	vector<double> PtCA;
	vector<double> PtCB;
	m_Chain->GetPt(Pos-1, PtPrevCA);
	m_Chain->GetPt(Pos, PtCA);
	m_Chain->GetPt(Pos+1, PtNextCA);

	vector<double> d1;
	vector<double> d2;
	Sub_Vecs(PtCA, PtPrevCA, d1);
	Sub_Vecs(PtCA, PtNextCA, d2);

	vector<double> VecPAB;
	Add_Vecs(d1, d2, VecPAB);
	NormalizeVec(VecPAB);

	vector<uint> SpherePosVec;
	GetSphere(PtCA, Radius, SpherePosVec);

	const uint N = SIZE(SpherePosVec);
	vector<double> Pt2;
	vector<double> Vec12;
	for (uint i = 0; i < N; ++i)
		{
		uint Pos2 = SpherePosVec[i];
		if (Pos2 == Pos)
			continue;
		m_Chain->GetPt(Pos2, Pt2);
		Sub_Vecs(Pt2, PtCA, Vec12);
		double Theta = GetTheta_Vecs(VecPAB, Vec12);
		double Deg = degrees(Theta);
		if (Deg < 90)
			++NU;
		else
			++ND;
		}
	}

uint StructProf::GetTSB(uint Pos, double Radius) const
	{
	const uint L = m_Chain->GetSeqLength();
	if (Pos == 0 || Pos+1 >= L)
		return 0;

	vector<double> PtPrevCA;
	vector<double> PtNextCA;
	vector<double> PtCA;
	vector<double> PtCB;
	m_Chain->GetPt(Pos-1, PtPrevCA);
	m_Chain->GetPt(Pos, PtCA);
	m_Chain->GetPt(Pos+1, PtNextCA);

	vector<double> d1;
	vector<double> d2;
	Sub_Vecs(PtCA, PtPrevCA, d1);
	Sub_Vecs(PtCA, PtNextCA, d2);

	vector<double> VecPAB;
	Add_Vecs(d1, d2, VecPAB);
	NormalizeVec(VecPAB);

	vector<double> CenterPt;
	for (uint i = 0; i < 3; ++i)
		{
		double Coord = PtCA[i] + Radius*VecPAB[i];
		CenterPt.push_back(Coord);
		}

	vector<uint> SpherePosVec;
	GetSphere(CenterPt, Radius, SpherePosVec);

	const uint N = SIZE(SpherePosVec);
	return N;
	}

static void DoStructProfPos(FILE *f, const StructProf &SP, uint Pos)
	{
	if (f == 0)
		return;

	const PDBChain &Chain = *SP.m_Chain;
	const char *Label = Chain.m_Label.c_str();
	const char *Seq = Chain.m_Seq.c_str();
	char aa = Seq[Pos];
	char ss = Chain.m_SS[Pos];
	uint NU, ND;
	SP.GetHSE(Pos, 12.0, NU, ND);
	uint TSB = SP.GetTSB(Pos, 10.0);

	uint Pos_aD = Chain.GetMotifPos(A) + 3;
	uint Pos_bG = Chain.GetMotifPos(B) + 1;
	uint Pos_cD = Chain.GetMotifPos(C) + 3;
	asserta(Seq[Pos_aD] == 'D');
	asserta(Seq[Pos_bG] == 'G');
	asserta(Seq[Pos_cD] == 'D');

	double Dist_aD = Chain.GetDist(Pos, Pos_aD);
	double Dist_bG = Chain.GetDist(Pos, Pos_bG);
	double Dist_cD = Chain.GetDist(Pos, Pos_cD);

	static bool HdrDone = false;
	if (!HdrDone)
		{
		fprintf(f, "Label");
		fprintf(f, "\tPos");
		fprintf(f, "\taa");
		fprintf(f, "\tss");
		fprintf(f, "\tNU");
		fprintf(f, "\tND");
		fprintf(f, "\tTSB");
		fprintf(f, "\tDist_aD");
		fprintf(f, "\tDist_bG");
		fprintf(f, "\tDist_cD");
		fprintf(f, "\n");
		HdrDone = true;
		}

	fprintf(f, "%s", Label);
	fprintf(f, "\t%u", Pos+1);
	fprintf(f, "\t%c", aa);
	fprintf(f, "\t%c", ss);
	fprintf(f, "\t%u", NU);
	fprintf(f, "\t%u", ND);
	fprintf(f, "\t%u", TSB);
	fprintf(f, "\t%.1f", Dist_aD);
	fprintf(f, "\t%.1f", Dist_bG);
	fprintf(f, "\t%.1f", Dist_cD);
	fprintf(f, "\n");
	}

static void DoStructProf(FILE *f, CMPSearcher &CS,
 PDBChain &Chain)
	{
	if (f == 0)
		return;

	CS.Search(Chain);

	uint APos = UINT_MAX;
	uint BPos = UINT_MAX;
	uint CPos = UINT_MAX;
	double PalmScore = CS.GetPSSMStarts(APos, BPos, CPos);
	if (PalmScore <= 0)
		return;

	Chain.SetMotifPosVec(APos, BPos, CPos);

	Chain.PrintSeqCoords(g_fLog);
	StructProf SP;
	SP.SetChain(Chain);
	const uint L = Chain.GetSeqLength();
	for (uint Pos = 0; Pos < L; ++Pos)
		DoStructProfPos(f, SP, Pos);

	uint PosD = SP.FindMofifD_Hueuristics();
	uint PosE = SP.FindMofifE_Hueuristics(PosD);
	}

void cmd_struct_prof()
	{
	const string &InputFN = opt_struct_prof;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;

	CMP Prof;
	Prof.FromFile(ModelFileName);

	CMPSearcher CS;
	CS.SetProf(Prof);

	ChainReader CR;
	CR.Open(InputFN, false);

	PDBChain Chain;
	while (CR.GetNext(Chain))
		{
		Chain.SetSS();
		DoStructProf(g_ftsv, CS, Chain);
		}
	}

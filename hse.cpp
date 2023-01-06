#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include <map>

static void LogVec(const string &Name, const vector<double> &v)
	{
	asserta(SIZE(v) == 3);
	Log("  x=%8.3f,  y=%8.3f,  z=%8.3f  %s\n", v[X], v[Y], v[Z], Name.c_str());
	}

void ReadMotifCoords(
  vector<vector<uint> > &MotifCoordsVec,
  vector<string> &Labels,
  map<string, uint> &LabelToIndex);

void GetHSE(const PDBChain &Chain, uint Pos, double Radius,
  uint &NU, uint &ND)
	{
	NU = 0;
	ND = 0;
	if (Pos == 0 || Pos+1 >= SIZE(Chain.m_Seq))
		return;

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

	vector<uint> SpherePosVec;
	Chain.GetSphere(Pos, Radius, SpherePosVec);

	const uint N = SIZE(SpherePosVec);
	vector<double> Pt2;
	vector<double> Vec12;
	for (uint i = 0; i < N; ++i)
		{
		uint Pos2 = SpherePosVec[i];
		Chain.GetPt(Pos2, Pt2);
		Sub_Vecs(Pt2, PtCA, Vec12);
		double Theta = GetTheta_Vecs(VecPAB, Vec12);
		double Deg = degrees(Theta);
		if (Deg < 90)
			++NU;
		else
			++ND;
		}
	}

static void HSE(const PDBChain &Chain, double Radius)
	{
	const uint L = Chain.GetSeqLength();
	for (uint Pos = 1; Pos + 1 < L; ++Pos)
		{
		char c = Chain.m_Seq[Pos];

		uint NU, ND;
		GetHSE(Chain, Pos, Radius, NU, ND);
		if (g_ftsv != 0)
			{
			fprintf(g_ftsv, "%s", Chain.m_Label.c_str());
			fprintf(g_ftsv, "\t%u", Pos+1);
			fprintf(g_ftsv, "\t%c", c);
			fprintf(g_ftsv, "\t%u", NU);
			fprintf(g_ftsv, "\t%u", ND);
			fprintf(g_ftsv, "\n");
			}
		}
	}

static void HSETrain(const PDBChain &Chain, double Radius,
  uint APos, uint BPos, uint CPos, vector<uint> &NUs)
	{
	NUs.clear();

	for (uint Pos = APos; Pos < APos + AL; ++Pos)
		{
		uint NU, ND;
		GetHSE(Chain, Pos, Radius, NU, ND);
		NUs.push_back(NU);
		}

	for (uint Pos = BPos; Pos < BPos + BL; ++Pos)
		{
		uint NU, ND;
		GetHSE(Chain, Pos, Radius, NU, ND);
		NUs.push_back(NU);
		}

	for (uint Pos = CPos; Pos < CPos + CL; ++Pos)
		{
		uint NU, ND;
		GetHSE(Chain, Pos, Radius, NU, ND);
		NUs.push_back(NU);
		}
	}

void cmd_hse()
	{
	const string &FN = opt_hse;

	vector<PDBChain *> Chains;
	ReadChains(FN, Chains, true);
	uint FlankSize = 0;

	const uint ChainCount = SIZE(Chains);
	vector<double> Angles;
	for (uint i = 0; i < ChainCount; ++i)
		{
		ProgressStep(i, ChainCount, "Processing");
		PDBChain &Chain = *Chains[i];
		HSE(Chain, 12.0);
		}
	}

void cmd_hse_train()
	{
	const string &QueryFN = opt_hse_train;

	vector<vector<uint> > MotifCoordsVec;
	vector<string> Labels;
	map<string, uint> LabelToIndex;
	ReadMotifCoords(MotifCoordsVec, Labels, LabelToIndex);

	PDBChain Chain;

	ChainReader CR;
	CR.Open(QueryFN, false);
	double Radius = 12.0;
	if (optset_radius)
		Radius = opt_radius;
	uint NotFoundCount = 0;
	vector<vector<uint> > NUVec;
	vector<string> FoundLabels;
	const uint ML = AL + BL + CL;
	for (;;)
		{
		bool Ok = CR.GetNext(Chain);
		if (!Ok)
			break;

		const string &Label = Chain.m_Label;
		map<string, uint>::const_iterator p = LabelToIndex.find(Label);
		if (p == LabelToIndex.end())
			{
			if (NotFoundCount < 10)
				Log("Not found >%s\n", Label.c_str());
			++NotFoundCount;
			continue;
			}
		uint Index = p->second;

		asserta(Index < SIZE(MotifCoordsVec));
		const vector<uint> &MotifCoords = MotifCoordsVec[Index];
		asserta(SIZE(MotifCoords) == 3);

		uint APos = MotifCoords[0];
		uint BPos = MotifCoords[1];
		uint CPos = MotifCoords[2];

		Log("%u %u %u >%s\n", APos, BPos, CPos, Label.c_str());

		vector<uint> NUs;
		HSETrain(Chain, Radius, APos, BPos, CPos, NUs);
		asserta(SIZE(NUs) == ML);
		FoundLabels.push_back(Label);
		NUVec.push_back(NUs);
		}

	if (g_ftsv == 0)
		return;

	const uint N = SIZE(FoundLabels);
	asserta(SIZE(NUVec) == N);
	fprintf(g_ftsv, "Pos");
	for (uint i = 0; i < N; ++i)
		fprintf(g_ftsv, "\t%s", FoundLabels[i].c_str());
	fprintf(g_ftsv, "\n");

	for (uint Pos = 0; Pos < ML; ++Pos)
		{
		fprintf(g_ftsv, "%u", Pos);
		for (uint i = 0; i < N; ++i)
			fprintf(g_ftsv, "\t%u", NUVec[i][Pos]);
		fprintf(g_ftsv, "\n");
		}
	}

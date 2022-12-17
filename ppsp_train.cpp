#include "myutils.h"
#include "pdbchain.h"
#include "rdrpmodel.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "ppsp.h"

static uint g_DoneCount;
static uint g_HitCount;

//static uint GetSeqPos(uint i, uint APos, uint BPos, uint CPos)
//	{
//	if (i >= AL + BL)
//		return CPos + i - (AL + BL);
//	if (i >= AL)
//		return BPos + i - AL;
//	return APos + i;
//	}
//
//static bool GetDistMx(const PDBChain &Q, 
//  uint APos, uint BPos, uint CPos,
//  vector<vector<double> > &DistMx)
//	{
//	DistMx.clear();
//	if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
//		return false;
//
//	const uint N = AL + BL + CL;
//	DistMx.resize(N);
//	for (uint i = 0; i < N; ++i)
//		DistMx[i].resize(N, DBL_MAX);
//
//	for (uint i = 0; i < N; ++i)
//		{
//		uint SeqPosi = GetSeqPos(i, APos, BPos, CPos);
//		for (uint j = 0; j < N; ++j)
//			{
//			uint SeqPosj = GetSeqPos(j, APos, BPos, CPos);
//			double d = Q.GetDist(SeqPosi, SeqPosj);
//			DistMx[i][j] = d;
//			}
//		}
//	return true;
//	}

static void GetMeanStdDev(
 const vector<vector<vector<double> > > &DistMxVec,
  uint i, uint j, double &Mean, double &StdDev)
	{
	const uint N = SIZE(DistMxVec);
	double Sum = 0;
	for (uint k = 0; k < N; ++k)
		{
		asserta(i < SIZE(DistMxVec[k]));
		asserta(j < SIZE(DistMxVec[k][i]));
		double d = DistMxVec[k][i][j];
		Sum += d;
		}
	Mean = double(Sum)/N;

	double Sumd2 = 0;
	for (uint k = 0; k < N; ++k)
		{
		double d = DistMxVec[k][i][j];
		double d2 = (d - Mean)*(d - Mean);
		Sumd2 += d2;
		}
	StdDev = (double) sqrt(Sumd2/N);
	}

void cmd_ppsp_train()
	{
	const string &QueryFN = opt_ppsp_train;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;
	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

	PDBChain Q;
	RdRpSearcher RS;
	RS.Init(Model);

	ChainReader CR;
	CR.m_SaveAtoms = true;
	CR.Open(QueryFN, false);

	vector<vector<vector<double> > > DistMxVec;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

		if (++g_DoneCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}

		const string &QSeq = Q.m_Seq;
		const string &QLabel = Q.m_Label;

		RS.Search(QLabel, QSeq);
		RS.WriteOutput();

		uint APos = RS.GetMotifPos(0);
		uint BPos = RS.GetMotifPos(1);
		uint CPos = RS.GetMotifPos(2);

		vector<vector<double> > DistMx;
		Ok = PPSP::GetDistMx(Q, APos, BPos, CPos, DistMx);
		if (!Ok)
			continue;
		++g_HitCount;
		DistMxVec.push_back(DistMx);
		}

	PPSP Prof;

	for (uint i = 0; i < PPSPL; ++i)
		{
		Prof.m_Means[i][i] = 0;
		Prof.m_StdDevs[i][i] = 0;
		for (uint j = 0; j < i; ++j)
			{
			double Mean, StdDev;
			GetMeanStdDev(DistMxVec, i, j, Mean, StdDev);

			Prof.m_Means[i][j] = Mean;
			Prof.m_Means[j][i] = Mean;

			Prof.m_StdDevs[i][j] = StdDev;
			Prof.m_StdDevs[j][i] = StdDev;
			}
		}

	Progress("100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);

	Prof.ToFile(opt_output);
	}

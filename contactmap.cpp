#include "myutils.h"
#include "pdbchain.h"
#include "rdrpmodel.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"

static uint g_DoneCount;
static uint g_HitCount;

static uint GetSeqPos(uint i, uint APos, uint BPos, uint CPos)
	{
	if (i >= AL + BL)
		return CPos + i - (AL + BL);
	if (i >= AL)
		return BPos + i - AL;
	return APos + i;
	}

static void PPContactMap(const PDBChain &Q, 
  uint APos, uint BPos, uint CPos)
	{
	if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
		return;

	const uint N = AL + BL + CL;
	vector<vector<double> > Mx;
	Mx.resize(N);
	for (uint i = 0; i < N; ++i)
		Mx[i].resize(N, DBL_MAX);

	for (uint i = 0; i < N; ++i)
		{
		uint SeqPosi = GetSeqPos(i, APos, BPos, CPos);
		for (uint j = 0; j < N; ++j)
			{
			uint SeqPosj = GetSeqPos(j, APos, BPos, CPos);
			double d = Q.GetDist(SeqPosi, SeqPosj);
			Mx[i][j] = d;
			}
		}
#if 0
	{
	Log("\n");
	Log(">%s\n", Q.m_Label.c_str());
	for (uint i = 0; i < N; ++i)
		{
		for (uint j = 0; j <= i; ++j)
			Log(" %8.3g", Mx[i][j]);
		Log("\n");
		}
	}
#endif
	FILE *f = g_fppcontactmap_tsv;
	if (f != 0)
		{
		const char *Label = Q.m_Label.c_str();
		for (uint i = 0; i < N; ++i)
			{
			fprintf(f, "%s\t%u", Label, i);
			for (uint j = 0; j <= i; ++j)
				fprintf(f, "\t%.3g", Mx[i][j]);
			fprintf(f, "\n");
			}
		}
	}

static void Thread(ChainReader &CR, vector<PDBChain> &Qs,  
  vector<RdRpSearcher> &RSs)
	{
	uint ThreadIndex = GetThreadIndex();
	PDBChain &Q = Qs[ThreadIndex];
	RdRpSearcher &RS = RSs[ThreadIndex];
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			return;
#pragma omp critical
		{
		if (++g_DoneCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}
		}

		const string &QSeq = Q.m_Seq;
		const string &QLabel = Q.m_Label;

		RS.Search(QLabel, QSeq);
		RS.WriteOutput();

		uint APos = RS.GetMotifPos(0);
		uint BPos = RS.GetMotifPos(1);
		uint CPos = RS.GetMotifPos(2);

		PPContactMap(Q, APos, BPos, CPos);
		}
	}

void cmd_ppcontactmap()
	{
	const string &QueryFN = opt_ppcontactmap;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;
	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

	uint ThreadCount = GetRequestedThreadCount();

	vector<vector<PDBChain *> > QVecs(ThreadCount);
	
	uint ThreadFinishedCount = 0;
	uint DoneCount = 0;
	uint HitCount = 0;
	vector<PDBChain> Qs(ThreadCount);
	vector<bool> ThreadDone(ThreadCount);
	vector<RdRpSearcher> RSs(ThreadCount);
	for (uint i = 0; i < ThreadCount; ++i)
		RSs[i].Init(Model);

	ChainReader CR;
	CR.m_SaveAtoms = true;
	CR.Open(QueryFN, false);

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, Qs, RSs);

	Progress("100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);
	}

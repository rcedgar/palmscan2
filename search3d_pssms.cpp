#include "myutils.h"
#include "pdbchain.h"
#include "rdrpmodel.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"

static uint g_DoneCount;
static uint g_HitCount;

static bool SplitB(const PDBChain &Q, uint APos, uint BPos, uint CPos)
	{
	if (!optset_splitb)
		return false;
	if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
		return false;

	uint L = Q.GetSeqLength();
	uint Lo1 = 0;
	if (BPos > 400)
		Lo1 = BPos - 400;
	uint Hi1 = BPos + 14;

	uint Lo2 = BPos; 
	uint Hi2 = BPos + 400;
	if (Hi2 >= L)
		Hi2 = L - 1;

	PDBChain Split1;
	PDBChain Split2;
	Q.GetRange(Lo1, Hi1, Split1);
	Q.GetRange(Lo2, Hi2, Split2);

	string s;
	for (uint i = 0; i < SIZE(Q.m_Label); ++i)
		{
		char c = Q.m_Label[i];
		if (!isspace(c))
			s += c;
		}
	const char *Label = s.c_str();
	string FileName1;
	string FileName2;
	Ps(FileName1, "%s.split1.pdb", Label);
	Ps(FileName2, "%s.split2.pdb", Label);
	Split1.ToPDB(FileName1);
	Split2.ToPDB(FileName2);
	return true;
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

		if (optset_splitb)
#pragma omp critical
			{
			bool Ok = SplitB(Q, APos, BPos, CPos);
			if (Ok)
				exit(0);
			}

#pragma omp critical
		{
		if (RS.m_TopPalmHit.m_Score > 0)
			++g_HitCount;
		}
		}
	}

void cmd_search3d_pssms()
	{
	const string &QueryFN = opt_search3d_pssms;

	RdRpModel Model;
	GetRdrpModel(Model);

	//RdRpSearcher::InitOutput();
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
	CR.Open(QueryFN);

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, Qs, RSs);

	Progress("100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);
	}

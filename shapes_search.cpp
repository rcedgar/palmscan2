#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "shapesearcher.h"
#include "motifsettings.h"
#include <time.h>

#if 0
static uint g_DoneCount;
static uint g_HitCount;

static bool Search1(const PDBChain &Q, ShapeSearcher &SS)
	{
	SS.ClearSearch();
	SS.SetQuery(Q);
////	SS.Search();
	bool IsHit = (SS.m_Score >= 0.5);
	if (opt_misses || IsHit)
		SS.ToTsv(g_ftsv);
	if (IsHit)
		{
		string LoadName;
		if (optset_loaddir)
			LoadName = string(opt_loaddir) + "/";
		LoadName += Q.m_Label;
		LoadName += ".pdb";
		SS.ToPml(g_fpml, LoadName.c_str());
		}
	return IsHit;
	}

static void Thread(ChainReader &CR, const Shapes &S)
	{
	ShapeSearcher SS;
	SS.Init(S);

	PDBChain Q;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

#pragma omp critical
		{
		if (++g_DoneCount%1000 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}
		}

		const char *Seq = Q.m_Seq.c_str();
		const string &Label = Q.m_Label;

		vector<uint> PosVec;
		vector<string> ShapeSeqs;
		bool IsHit = Search1(Q, SS);
		if (IsHit)
			++g_HitCount;
		}
	}
#endif 

void cmd_shapes_search()
	{
	Die("shapes_search");
#if 0
	const string &QueryFN = opt_shapes_search;
	Shapes S;
	S.InitFromCmdLine();

	ChainReader CR;
	CR.Open(QueryFN);

	uint ThreadCount = GetRequestedThreadCount();

	time_t tStart = time(0);
#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, S);
	time_t tEnd = time(0);
	uint Secs = uint(tEnd - tStart);
	if (Secs <= 0)
		Secs = 1;
	double Throughput = double(g_DoneCount)/(Secs*ThreadCount);
	ProgressLog("%u/%u hits (%.3g%%), %s secs (%u threads, %.1f/sec/thread)\n",
	  g_HitCount, g_DoneCount, GetPct(g_HitCount, g_DoneCount),
	  IntToStr(Secs), ThreadCount, Throughput);
#endif
	}

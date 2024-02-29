#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "shapesearcher.h"
#include "motifsettings.h"
#include <time.h>

static uint g_DoneCount;
static uint g_HitCount;

void AItoms(const PDBChain &Q, uint PosA, uint PosB, uint PosC);

static bool Search1(const PDBChain &Q, ShapeSearcher &SS)
	{
	SS.SearchDom(Q);
	SS.SetFinalScore();
	SS.SetClass();
	bool Hit = SS.IsHit();
	if (opt_misses || Hit)
		{
		SS.ToTsv(g_ftsv);

		uint PosA = SS.GetShapePos(SS.m_ShapeIndexA);
		uint PosB = SS.GetShapePos(SS.m_ShapeIndexB);
		uint PosC = SS.GetShapePos(SS.m_ShapeIndexC);
		AItoms(Q, PosA, PosB, PosC);
		}
	return Hit;
	}

static void Thread(ChainReader &CR, const Shapes &S)
	{
	ShapeSearcher SS;
	SS.Init(S);
	if (GetThreadIndex() == 0)
		SS.LogParams();

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
#pragma omp critical
			{
			++g_HitCount;
			}
		}
	}

void cmd_aitoms()
	{
	const string &QueryFN = g_Arg1;
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
	}

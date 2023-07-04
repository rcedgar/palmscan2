#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "shapesearcher.h"
#include <time.h>

static uint g_DoneCount;
static uint g_HitCount;

void ShapeSearcher::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	static bool HdrDone = false;
#pragma omp critical
	{
	//Log("\n");
	//Log("# score %.3g\n", m_ScoreABC);
	//Log("fetch %s\n", m_Query->m_Label.c_str());
//	ToPmlABC(g_fLog);

	if (!HdrDone)
		{
		HdrDone = true;
		fprintf(f, "Label");
		fprintf(f, "\tpalm_score");
		fprintf(f, "\tpp_score");
		fprintf(f, "\tgate");
		for (uint i = 0; i < m_ShapeCount; ++i)
			{
			const char *ShapeName = GetShapeName(i);
			fprintf(f, "\t%s_pos", ShapeName);
			fprintf(f, "\t%s_seq", ShapeName);
			fprintf(f, "\t%s_score", ShapeName);
			}
		fprintf(f, "\n");
		}

	char Gate = GetGate();
	fprintf(f, "%s", m_Query->m_Label.c_str());
	fprintf(f, "\t%.3g", m_PalmScore);
	fprintf(f, "\t%.3g", m_ScoreABC);
	fprintf(f, "\t%c", Gate);
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		double Score = m_ShapeScores[i];
		string Seq;
		GetShapeSeq(i, Seq);
		uint Pos = m_ShapePosVec[i];
		if (Pos == UINT_MAX)
			fprintf(f, "\t.");
		else
			fprintf(f, "\t%u", Pos+1);
		fprintf(f, "\t%s", Seq.c_str());
		fprintf(f, "\t%.3g", Score);
		}
	fprintf(f, "\n");
	}
	}

static bool Search1(const PDBChain &Q, ShapeSearcher &SS)
	{
	SS.SearchPalm(Q);
	const uint ShapeCount = SS.GetShapeCount();
	asserta(SIZE(SS.m_ShapePosVec) == ShapeCount);
	asserta(SIZE(SS.m_ShapeScores) == ShapeCount);

	uint IXA = SS.m_ShapeIndexA;
	uint IXB = SS.m_ShapeIndexB;
	uint IXC = SS.m_ShapeIndexC;
	bool IsHit = (SS.m_ScoreABC >= SS.m_MinScoreABC);
	if (opt_misses || IsHit)
		SS.ToTsv(g_ftsv);
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

void cmd_shapes_search()
	{
	const string &QueryFN = opt_shapes_search;

	time_t tStart = time(0);
	if (!optset_shapes)
		Die("Must specify -shapes");
	const string &ShapesFileName = opt_shapes;

	Shapes S;
	S.FromFile(ShapesFileName);

	ChainReader CR;
	CR.Open(QueryFN);

	uint ThreadCount = GetRequestedThreadCount();

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, S);
	time_t tEnd = time(0);
	uint Secs = uint(tEnd - tStart);
	if (Secs <= 0)
		Secs = 1;
	double Throughput = double(g_DoneCount)/(Secs*ThreadCount);
	ProgressLog("%u/%u hits (%.3g%%), %s secs (%u threads, %.1f/ sec/ thread)\n",
	  g_HitCount, g_DoneCount, GetPct(g_HitCount, g_DoneCount),
	  IntToStr(Secs), ThreadCount, Throughput);
	}

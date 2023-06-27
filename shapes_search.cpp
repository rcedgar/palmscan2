#include "myutils.h"
#include "pdbchain.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
//#include "cmpsearcher.h"
#include "shapesearcher.h"
#include <time.h>

static uint g_DoneCount;
static uint g_HitCount;

static bool Search1(const PDBChain &Q, ShapeSearcher &SS)
	{
	SS.SearchPalm(Q);
	const uint ShapeCount = SS.GetShapeCount();
	asserta(SIZE(SS.m_ShapePosVec) == ShapeCount);
	asserta(SIZE(SS.m_ShapeScores) == ShapeCount);

	uint IXA = SS.m_ShapeIndexA;
	uint IXB = SS.m_ShapeIndexB;
	uint IXC = SS.m_ShapeIndexC;
	bool IsHit = (SS.m_ShapePosVec[IXA] != UINT_MAX &&
	  SS.m_ShapePosVec[IXB] != UINT_MAX
	  && SS.m_ShapePosVec[IXC] != UINT_MAX);

	if (g_ftsv != 0)
		{
		static bool HdrDone = false;
#pragma omp critical
			{
			if (!HdrDone)
				{
				HdrDone = true;
				fprintf(g_ftsv, "Label");
				for (uint i = 0; i < ShapeCount; ++i)
					{
					const char *ShapeName = SS.GetShapeName(i);
					fprintf(g_ftsv, "\t%s_pos", ShapeName);
					fprintf(g_ftsv, "\t%s_seq", ShapeName);
					fprintf(g_ftsv, "\t%s_score", ShapeName);
					}
				fprintf(g_ftsv, "\n");
				}

			fprintf(g_ftsv, "%s", Q.m_Label.c_str());
			for (uint i = 0; i < ShapeCount; ++i)
				{
				double Score = SS.m_ShapeScores[i];
				string Seq;
				SS.GetShapeSeq(i, Seq);
				uint Pos = SS.m_ShapePosVec[i];
				if (Pos == UINT_MAX)
					fprintf(g_ftsv, "\t.");
				else
					fprintf(g_ftsv, "\t%u", Pos);
				fprintf(g_ftsv, "\t%s", Seq.c_str());
				fprintf(g_ftsv, "\t%.3g", Score);
				}
			fprintf(g_ftsv, "\n");
			}
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
		if (++g_DoneCount%100 == 0)
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
	ProgressLog("%u done, %u hits, %s secs (%u threads, %.1f/ sec/ thread)\n",
	  g_DoneCount, g_HitCount, IntToStr(Secs), ThreadCount, Throughput);
	}

#include "myutils.h"
#include "pdbchain.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "cmpsearcher.h"
#include "shapesearcher.h"
#include <time.h>

static uint g_DoneCount;
static uint g_HitCount;

static void Search1(const PDBChain &Q, ShapeSearcher &SS,
  uint PosA, uint PosB, uint PosC)
	{

	SS.SetQuery(Q, PosA, PosB, PosC);

	uint SSPosA, SSPosB, SSPosC;
	SS.SearchABC(SSPosA, SSPosB, SSPosC);

	const uint AL = SS.GetShapeLength(SS.m_ShapeIndexA);
	const uint BL = SS.GetShapeLength(SS.m_ShapeIndexB);
	const uint CL = SS.GetShapeLength(SS.m_ShapeIndexC);
	string SeqA = ".";
	string SeqB = ".";
	string SeqC = ".";
	if (SSPosA != UINT_MAX)
		Q.GetSubSeq(SSPosA, AL, SeqA);
	if (SSPosB != UINT_MAX)
		Q.GetSubSeq(SSPosB, BL, SeqB);
	if (SSPosC != UINT_MAX)
		Q.GetSubSeq(SSPosC, CL, SeqC);

	Log(">%s", Q.m_Label.c_str());
	Log(" A%u/%u=%s  B %u/%u=%s  C %u/%u=%s",
	  PosA, SSPosA, SeqA.c_str(),
	  PosB, SSPosB, SeqB.c_str(),
	  PosC, SSPosC, SeqC.c_str());
	if (PosA == SSPosA && PosB == SSPosB && PosC == SSPosC)
		Log(" same\n");
	else
		Log(" **DIFF**\n");

#if 0
	const uint ShapeCount = SS.GetShapeCount();

	vector<uint> PosVec;
	PosVec.push_back(UINT_MAX); // F
	//PosVec.push_back(UINT_MAX); // H
	PosVec.push_back(PosA);
	PosVec.push_back(PosB);
	PosVec.push_back(PosC);
	PosVec.push_back(UINT_MAX); // D
	PosVec.push_back(UINT_MAX); // E
	asserta(SIZE(PosVec) == ShapeCount);

	uint MinLo, MaxHi;
	uint QL = Q.GetSeqLength();
	Log("\n");
	Log(">%s\n", Q.m_Label.c_str());
//	Q.WriteSeqWithCoords(g_fLog);
	vector<uint> HitPosVec;
	vector<double> HitScores;

	for (uint k = 0; k < 3; ++k)
		{
		uint SIX = k + 2;
		if (k == 1)
			SS.SearchShapeSelf(SIX, 0.7, 20, QL-20, 'G', 1, HitPosVec, HitScores);
		else
			SS.SearchShapeSelf(SIX, 0.7, 20, QL-20, 'D', 3, HitPosVec, HitScores);
		const uint n = SIZE(HitPosVec);
		uint SL = SS.GetShapeLength(SIX);
		Log ("SIX [%u]", n);
		for (uint i = 0; i < n; ++i)
			{
			uint Pos = HitPosVec[i];
			string ShapeSeq;
			Q.GetSubSeq(Pos, SL, ShapeSeq);
			Log(" %u(%.4f)%s", Pos, HitScores[i], ShapeSeq.c_str());
			}
		Log("\n");
		}

	Log("A=%u B=%u C=%u QL=%u\n", PosA, PosB, PosC, QL);
	vector<string> ShapeSeqs;
	for (uint ShapeIndex = 0; ShapeIndex < 7; ++ShapeIndex)
		{
		string TopShapeSeq = ".";
		double TopScore = 0;
		SS.GetMinLoMaxHi(ShapeIndex, PosVec, MinLo, MaxHi);
		//SS.SearchShapeSelf(ShapeIndex, 0.5, MinLo, MaxHi, 0, UINT_MAX,
		//	HitPosVec, HitScores);
		SS.SearchShape(ShapeIndex, PosVec, 0.5, MinLo, MaxHi, 0, UINT_MAX,
			HitPosVec, HitScores);
		Log(">> Shape %s  MinLo %5u  MaxHi  %5u: ",
			SS.GetShapeName(ShapeIndex), MinLo, MaxHi);
		const uint n = SIZE(HitPosVec);
		uint SL = SS.GetShapeLength(ShapeIndex);
		Log ("[%u]", n);
		for (uint i = 0; i < n; ++i)
			{
			uint Pos = HitPosVec[i];
			string ShapeSeq;
			Q.GetSubSeq(Pos, SL, ShapeSeq);
			double Score = HitScores[i];
			if (Score > TopScore)
				{
				TopShapeSeq = ShapeSeq;
				TopScore = Score;
				}
			Log(" %u(%.4f)%s", Pos, Score, ShapeSeq.c_str());
			}
		Log("\n");

		ShapeSeqs.push_back(TopShapeSeq);
		}
	
	if (g_ftsv == 0)
		return;

	fprintf(g_ftsv, "%s", Q.m_Label.c_str());
	fprintf(g_ftsv, "\t%u", PosA+1);
	fprintf(g_ftsv, "\t%u", PosB+1);
	fprintf(g_ftsv, "\t%u", PosC+1);
	fprintf(g_ftsv, "\t%s", (PosA < PosC ? "ABC" : "CAB"));
	for (uint i = 0; i < ShapeCount; ++i)
		fprintf(g_ftsv, "\t%s", ShapeSeqs[i].c_str());
	fprintf(g_ftsv, "\n");
#endif
	}


static void Thread(ChainReader &CR, const Shapes &S, const CMP &Prof)
	{
	CMPSearcher CS;
	CS.SetProf(Prof);

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
		uint PosA = UINT_MAX;
		uint PosB = UINT_MAX;
		uint PosC = UINT_MAX;
		CS.Search(Q);
		double PalmScore = CS.GetPosABC(PosA, PosB, PosC);
		if (PalmScore == 0)
			continue;
		Search1(Q, SS, PosA, PosB, PosC);
		++g_HitCount;
		}
	}

void cmd_shapes_search()
	{
	const string &QueryFN = opt_shapes_search;

	time_t tStart = time(0);
	if (!optset_model)
		Die("Must specify -model");
	if (!optset_shapes)
		Die("Must specify -shapes");
	const string &ModelFileName = opt_model;
	const string &ShapesFileName = opt_shapes;

	CMP Prof;
	Prof.FromFile(ModelFileName);

	Shapes S;
	S.FromFile(ShapesFileName);

	ChainReader CR;
	CR.Open(QueryFN);

	uint ThreadCount = GetRequestedThreadCount();

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, S, Prof);
	time_t tEnd = time(0);
	uint Secs = uint(tEnd - tStart);
	if (Secs <= 0)
		Secs = 1;
	double Throughput = double(g_DoneCount)/(Secs*ThreadCount);
	ProgressLog("%u done, %u hits, %s secs (%u threads, %.1f/ sec/ thread)\n",
	  g_DoneCount, g_HitCount, IntToStr(Secs), ThreadCount, Throughput);
	}

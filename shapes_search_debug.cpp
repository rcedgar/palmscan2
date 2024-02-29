#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "shapesearcher.h"
#include "rdrpsearcher.h"
#include "rdrpmodel.h"

static const double MIN_PSSM_PALM_SCORE = 15.0;

static bool Search1(const PDBChain &Q,
  RdRpSearcher &RS, ShapeSearcher &SS)
	{
	const string &Label = Q.m_Label;
	const string &Seq = Q.m_Seq;
	const uint L = SIZE(Seq);

	RS.Search(Label, Seq);
	if (RS.m_TopPalmHit.m_Score < MIN_PSSM_PALM_SCORE)
		{
		Log(">%s No PSSM hit score %.3g\n",
		  Label.c_str(), RS.m_TopPalmHit.m_Score);
		return false;
		}

	string PSSM_A, PSSM_B, PSSM_C;
	RS.GetShapesTrainABC(PSSM_A, PSSM_B, PSSM_C);

	size_t stPosA = Seq.find(PSSM_A);
	size_t stPosB = Seq.find(PSSM_B);
	size_t stPosC = Seq.find(PSSM_C);
	asserta(stPosA < L);
	asserta(stPosB < L);
	asserta(stPosC < L);
	if (stPosC < stPosA)
		{
		Log("   >%s permuted\n", Label.c_str());
		return false;
		}

	uint PSSM_PosA = uint(stPosA);
	uint PSSM_PosB = uint(stPosB);
	uint PSSM_PosC = uint(stPosC);
	Log(">%s PSSM A:%u:%s B:%u:%s C:%u:%s\n",
	  Label.c_str(),
	  PSSM_PosA,
	  PSSM_A.c_str(),
	  PSSM_PosB,
	  PSSM_B.c_str(),
	  PSSM_PosC,
	  PSSM_C.c_str());

	vector<uint> ABCIndexes;
	ABCIndexes.push_back(SS.m_ShapeIndexA);
	ABCIndexes.push_back(SS.m_ShapeIndexB);
	ABCIndexes.push_back(SS.m_ShapeIndexC);

	vector<uint> PosVec(SS.m_ShapeCount, UINT_MAX);
	PosVec[SS.m_ShapeIndexA] = PSSM_PosA;
	PosVec[SS.m_ShapeIndexB] = PSSM_PosB;
	PosVec[SS.m_ShapeIndexC] = PSSM_PosC;

	SS.SetQuery(Q);
	double ScoreA = SS.GetSelfScore(SS.m_ShapeIndexA, PSSM_PosA);
	double ScoreB = SS.GetSelfScore(SS.m_ShapeIndexB, PSSM_PosB);
	double ScoreC = SS.GetSelfScore(SS.m_ShapeIndexC, PSSM_PosC);
	double ScoreABC = SS.GetScoreShapes(PosVec);
	Log("  ScoreABC PSSM pos A=%.4g B=%.4g C=%.4g ABC=%.4g\n",
	  ScoreA, ScoreB, ScoreC, ScoreABC);

	SS.SearchABC(true);
	const uint ShapeCount = SS.GetShapeCount();
	asserta(SIZE(SS.m_ShapePosVec) == ShapeCount);
	asserta(SIZE(SS.m_ShapeScores) == ShapeCount);

	uint IXA = SS.m_ShapeIndexA;
	uint IXB = SS.m_ShapeIndexB;
	uint IXC = SS.m_ShapeIndexC;
	bool IsHit = (SS.m_ABCScore >= SS.m_MinABCScore);
	if (opt_misses || IsHit)
		SS.ToTsv(g_ftsv);
	return IsHit;
	}

void cmd_shapes_search_debug()
	{
	const string &QueryFN = opt_shapes_search_debug;

	Shapes S;
	S.InitFromCmdLine();
	ShapeSearcher SS;
	SS.Init(S);

	ChainReader CR;
	CR.Open(QueryFN);

	RdRpModel Model;
	GetRdrpModel(Model);

	RdRpSearcher RS;
	RS.Init(Model);

	PDBChain Q;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

		const char *Seq = Q.m_Seq.c_str();
		const string &Label = Q.m_Label;

		vector<uint> PosVec;
		vector<string> ShapeSeqs;
		bool IsHit = Search1(Q, RS, SS);
		ProgressLog("%s hit %c\n", Q.m_Label.c_str(), yon(IsHit));
		}
	}

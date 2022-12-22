#include "myutils.h"
#include "pdbchain.h"
#include "rdrpmodel.h"
#include "rdrpsearcher.h"
#include "ppspsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "ppsp.h"

static uint g_DoneCount;
static uint g_HitCount;

static double g_MinScoreA = DBL_MAX;
static double g_MinScoreB = DBL_MAX;
static double g_MinScoreC = DBL_MAX;

static double g_MinScoreAB = DBL_MAX;
static double g_MinScoreBC = DBL_MAX;
static double g_MinScoreAC = DBL_MAX;

static double g_MinScore3 = DBL_MAX;

static void Train(ChainReader &CR, RdRpSearcher &RS, PPSP &Prof)
	{
	g_DoneCount = 0;
	g_HitCount = 0;

	g_MinScoreA = DBL_MAX;
	g_MinScoreB = DBL_MAX;
	g_MinScoreC = DBL_MAX;

	g_MinScoreAB = DBL_MAX;
	g_MinScoreBC = DBL_MAX;
	g_MinScoreAC = DBL_MAX;

	g_MinScore3 = DBL_MAX;

	PDBChain Q;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

		if (++g_DoneCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("Train %s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}

		const string &QSeq = Q.m_Seq;
		const string &QLabel = Q.m_Label;

		RS.Search(QLabel, QSeq);
		RS.WriteOutput();

		uint APos = RS.GetMotifPos(0);
		uint BPos = RS.GetMotifPos(1);
		uint CPos = RS.GetMotifPos(2);
		if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
			continue;

		++g_HitCount;
		Prof.TrainChain(Q, APos, BPos, CPos);
		}
	Progress("Train 100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);
	Prof.FinalizeTrain();
	}

static void Write1(const char *Method, const PDBChain &Q, const PPSP &Prof,
  uint APos, uint BPos, uint CPos)
	{
	if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
		{
		Log("%6.6s", ".");
		Log("  %4.4s", ".");
		Log("  %4.4s", ".");
		Log("  %4.4s", ".");
		Log("  %6.6s", ".");
		Log("  %6.6s", ".");
		Log("  %6.6s", ".");
		Log("  %6.6s", ".");
		Log("  %6.6s", ".");
		Log("  %6.6s", ".");
		Log("  %s", Method);
		Log("  >%s\n", Q.m_Label.c_str());
		return;
		}

	double ScoreA = Prof.GetScoreA(Q, APos);
	double ScoreB = Prof.GetScoreB(Q, BPos);
	double ScoreC = Prof.GetScoreC(Q, CPos);

	double ScoreAB = Prof.GetScoreAB(Q, APos, BPos);
	double ScoreBC = Prof.GetScoreBC(Q, BPos, CPos);
	double ScoreAC = Prof.GetScoreAC(Q, APos, CPos);

	g_MinScoreA = min(ScoreA, g_MinScoreA);
	g_MinScoreB = min(ScoreB, g_MinScoreB);
	g_MinScoreC = min(ScoreC, g_MinScoreC);

	g_MinScoreAB = min(ScoreAB, g_MinScoreAB);
	g_MinScoreBC = min(ScoreBC, g_MinScoreBC);
	g_MinScoreAC = min(ScoreAC, g_MinScoreAC);

	double Score3 = 0;

	Score3 += ScoreA;
	Score3 += ScoreB;
	Score3 += ScoreC;

	Score3 += ScoreAB;
	Score3 += ScoreBC;
	Score3 += ScoreAC;

	Score3 /= 6;

	g_MinScore3 = min(Score3, g_MinScore3);

	Log("%6.4f", Score3);
	Log("  %4u", APos);
	Log("  %4u", BPos);
	Log("  %4u", CPos);
	Log("  %6.4f", ScoreA);
	Log("  %6.4f", ScoreB);
	Log("  %6.4f", ScoreC);
	Log("  %6.4f", ScoreAB);
	Log("  %6.4f", ScoreBC);
	Log("  %6.4f", ScoreAC);
	Log("  %s", Method);
	Log("  >%s\n", Q.m_Label.c_str());
	}

static void Test(ChainReader &CR, RdRpSearcher *RS, 
  PPSPSearcher *PS, PPSP &Prof)
	{
	g_DoneCount = 0;
	g_HitCount = 0;

	PDBChain Q;

	Log("%6.6s", "\nScore");
	Log("  %4.4s", "A");
	Log("  %4.4s", "B");
	Log("  %4.4s", "C");
	Log("  %6.6s", "A");
	Log("  %6.6s", "B");
	Log("  %6.6s", "C");
	Log("  %6.6s", "AB");
	Log("  %6.6s", "BC");
	Log("  %6.6s", "AC");
	Log("  %s\n", "Label");
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

		if (++g_DoneCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("Test %s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}

		const string &QSeq = Q.m_Seq;
		const string &QLabel = Q.m_Label;

		uint APos_RS = UINT_MAX;
		uint BPos_RS = UINT_MAX;
		uint CPos_RS = UINT_MAX;
		if (RS != 0)
			{
			RS->Search(QLabel, QSeq);
			RS->WriteOutput();

			APos_RS = RS->GetMotifPos(0);
			BPos_RS = RS->GetMotifPos(1);
			CPos_RS = RS->GetMotifPos(2);
			Write1("PSSM", Q, Prof, APos_RS, BPos_RS, CPos_RS);
			}

		uint APos_PS = UINT_MAX;
		uint BPos_PS = UINT_MAX;
		uint CPos_PS = UINT_MAX;
		if (PS != 0)
			{
			PS->Search(Q);
			PS->GetPSSMStarts(APos_PS, BPos_PS, CPos_PS);
			const char *Method = "PPSP";
			if (RS != 0)
				{
				if (APos_PS == APos_RS
				 && BPos_PS == BPos_RS
				 && CPos_PS == CPos_RS)
					Method = "SAME";
				else
					Method = "*DIF";
				}
			Write1(Method, Q, Prof, APos_PS, BPos_PS, CPos_PS);
			}

		if (APos_PS != UINT_MAX || APos_RS != UINT_MAX)
			++g_HitCount;
		if (RS != 0 && PS != 0)
			Log("\n");
		}

	Log("\n");
	if (g_HitCount > 0)
		{
		Log("%6.4f", g_MinScore3);
		Log("  %4.4s", ".");
		Log("  %4.4s", ".");
		Log("  %4.4s", ".");
		Log("  %6.4f", g_MinScoreA);
		Log("  %6.4f", g_MinScoreB);
		Log("  %6.4f", g_MinScoreC);
		Log("  %6.4f", g_MinScoreAB);
		Log("  %6.4f", g_MinScoreBC);
		Log("  %6.4f", g_MinScoreAC);
		Log("  Minimum\n");
		}

	Progress("Test 100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);
	Prof.FinalizeTrain();

	Log("%u / %u %s hits\n",
	  g_HitCount, g_DoneCount, (RS == 0 ? "PPSP" : "PSSM"));
	}

void cmd_ppsp_train_pssm()
	{
	const string &QueryFN = opt_ppsp_train_pssm;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;
	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

	PDBChain Q;
	RdRpSearcher RS;
	RS.Init(Model);

	ChainReader CR;
	CR.Open(QueryFN, false);

	PPSP Prof;
	Train(CR, RS, Prof);
	Prof.ToFile(opt_output);

	CR.Clear();
	CR.Open(QueryFN, false);

	PPSPSearcher PS;
	PS.m_Prof = &Prof;

	CR.Clear();
	CR.Open(QueryFN, false);
	Test(CR, &RS, &PS, Prof);
	}

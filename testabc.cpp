#include "myutils.h"
#include "shapesearcher.h"
#include "motifsettings.h"

static int MAX_SHIFT = 4;

static uint g_NoRefCount;
static uint g_NoPredCount;
static uint g_AgreeCount;
static uint g_DisagreeCount;
static uint g_RefPermutedCount;
static uint g_PredScoreHigherCount;
static uint g_PredScoreLowerCount;
static uint g_PredScoreTooLowCount;
static uint g_RevHitsCount;

void ShapeSearcher::LogShapes(const vector<uint> &ShapeIndexes,
  const vector<string> &ShapeSeqs) const
	{
	const uint N = SIZE(ShapeIndexes);
	asserta(SIZE(ShapeSeqs) == N);

	Log(">%s", m_Query->m_Label.c_str());
	const string &Seq = m_Query->m_Seq;
	vector<uint> PosVec;
	for (uint i = 0; i < N; ++i)
		{
		uint ShapeIndex = ShapeIndexes[i];
		const string &ShapeSeq = ShapeSeqs[i];
		size_t stPos = Seq.find(ShapeSeq);
		uint Pos = (stPos == string::npos ? UINT_MAX : uint(stPos));
		const char *ShapeName = GetShapeName(ShapeIndex);
		Log(" %s:%s", ShapeName, ShapeSeq.c_str());
		if (Pos != UINT_MAX)
			{
			double Score = GetSelfScore(ShapeIndex, Pos);
			Log("(%.3f)", Score);
			}
		}
	Log("\n");
	}

void ShapeSearcher::LogShape(const string &Msg, uint ShapeIndex, uint Pos) const
	{
	Log("%s Shape_", Msg.c_str());
	if (ShapeIndex == UINT_MAX)
		{
		Log("*  ");
		return;
		}
	Log("%s", GetShapeName(ShapeIndex));
	if (Pos == UINT_MAX)
		{
		Log(" pos=.  ", Pos);
		return;
		}
	Log(" pos=%u", Pos);
	string ShapeSeq;
	uint ShapeLength = GetShapeLength(ShapeIndex);
	GetSubSeq(Pos, ShapeLength, ShapeSeq);
	Log(" %s", ShapeSeq.c_str());
	double Score = GetSelfScore(ShapeIndex, Pos);
	Log(" score %.3g", Score);
	Log("  ");
	}

void ShapeSearcher::TestABC(const Shapes &S,
  const vector<PDBChain *> &Chains,
  vector<vector<string> > &MotifSeqsVec,
  double MinPredScore)
	{
	const uint N = SIZE(Chains);
	asserta(SIZE(MotifSeqsVec) == N);

	g_AgreeCount = 0;
	ShapeSearcher SS;
	SS.Init(S);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Test ABC");
		SS.TestABC1(*Chains[i], MotifSeqsVec[i], MinPredScore);
		}

	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Test ABC reverse");
		const PDBChain &Chain = *Chains[i];
		PDBChain RevChain;
		Chain.GetReverse(RevChain);
		SS.SetQuery(RevChain);
		SS.SearchABC();
		if (SS.m_ABCScore >= MinPredScore)
			++g_RevHitsCount;
		}

	ProgressLog("%10u  Sequences (min score %.3g)\n", N, MinPredScore);
#define X(x)	ProgressLog("%10u  (%5.1f%%)  %s\n", g_##x##Count, GetPct(g_##x##Count, N), #x)
	X(Agree);
	X(Disagree);
	X(PredScoreHigher);
	X(PredScoreLower);
	X(RefPermuted);
	X(NoRef);
	X(NoPred);
	X(PredScoreTooLow);
	X(RevHits);
#undef X
	}

static bool ShiftOk(uint RefPos, uint PredPos)
	{
	int Shift = abs(int(RefPos) - int(PredPos));
	return Shift <= MAX_SHIFT;
	}

void ShapeSearcher::TestABC1(const PDBChain &Chain,
  const vector<string> &RefMotifSeqs, double MinPredScore)
	{
	const char *Label = Chain.m_Label.c_str();
	const string &Seq = Chain.m_Seq;
	size_t stPosA = Seq.find(RefMotifSeqs[m_ShapeIndexA]);
	size_t stPosB = Seq.find(RefMotifSeqs[m_ShapeIndexB]);
	size_t stPosC = Seq.find(RefMotifSeqs[m_ShapeIndexC]);

	uint RefPosA = (stPosA == string::npos ? UINT_MAX : uint(stPosA));
	uint RefPosB = (stPosB == string::npos ? UINT_MAX : uint(stPosB));
	uint RefPosC = (stPosC == string::npos ? UINT_MAX : uint(stPosC));

	if (RefPosC != UINT_MAX && RefPosA != UINT_MAX && RefPosC < RefPosA)
		{
		++g_RefPermutedCount;
		return;
		}

	SetQuery(Chain);

	if (RefPosA == UINT_MAX || RefPosB == UINT_MAX || RefPosC == UINT_MAX)
		{
		Log("No ref>%s\n", Chain.m_Label.c_str());
		++g_NoRefCount;
		return;
		}

	SearchABC();
	if (m_ABCScore == 0)
		{
		Log("No pred (zero score) >%s\n", Label);
		++g_NoPredCount;
		return;
		}

	if (m_ABCScore < MinPredScore)
		{
		Log("Pred score too low %6.4f >%s\n", m_ABCScore, Label);
		++g_PredScoreTooLowCount;
		return;
		}

	uint PredPosA = m_ShapePosVec[m_ShapeIndexA];
	uint PredPosB = m_ShapePosVec[m_ShapeIndexB];
	uint PredPosC = m_ShapePosVec[m_ShapeIndexC];

	{
	string A, B, C;
	GetShapeSeq(m_ShapeIndexA, A);
	GetShapeSeq(m_ShapeIndexB, B);
	GetShapeSeq(m_ShapeIndexC, C);
	CheckA(A);
	CheckB(B);
	CheckC(C);
	}

	if (ShiftOk(PredPosA, RefPosA) &&
	  ShiftOk(PredPosB, RefPosB) &&
	  ShiftOk(PredPosC,RefPosC))
		{
		++g_AgreeCount;
		return;
		}

	if (PredPosA == UINT_MAX || PredPosB == UINT_MAX || PredPosC == UINT_MAX)
		{
		++g_NoPredCount;
		return;
		}

	vector<uint> ShapeIndexes;
	ShapeIndexes.push_back(m_ShapeIndexA);
	ShapeIndexes.push_back(m_ShapeIndexB);
	ShapeIndexes.push_back(m_ShapeIndexC);

	vector<uint> PredPosVec(m_ShapeCount, UINT_MAX);
	PredPosVec[m_ShapeIndexA] = PredPosA;
	PredPosVec[m_ShapeIndexB] = PredPosB;
	PredPosVec[m_ShapeIndexC] = PredPosC;
	double PredScore = GetScoreShapes(PredPosVec);

	vector<uint> RefPosVec(m_ShapeCount, UINT_MAX);
	RefPosVec[m_ShapeIndexA] = RefPosA;
	RefPosVec[m_ShapeIndexB] = RefPosB;
	RefPosVec[m_ShapeIndexC] = RefPosC;
	double RefScore = GetScoreShapes(RefPosVec);

	if (m_ABCScore > RefScore)
		{
		++g_PredScoreHigherCount;
		return;
		}
	if (m_ABCScore < RefScore)
		++g_PredScoreLowerCount;

	++g_DisagreeCount;

	Log(">%s disagree ref %.3g, pred %.3g\n",
	  m_Query->m_Label.c_str(), RefScore, m_ABCScore);
#if 0
	{
	LogShape("RefA", m_ShapeIndexA, RefPosA);
	if (PredPosA != RefPosA)
		LogShape("PrdA", m_ShapeIndexA, PredPosA);
	Log("\n");

	LogShape("RefB", m_ShapeIndexB, RefPosB);
	if (PredPosB != RefPosB)
		LogShape("PrdB", m_ShapeIndexB, PredPosB);
	Log("\n");

	LogShape("RefC", m_ShapeIndexC, RefPosC);
	if (PredPosC != RefPosC)
		LogShape("PrdC", m_ShapeIndexC, PredPosC);
	Log("\n");

	Log("reinitialize\n");
	Log("load d:/int/ictv_rdrp_structures/pdb/%s.pdb\n",
	  m_Query->m_Label.c_str());
	Log("color gray70\n");

	Log("select pepseq %s\n", RefMotifSeqs[0].c_str());
	Log("color tv_blue, sele\n");
	if (PredPosA != RefPosA)
		{
		string PredSeqA;
		uint LA = GetShapeLength(m_ShapeIndexA);
		GetSubSeq(PredPosA, LA, PredSeqA);
		Log("select pepseq %s\n", PredSeqA.c_str());
		Log("color slate, sele\n");
		}

	Log("select pepseq %s\n", RefMotifSeqs[1].c_str());
	Log("color tv_green, sele\n");
	if (PredPosB != RefPosB)
		{
		string PredSeqB;
		uint LB = GetShapeLength(m_ShapeIndexB);
		GetSubSeq(PredPosB, LB, PredSeqB);
		Log("select pepseq %s\n", PredSeqB.c_str());
		Log("color splitpea, sele\n");
		}

	Log("select pepseq %s\n", RefMotifSeqs[2].c_str());
	Log("color tv_red, sele\n");
	if (PredPosC != RefPosC)
		{
		string PredSeqC;
		uint LC = GetShapeLength(m_ShapeIndexC);
		GetSubSeq(PredPosC, LC, PredSeqC);
		Log("select pepseq %s\n", PredSeqC.c_str());
		Log("color salmon, sele\n");
		}
	Log("deselect\n");
	}
#endif
	}

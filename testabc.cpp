#include "myutils.h"
#include "shapesearcher.h"

static uint g_AgreeCount;
static uint g_DisagreeCount;
static uint g_RefPermutedCount;

void ShapeSearcher::TestABC(const Shapes &S,
  const vector<PDBChain *> &Chains,
  vector<vector<string> > &MotifSeqsVec)
	{
	const uint N = SIZE(Chains);
	asserta(SIZE(MotifSeqsVec) == N);

	g_AgreeCount = 0;

	ShapeSearcher SS;
	SS.Init(S);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Test ABC");
		SS.TestABC1(*Chains[i], MotifSeqsVec[i]);
		}

	ProgressLog("Agree %u, disagree %u, permuted %u\n",
	  g_AgreeCount, g_DisagreeCount, g_RefPermutedCount);
	}

void ShapeSearcher::TestABC1(const PDBChain &Chain,
  const vector<string> &MotifSeqs)
	{
	const string &Seq = Chain.m_Seq;
	size_t stPosA = Seq.find(MotifSeqs[m_ShapeIndexA]);
	size_t stPosB = Seq.find(MotifSeqs[m_ShapeIndexB]);
	size_t stPosC = Seq.find(MotifSeqs[m_ShapeIndexC]);
	if (stPosA == string::npos)
		return;
	if (stPosB == string::npos)
		return;
	if (stPosC == string::npos)
		return;

	if (stPosC < stPosA)
		{
		++g_RefPermutedCount;
		return;
		}

	uint RefPosA = uint(stPosA);
	uint RefPosB = uint(stPosB);
	uint RefPosC = uint(stPosC);

	SetQuery(Chain);

	uint PosA, PosB, PosC;
	SearchABC(PosA, PosB, PosC);

	if (PosA == RefPosA && PosB == RefPosB && PosC == RefPosC)
		{
		++g_AgreeCount;
		return;
		}
	++g_DisagreeCount;

	string RefSeqA;
	string RefSeqB;
	string RefSeqC;

	string SeqA;
	string SeqB;
	string SeqC;

	uint AL = GetShapeLength(m_ShapeIndexA);
	uint BL = GetShapeLength(m_ShapeIndexB);
	uint CL = GetShapeLength(m_ShapeIndexC);

	Chain.GetSubSeq(RefPosA, AL, RefSeqA);
	Chain.GetSubSeq(RefPosB, BL, RefSeqB);
	Chain.GetSubSeq(RefPosC, CL, RefSeqC);

	Chain.GetSubSeq(PosA, AL, SeqA);
	Chain.GetSubSeq(PosB, BL, SeqB);
	Chain.GetSubSeq(PosC, CL, SeqC);

	Log("\n");
	Log(">%s\n", Chain.m_Label.c_str());
	Log("A: %4u %4u  %s  %s", RefPosA, PosA, RefSeqA.c_str(), SeqA.c_str());
	if (RefPosA != PosA)
		Log(" << DIFF");
	Log("\n");

	Log("B: %4u %4u  %s  %s", RefPosB, PosB, RefSeqB.c_str(), SeqB.c_str());
	if (RefPosB != PosB)
		Log(" << DIFF");
	Log("\n");

	Log("C: %4u %4u  %s  %s", RefPosC, PosC, RefSeqC.c_str(), SeqC.c_str());
	if (RefPosC != PosC)
		Log(" << DIFF");
	Log("\n");
	}

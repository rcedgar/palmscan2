#include "myutils.h"
#include "shapesearcher.h"

static uint g_AgreeCount;
static uint g_DisagreeCount;
static uint g_RefPermutedCount;

static bool Cmp1(ShapeSearcher &SS, const PDBChain &Chain,
  uint ShapeIndex, uint RefPosX, uint PredPosX)
	{
	if (RefPosX == UINT_MAX && PredPosX == UINT_MAX)
		return true;

	string RefSeqX = ".";
	string PredSeqX = ".";

	uint XL = SS.GetShapeLength(ShapeIndex);

	const char *Label = Chain.m_Label.c_str();
	if (RefPosX != UINT_MAX)
		Chain.GetSubSeq(RefPosX, XL, RefSeqX);

	if (PredPosX != UINT_MAX)
		Chain.GetSubSeq(PredPosX, XL, PredSeqX);

	if (RefSeqX == PredSeqX)
		return true;

	Log(">%s  Motif-%u", Label, ShapeIndex);
	if (RefPosX == UINT_MAX)
		Log(" ref=.");
	else
		Log("  ref=%s(%u)", RefSeqX.c_str(), RefPosX);

	if (PredPosX == UINT_MAX)
		Log(" pred=.");
	else
		Log("  pred=%s(%u)", PredSeqX.c_str(), PredPosX);

	if (RefPosX == UINT_MAX && PredPosX != UINT_MAX)
		Log(" << found");

	int Shift = 999;
	if (RefPosX != UINT_MAX && PredPosX != UINT_MAX)
		{
		Shift = int(PredPosX) - int(RefPosX);
		Log(" << shift %+d", Shift);
		}
	if (RefPosX != UINT_MAX && PredPosX == UINT_MAX)
		Log(" << NOT FOUND");


	Log("\n");

	if (abs(Shift) <= 3)
		return true;
	if (RefPosX == UINT_MAX)
		return true;
	return false;
	}

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

	uint RefPosA = (stPosA == string::npos ? UINT_MAX : uint(stPosA));
	uint RefPosB = (stPosB == string::npos ? UINT_MAX : uint(stPosB));
	uint RefPosC = (stPosC == string::npos ? UINT_MAX : uint(stPosC));

	if (RefPosC != UINT_MAX && RefPosA != UINT_MAX && RefPosC < RefPosA)
		{
		++g_RefPermutedCount;
		return;
		}

	SetQuery(Chain);

	uint PredPosA, PredPosB, PredPosC;
	SearchABC(PredPosA, PredPosB, PredPosC);

	if (PredPosA == RefPosA && PredPosB == RefPosB && PredPosC == RefPosC)
		{
		Log(">%s agree ABC\n", Chain.m_Label.c_str());
		++g_AgreeCount;
		return;
		}

	bool OkA = Cmp1(*this, Chain, m_ShapeIndexA, RefPosA, PredPosA);
	bool OkB = Cmp1(*this, Chain, m_ShapeIndexB, RefPosB, PredPosB);
	bool OkC = Cmp1(*this, Chain, m_ShapeIndexC, RefPosC, PredPosC);
	if (!OkA || !OkB || !OkC)
		++g_DisagreeCount;
	}

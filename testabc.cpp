#include "myutils.h"
#include "shapesearcher.h"

static uint g_AgreeCount;
static uint g_DisagreeCount;
static uint g_RefPermutedCount;

static void LogDiff(ShapeSearcher &SS, const PDBChain &Chain,
  uint ShapeIndex, uint RefPosX, uint PredPosX)
	{
	if (RefPosX == UINT_MAX && PredPosX == UINT_MAX)
		return;

	string RefSeqX = ".";
	string PredSeqX = ".";

	uint XL = SS.GetShapeLength(ShapeIndex);

	const char *Label = Chain.m_Label.c_str();
	if (RefPosX != UINT_MAX)
		Chain.GetSubSeq(RefPosX, XL, RefSeqX);

	if (PredPosX != UINT_MAX)
		Chain.GetSubSeq(PredPosX, XL, PredSeqX);

	if (RefSeqX == PredSeqX)
		return;

	Log("%s:", SS.GetShapeName(ShapeIndex));
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
		if (abs(Shift) <= 3)
			Log(" ok");
		else
			Log(" *ERROR*");
		}
	if (RefPosX != UINT_MAX && PredPosX == UINT_MAX)
		Log(" << NOT FOUND");

	Log("\n");
	}

static bool CmpOk(ShapeSearcher &SS, const PDBChain &Chain,
  uint ShapeIndex, uint RefPosX, uint PredPosX)
	{
	if (RefPosX == PredPosX)
		return true;

	if (RefPosX == UINT_MAX && PredPosX != UINT_MAX)
		return true;

	if (RefPosX != UINT_MAX && PredPosX == UINT_MAX)
		return false;

	asserta(RefPosX != UINT_MAX && PredPosX != UINT_MAX);
	int Shift = int(PredPosX) - int(RefPosX);
	if (abs(Shift) <= 3)
		return true;
	else
		return false;
	}

static bool Cmp1(ShapeSearcher &SS, const PDBChain &Chain,
  uint ShapeIndex, uint RefPosX, uint PredPosX)
	{
	bool Ok = CmpOk(SS, Chain, ShapeIndex, RefPosX, PredPosX);
	if (Ok)
		return true;
	LogDiff(SS, Chain, ShapeIndex, RefPosX, PredPosX);
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
#if 0
	{
	uint MinDistAB, MaxDistAB;
	uint MinDistBC, MaxDistBC;
	SS.GetDistRange(SS.m_ShapeIndexA, SS.m_ShapeIndexB, MinDistAB, MaxDistAB);
	SS.GetDistRange(SS.m_ShapeIndexB, SS.m_ShapeIndexC, MinDistBC, MaxDistBC);

	MinDistAB *= 0.8;
	MaxDistAB *= 1.2;

	MinDistBC *= 0.8;
	MaxDistBC *= 1.2;
	Log("DistAB = %u .. %u, BC = %u .. %u\n",
	  MinDistAB, MaxDistAB, MinDistBC, MaxDistBC);
	}
#endif // 0
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
		Log("Permuted >%s\n", Chain.m_Label.c_str());
		++g_RefPermutedCount;
		return;
		}

	SetQuery(Chain);

	SearchABC();
	uint PredPosA = m_ShapePosVec[m_ShapeIndexA];
	uint PredPosB = m_ShapePosVec[m_ShapeIndexB];
	uint PredPosC = m_ShapePosVec[m_ShapeIndexC];

	if (PredPosA == RefPosA && PredPosB == RefPosB && PredPosC == RefPosC)
		{
		Log(">%s agree ABC\n", Chain.m_Label.c_str());
		++g_AgreeCount;
		return;
		}

	Log("\n");
	Log(">%s\n", Chain.m_Label);

	if (
	  PredPosA != UINT_MAX &&
	  PredPosB != UINT_MAX &&
	  PredPosC != UINT_MAX &&
	  RefPosA != UINT_MAX &&
	  RefPosB != UINT_MAX &&
	  RefPosC != UINT_MAX)
		{
		vector<uint> ShapeIndexes;
		ShapeIndexes.push_back(m_ShapeIndexA);
		ShapeIndexes.push_back(m_ShapeIndexB);
		ShapeIndexes.push_back(m_ShapeIndexC);

		vector<uint> RefPosVec;
		RefPosVec.push_back(RefPosA);
		RefPosVec.push_back(RefPosB);
		RefPosVec.push_back(RefPosC);
		double RefScore = GetScoreShapes(ShapeIndexes, RefPosVec);

		const string &Seq = m_Query->m_Seq;
		double RefScoreSelfA = GetSelfScore(m_ShapeIndexA, RefPosA);
		double RefScoreSelfB = GetSelfScore(m_ShapeIndexB, RefPosB);
		double RefScoreSelfC = GetSelfScore(m_ShapeIndexC, RefPosC);

		char cA = Seq[RefPosA + m_OffsetABCs[0]];
		char cB = Seq[RefPosB + m_OffsetABCs[1]];
		char cC = Seq[RefPosC + m_OffsetABCs[2]];

		Log("RefA pos %u char %c self %.3g\n", RefPosA, cA, RefScoreSelfA);
		Log("RefB pos %u char %c self %.3g\n", RefPosB, cB, RefScoreSelfB);
		Log("RefC pos %u char %c self %.3g\n", RefPosC, cC, RefScoreSelfC);

		double PredScoreSelfA = GetSelfScore(m_ShapeIndexA, PredPosA);
		double PredScoreSelfB = GetSelfScore(m_ShapeIndexB, PredPosB);
		double PredScoreSelfC = GetSelfScore(m_ShapeIndexC, PredPosC);
		Log("PredA pos %u self %.3g\n", PredPosA, PredScoreSelfA);
		Log("PredB pos %u self %.3g\n", PredPosB, PredScoreSelfB);
		Log("PredC pos %u self %.3g\n", PredPosC, PredScoreSelfC);

		vector<uint> PredPosVec;
		PredPosVec.push_back(PredPosA);
		PredPosVec.push_back(PredPosB);
		PredPosVec.push_back(PredPosC);
		double PredScore = GetScoreShapes(ShapeIndexes, PredPosVec);
		Log("RefScore %.3g AB=%u BC=%u, PredScore %.3g\n",
		  RefScore, RefPosB - RefPosA, RefPosC - RefPosB, PredScore);
		}

	bool OkA = Cmp1(*this, Chain, m_ShapeIndexA, RefPosA, PredPosA);
	bool OkB = Cmp1(*this, Chain, m_ShapeIndexB, RefPosB, PredPosB);
	bool OkC = Cmp1(*this, Chain, m_ShapeIndexC, RefPosC, PredPosC);
	if (!OkA || !OkB || !OkC)
		++g_DisagreeCount;
	}

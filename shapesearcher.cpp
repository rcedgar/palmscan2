#include "myutils.h"
#include "shapesearcher.h"
#include "motifsettings.h"
#include "quarts.h"
#include "sort.h"

#define TRACE	0

static uint g_Queries;
static uint g_Hits;
static uint g_Misses;
static uint g_Miss_LowABCScore;
static uint g_Miss_LowPalmScore;
//static uint g_Miss_HighLEFPPM;
static uint g_Miss_MissingMotif;
static uint g_PermutedHits;
static uint g_RdRpPlus;
static uint g_RdRpMinus;
static uint g_RdRpLike;
static uint g_TPCount;
static uint g_FPCount;
static uint g_TNCount;
static uint g_FNCount;

static vector<uint> g_CalibratePosScoreToCount(101);
static vector<uint> g_CalibrateNegScoreToCount(101);
static vector<double> g_CalibrateNegScores;

static void IncStat(uint &i)
	{
#pragma omp critical
		{
		++i;
		}
	}

//void ShapeSearcher::SetLEFPPM()
//	{
//	double GetExpValue(double Value, double Mean, double StdDev,
//	  double LogDBSize);
//
//	m_LEFPPM = GetExpValue(m_FinalScore, m_MeanFinalScore,
//	  m_StdDevFinalScore, m_Log10DBSize);
//	}

void ShapeSearcher::SetFinalScore()
	{
	m_FinalScore = m_DomScore;
	}

void ShapeSearcher::CalibrateAdd(bool Hit) const
	{
	const string &Label = m_Query->m_Label;
	double Score = m_FinalScore;
	asserta(Score >= 0 && Score <= 1);
	uint iScore = uint(Score*100);
	asserta(iScore >= 0 && iScore <= 100);
	bool IsPos = StartsWith(Label, "pos.");
	bool IsNeg = StartsWith(Label, "neg.");
	if (IsPos)
		{
		if (Hit)
			IncStat(g_TPCount);
		else
			IncStat(g_FNCount);
		}
	else
		{
		if (Hit)
			IncStat(g_FPCount);
		else
			IncStat(g_TNCount);
		}

#pragma omp critical
		{
		if (IsPos)
			g_CalibratePosScoreToCount[iScore] += 1;
		else
			{
			g_CalibrateNegScoreToCount[iScore] += 1;
			if (Score > 0)
				g_CalibrateNegScores.push_back(Score);
			}
		}
	}

void ShapeSearcher::CalibrateWrite()
	{
	{
	FILE *f = CreateStdioFile("calibrate_neg_scores.txt");
	for (uint i = 0; i < SIZE(g_CalibrateNegScores); ++i)
		fprintf(f, "%.6g\n", g_CalibrateNegScores[i]);
	CloseStdioFile(f);
	}

	QuartsDouble Q;
	GetQuartsDouble(g_CalibrateNegScores, Q);
	Log("Neg scores: mean %.6g, stddev %.6g\n", Q.Avg, Q.StdDev);
	Q.LogMe();

	FILE *f = CreateStdioFile("calibrate.tsv");
	asserta(SIZE(g_CalibratePosScoreToCount) == 101);
	asserta(SIZE(g_CalibrateNegScoreToCount) == 101);
	uint NT = 0;
	uint NF = 0;
	for (int i = 0; i <= 100; ++i)
		{
		NT += g_CalibratePosScoreToCount[i];
		NF += g_CalibrateNegScoreToCount[i];
		}

	uint TotT = 0;
	uint TotF = 0;
	vector<double> Ts;
	vector<bool> Dones;
#define T(x)	Ts.push_back(x); Dones.push_back(false);
	T(1e-3);
	T(1e-4);
	T(1e-5);
	T(1e-6);
	T(1e-7);
	T(1e-8);
#undef T
	const uint TCount = SIZE(Ts);
	asserta(SIZE(Dones) == TCount);

	for (int i = 100; i > 0; --i)
		{
		uint nt = g_CalibratePosScoreToCount[i];
		uint nf = g_CalibrateNegScoreToCount[i];
		TotT += nt;
		TotF += nf;
		double FractF = NF == 0 ? 0 : double(nf)/NF;
		double FractT = NT == 0 ? 0 : double(nt)/NT;
		for (uint k = 0; k < TCount; ++k)
			{
			double T = Ts[k];
			if (!Dones[k] && FractF >= T)
				{
				Log("Err %8.3e  Score %.4f  FractT %.4g\n",
				  T, i/100.0, FractT);
				Dones[k] = true;
				}
			}
		fprintf(f, "%u\t%u\t%u\t%.5g\t%.5g\n", i, nt, nf, FractF, FractT);
		}
	CloseStdioFile(f);
	}

void ShapeSearcher::StatsToFev(FILE *f)
	{
	if (f == 0)
		return;
	string CmdLine;
	GetCmdLine(CmdLine);

	g_Misses = g_Queries - g_Hits;
	double PctHits = GetPct(g_Hits, g_Queries);
	double PctMisses = GetPct(g_Misses, g_Queries);

#define X(x)	fprintf(f, "\t%s=%u", #x, g_##x);
	X(Queries);
	X(Hits);
	X(Misses);
	X(RdRpPlus);
	X(RdRpMinus);
	X(RdRpLike);
	X(PermutedHits);
	X(Miss_LowABCScore);
	X(Miss_LowPalmScore);
	//X(Miss_HighLEFPPM);
	X(Miss_MissingMotif);
#undef X

#define X(x)	fprintf(f, "\t%s=%.4g", #x, x);
	X(PctHits);
	X(PctMisses);
#undef X	
	fprintf(f, "\tcmdline=%s", CmdLine.c_str());
	fprintf(f, "\n");
	}

void ShapeSearcher::LogStats()
	{
	g_Misses = g_Queries - g_Hits;
#define X(x)	Log("%10u  %s\n", g_##x, #x);
	X(Hits);
	X(Misses);
	X(PermutedHits);
	X(Miss_LowABCScore);
	X(Miss_LowPalmScore);
	//X(Miss_HighLEFPPM);
	X(Miss_MissingMotif);
	if (opt_calibrate)
		{
		X(TPCount);
		X(TNCount);
		X(FPCount);
		X(FNCount);
		}
#undef X
	}

void ShapeSearcher::SetShapeIndexesABC()
	{
	m_ShapeIndexA = UINT_MAX;
	m_ShapeIndexB = UINT_MAX;
	m_ShapeIndexC = UINT_MAX;
	const vector<string> &Names = m_Shapes->m_Names;
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		if (Names[i] == "A")
			m_ShapeIndexA = i;
		else if (Names[i] == "B")
			m_ShapeIndexB = i;
		else if (Names[i] == "C")
			m_ShapeIndexC = i;
		}
	}

void ShapeSearcher::Init(const Shapes &S)
	{
	m_Shapes = &S;
	m_ShapeCount = S.GetShapeCount();

	SetShapeIndexesABC();
	SetParamOpts();
	}

void ShapeSearcher::SetParamOpts()
	{
#define S(x, y)	if (!optset_##x) { optset_##x = true; opt_##x = (y); }
	if (optset_lowerrors)
		{
		S(minselfscorepp, 0.55);
		S(minselfscorenonpp, 0.55);
		S(minppscore, 0.6);
		S(minpalmscore, 0.69);
		//S(maxlefppm, 1);
		S(searchmfs, "*");
		S(requiremfs, "JABCD");
		}
	else if (optset_calibrate)
		{
		S(minselfscorepp, 0.4);
		S(minselfscorenonpp, 0.4);
		S(minppscore, 0.2);
		S(minpalmscore, 0.2);
		//S(maxlefppm, 10);
		S(searchmfs, "*");
		S(requiremfs, "*");
		}
	else
		{
	// default is -sensitive
		if (!optset_sensitive)
			{
			optset_sensitive = true;
			opt_sensitive = true;
			}
		S(minselfscorepp, 0.50);
		S(minselfscorenonpp, 0.50);
		S(minppscore, 0.55);
		S(minpalmscore, 0.69);
		//S(maxlefppm, 1000);
		S(searchmfs, "*");
		S(requiremfs, "ABC");
		}
#undef S

	m_MinABCScore = opt_minppscore;
	m_MinDomScore = opt_minpalmscore;
	//m_MaxLEFPPM = opt_maxlefppm;
	m_MinSelfScoreABC = opt_minselfscorepp;
	m_MinSelfScoreNonABC = opt_minselfscorenonpp;

	IncludesStrToBools("-searchmfs", opt_searchmfs,  m_SearchShapes);
	IncludesStrToBools("-requiremfs", opt_requiremfs, m_RequireShapes);
	asserta(m_ShapeIndexA < m_ShapeCount);
	asserta(m_ShapeIndexB < m_ShapeCount);
	asserta(m_ShapeIndexC < m_ShapeCount);
	if (!m_SearchShapes[m_ShapeIndexA] ||
		!m_SearchShapes[m_ShapeIndexB] ||
		!m_SearchShapes[m_ShapeIndexC])
		Die("Bad -searchmfs, must include ABC");

	vector<bool> BoolsABC;
	IncludesStrToBools("BoolsABC", "ABC", BoolsABC);
	m_SearchABCOnly = false;
	if (m_SearchShapes == BoolsABC)
		m_SearchABCOnly = true;
	}

void ShapeSearcher::LogParams() const
	{
	Log("Motifs: ");
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		if (i > 0)
			Log(", ");
		Log("%s", GetShapeName(i));
		Log("/%u", GetShapeLength(i));
		}
	Log("\n");

	string SearchStr;
	string RequireStr;
	IncludesBoolsToStr(m_SearchShapes, SearchStr);
	IncludesBoolsToStr(m_RequireShapes, RequireStr);

	Log("%12.12s  Search motifs\n", SearchStr.c_str());
	Log("%12.12s  Require motifs\n", RequireStr.c_str());
	Log("%12c  SearchABCOnly\n", yon(m_SearchABCOnly));

#define Sf(x)	Log("%12.4g  %s\n", m_##x, #x);
	Sf(MinABCScore);
	//Sf(MaxLEFPPM);
	Sf(MinSelfScoreNonABC);
	Sf(MinSelfScoreABC);
	Sf(Sigmas);
#undef Sf

#define Si(x)	Log("%12u  %s\n", m_##x, #x);
	Si(MaxTopHitCount);
	Si(MaxTopHitCountABC);
	Si(ShapeCount);
#undef Si
	}

void ShapeSearcher::SetClass()
	{
	m_Class = ".";
	char Gate = GetGate();
	if (m_ABCScore >= m_MinABCScore)
		{
		switch (Gate)
			{
		case 'D':
		case 'S':
		case 'T':
		case 'G':
		case 'C':
			m_Class = "RdRp+";
			IncStat(g_RdRpPlus);
			break;

		case 'Y':
		case 'F':
			m_Class = "RdRp-";
			IncStat(g_RdRpMinus);
			break;

		default:
			m_Class = "RdRp_like";
			IncStat(g_RdRpLike);
			break;
			}
		}
	}

bool ShapeSearcher::IsHit() const
	{
	if (m_Query == 0)
		return false;
	IncStat(g_Queries);
	if (m_ABCScore < m_MinABCScore)
		{
		IncStat(g_Miss_LowABCScore);
		return false;
		}
	if (m_DomScore < m_MinDomScore)
		{
		IncStat(g_Miss_LowPalmScore);
		return false;
		}
	//if (m_LEFPPM > m_MaxLEFPPM)
	//	{
	//	IncStat(g_Miss_HighLEFPPM);
	//	return false;
	//	}
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		if (m_RequireShapes[i] && m_ShapePosVec[i] == UINT_MAX)
			{
			IncStat(g_Miss_MissingMotif);
			return false;
			}
		}
	if (m_Permuted)
		IncStat(g_PermutedHits);
	IncStat(g_Hits);
	return true;
	}

void ShapeSearcher::SetQuery(const PDBChain &Query)
	{
	ClearSearch();
	m_Query = &Query;

	const uint ShapeCount = GetShapeCount();
	m_ShapePosVec.resize(ShapeCount);
	m_ShapeScores.resize(ShapeCount);
	for (uint i = 0; i < ShapeCount; ++i)
		{
		m_ShapePosVec[i] = UINT_MAX;
		m_ShapeScores[i] = 0;
		}
	}

void ShapeSearcher::GetFoundMotifsStr(string &Str) const
	{
	Str.resize(0);
	if (SIZE(m_ShapePosVec) != m_ShapeCount)
		{
		Str = "*";
		return;
		}
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		bool Found = (m_ShapePosVec[i] != UINT_MAX);
		const string Name = (string) GetShapeName(i);
		if (Found)
			Str += Name;
		else
			{
			Str += "-";
			for (uint j = 1; j < SIZE(Name); ++j)
				Str += ".";
			}
		}

	if (m_Permuted)
		{
		size_t n = Str.find("ABC");
		if (n == string::npos)
			Str += "-permuted";
		else
			{
			Str[n] = 'C';
			Str[n+1] = 'A';
			Str[n+2] = 'B';
			}
		}
	}

char ShapeSearcher::GetGate() const
	{
	if (m_Query == 0)
		return '.';
	if (m_ShapeIndexA == UINT_MAX)
		return '.';
	asserta(m_ShapeIndexA < SIZE(m_ShapePosVec));
	if (m_PosA == UINT_MAX)
		return '.';
	const string &Seq = m_Query->m_Seq;
	uint Pos = m_PosA + g_OffAd + 5;
	if (Pos >= SIZE(Seq))
		return '.';
	char Gate = Seq[Pos];
	return Gate;
	}

void ShapeSearcher::GetSubSeq(uint Pos, uint n, string &Seq) const
	{
	m_Query->GetSubSeq(Pos, n, Seq);
	}

void ShapeSearcher::GetA(string &A) const
	{
	GetShapeSeq(m_ShapeIndexA, A);
	}

void ShapeSearcher::GetB(string &B) const
	{
	GetShapeSeq(m_ShapeIndexB, B);
	}

void ShapeSearcher::GetC(string &C) const
	{
	GetShapeSeq(m_ShapeIndexC, C);
	}

void ShapeSearcher::GetShapeSeq(uint ShapeIndex, string &Seq) const
	{
	Seq = ".";
	if (ShapeIndex == UINT_MAX)
		return;
	uint ShapeLength = GetShapeLength(ShapeIndex);
	asserta(ShapeIndex < SIZE(m_ShapePosVec));
	uint Pos = m_ShapePosVec[ShapeIndex];
	m_Query->GetSubSeq(Pos, ShapeLength, Seq);
	}

// dist = start_j - start_i
void ShapeSearcher::GetDistRange(uint ShapeIndexi, uint ShapeIndexj,
  uint &MinDist, uint &MaxDist) const
	{
	asserta(ShapeIndexj > ShapeIndexi);

	MinDist = 0;
	MaxDist = 0;
	for (uint ShapeIndex = ShapeIndexi; ShapeIndex < ShapeIndexj;
	  ++ShapeIndex)
		{
		MinDist += m_Shapes->m_MinNeighborDists[ShapeIndex];
		MaxDist += m_Shapes->m_MaxNeighborDists[ShapeIndex];
		}
	}

double ShapeSearcher::GetScoreShapes(const vector<uint> &PosVec) const
	{
	double Sum = 0;
	uint Pairs = 0;
	uint FoundCount = 0;
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		uint Posi = PosVec[i];
		if (Posi == UINT_MAX)
			continue;

		for (uint j = 0; j <= i; ++j)
			{
			uint Posj = PosVec[j];
			if (Posj == UINT_MAX)
				continue;

			double Score =
			  GetScoreShapePair(i, j, Posi, Posj);
			Sum += Score;
			++Pairs;
			}
		}
	if (Pairs == 0)
		return 0;
	double Score = Sum/Pairs;
	return Score;
	}

void ShapeSearcher::LogPairwiseScores(const vector<uint> &PosVec) const
	{
	Log("\n");
	double Sum = 0;
	Log(">%s\n", m_Query->m_Label.c_str());

	uint Pairs = 0;
	uint FoundCount = 0;
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		uint Posi = PosVec[i];
		if (Posi == UINT_MAX)
			continue;
		const string Namei = GetShapeName(i);
		Log("%2s[%4u]: ", Namei.c_str(), Posi);
		for (uint j = 0; j <= i; ++j)
			{
			uint Posj = PosVec[j];
			if (Posj == UINT_MAX)
				continue;
			const string Namej = GetShapeName(j);

			double Score =
			  GetScoreShapePair(i, j, Posi, Posj);
			Sum += Score;
			Log(" %s(%.3g)", Namej.c_str(), Score);
			++Pairs;
			}
		Log("\n");
		}
	double Avg = Sum/Pairs;
	double Score = Avg;
	Log(" /%u = %.3g\n", Pairs, Score);
	}

double ShapeSearcher::GetScoreShapePair(uint ShapeIndex1, uint ShapeIndex2,
  uint Pos1, uint Pos2) const
	{
	if (ShapeIndex1 == UINT_MAX || ShapeIndex2 == UINT_MAX)
		return 0;
	if (Pos1 == UINT_MAX || Pos2 == UINT_MAX)
		return 0;
	double Sum = 0;
	uint L1 = GetShapeLength(ShapeIndex1);
	uint L2 = GetShapeLength(ShapeIndex2);
	uint QL = GetQL();
	if (Pos1 + L1 > QL)
		return 0;
	if (Pos2 + L2 > QL)
		return 0;
	uint n = 0;
	for (uint Offset1 = 0; Offset1 < L1; ++Offset1)
		{
		uint StartOffset2 = (ShapeIndex1 == ShapeIndex2 ? Offset1 + 1 : 0);
		for (uint Offset2 = StartOffset2; Offset2 < L2; ++Offset2)
			{
			Sum += GetScoreResiduePair(ShapeIndex1, ShapeIndex2,
			  Pos1, Pos2, Offset1, Offset2);
			++n;
			}
		}
	assert(n > 0);
	double Score = Sum/n;
	return Score;
	}

double ShapeSearcher::GetScoreResiduePair(uint ShapeIndex1, uint ShapeIndex2,
  uint Pos1, uint Pos2, uint Offset1, uint Offset2) const
	{
	bool Diag = (ShapeIndex1 == ShapeIndex2);
	double XS = (Diag ? 1.5 : 2);
	const PDBChain &Query = *m_Query;
	double Observed_d = Query.GetDist(Pos1+Offset1, Pos2+Offset2);
	Observed_d = fabs(Observed_d);
	double Mu = GetMeanDist3(ShapeIndex1, ShapeIndex2, Offset1, Offset2);
	Mu = fabs(Mu);
	double Sigma = GetStdDev3(ShapeIndex1, ShapeIndex2, Offset1, Offset2);
	double y = GetNormal(Mu, XS*Sigma, Observed_d);
	double Max = GetNormal(Mu, XS*Sigma, Mu);
	double Score = y/Max;
	return Score;
	}

double ShapeSearcher::GetSelfScore(uint ShapeIndex, uint Pos) const
	{
	double Score = GetScoreShapePair(ShapeIndex, ShapeIndex, Pos, Pos);
	return Score;
	}

double ShapeSearcher::GetScore(uint ShapeIndex, uint Pos,
  const vector<uint> &PosVec) const
	{
	const uint ShapeCount = GetShapeCount();
	asserta(SIZE(PosVec) == ShapeCount);
	double Sum = 0;
	uint n = 0;
	for (uint ShapeIndex2 = 0; ShapeIndex2 < ShapeCount; ++ShapeIndex2)
		{
		uint Pos2 = PosVec[ShapeIndex2];
		if (Pos2 == UINT_MAX)
			continue;
		++n;
		Sum += GetScoreShapePair(ShapeIndex, ShapeIndex2, Pos, Pos2);
		}
	asserta(n > 0);
	double Score = Sum/n;
	return Score;
	}

void ShapeSearcher::SearchShapeSelfTop(uint ShapeIndex,
  double MinScore, uint MaxHits, vector<uint> &HitPosVec) const
	{
	Die("SearchShapeSelfTop");
#if 0
	HitPosVec.clear();
	uint QL = GetQL();
	uint Length = GetShapeLength(ShapeIndex);
	if (QL < Length + 10)
		return;

	vector<uint> TmpHitPosVec;
	vector<double> TmpHitScores;
	SearchShapeSelf(ShapeIndex, m_MinSelfScore, 0, QL-Length, 0, UINT_MAX,
	  TmpHitPosVec, TmpHitScores);
	const uint N = SIZE(TmpHitPosVec);
	asserta(SIZE(TmpHitScores) == N);
	if (N == 0)
		return;
	vector<uint> Order(N);
	QuickSortOrderDesc(TmpHitScores.data(), N, Order.data());
	uint M = min(N, MaxHits);
	for (uint k = 0; k < M; ++k)
		{
		uint i = Order[k];
		HitPosVec.push_back(TmpHitPosVec[i]);
		}
#endif
	}

double ShapeSearcher::SearchShapeSelf(uint ShapeIndex, double MinScore,
  uint Lo, uint Hi, char Letter, uint LetterOffset,
  vector<uint> &HitPosVec, vector<double> &HitScores) const
	{
	if (Lo == UINT_MAX || Hi == UINT_MAX)
		return 0;
	HitPosVec.clear();
	HitScores.clear();
	uint L = GetShapeLength(ShapeIndex);
	uint QL = GetQL();
	if (Lo + L > QL)
		return 0;
	double TopScore = 0;
	uint Top = Hi;
	if (Top + L > QL)
		{
		if (QL < L)
			return 0;
		Top = QL - L;
		}
	for (uint Pos = Lo; Pos <= Top; ++Pos)
		{
		if (Letter != 0 && LetterOffset != UINT_MAX &&
		  m_Query->m_Seq[Pos+LetterOffset] != Letter)
			continue;
		double Score = GetSelfScore(ShapeIndex, Pos);
		if (Score < MinScore)
			continue;
		if (Score > TopScore)
			TopScore = Score;
		uint HitCount = SIZE(HitPosVec);
		if (HitCount > 0)
			{
			uint LastHit = HitPosVec[HitCount-1];
			if (Pos - LastHit < L)
				{
				double LastScore = HitScores[HitCount-1];
				if (Score > LastScore)
					{
					HitPosVec[HitCount-1] = Pos;
					HitScores[HitCount-1] = Score;
					}
				continue;
				}
			}
		HitPosVec.push_back(Pos);
		HitScores.push_back(Score);
		}
	return TopScore;
	}

void ShapeSearcher::SearchShapeTopHit(uint ShapeIndex,
  const vector<uint> &PosVec, double MinScore, uint Lo, uint Hi,
  char Letter, uint LetterOffset, uint &Pos, double &Score) const
	{
	Pos = UINT_MAX;
	Score = 0;
	if (Lo == UINT_MAX)
		return;
	asserta(Hi != UINT_MAX);

	Pos = UINT_MAX;
	Score = 0;

	vector<uint> HitPosVec;
	vector<double> Scores;
	SearchShape(ShapeIndex, PosVec, MinScore, Lo, Hi,
	  Letter, LetterOffset, HitPosVec, Scores);

	const uint N = SIZE(HitPosVec);
	for (uint i = 0; i < N; ++i)
		{
		if (Scores[i] > Score)
			{
			Pos = HitPosVec[i];
			Score = Scores[i];
			}
		}
	}

void ShapeSearcher::SearchShape(uint ShapeIndex, const vector<uint> &PosVec,
  double MinScore, uint Lo, uint Hi, char Letter, uint LetterOffset,
  vector<uint> &HitPosVec, vector<double> &HitScores) const
	{
	HitPosVec.clear();
	HitScores.clear();
	if (Lo == UINT_MAX)
		return;
	asserta(Hi != UINT_MAX);
	uint L = GetShapeLength(ShapeIndex);
	uint QL = GetQL();
	//asserta(Hi + L <= QL);
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		if (Letter != 0 && LetterOffset != UINT_MAX &&
		  m_Query->m_Seq[Pos+LetterOffset] != Letter)
			continue;
		double Score = GetScore(ShapeIndex, Pos, PosVec);
		if (Score < MinScore)
			continue;
		uint HitCount = SIZE(HitPosVec);
		if (HitCount > 0)
			{
			uint LastHit = HitPosVec[HitCount-1];
			if (Pos - LastHit < L)
				{
				double LastScore = HitScores[HitCount-1];
				if (Score > LastScore)
					{
					HitPosVec[HitCount-1] = Pos;
					HitScores[HitCount-1] = Score;
					}
				continue;
				}
			}
		HitPosVec.push_back(Pos);
		HitScores.push_back(Score);
		}
	}

void ShapeSearcher::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	static bool HdrDone = false;
#pragma omp critical
	{
	if (!HdrDone)
		{
		HdrDone = true;
		fprintf(f, "Label");
		fprintf(f, "\tClass");
		//fprintf(f, "\tLEFPPM");
		fprintf(f, "\tPalm_score");
		fprintf(f, "\tPP_score");
		fprintf(f, "\tMotifs");
		fprintf(f, "\tGate");
		for (uint i = 0; i < m_ShapeCount; ++i)
			{
			const char *ShapeName = GetShapeName(i);
			fprintf(f, "\t%s_pos", ShapeName);
			fprintf(f, "\t%s_seq", ShapeName);
			fprintf(f, "\t%s_score", ShapeName);
			}
		fprintf(f, "\n");
		}

	string FoundMotifsStr;
	GetFoundMotifsStr(FoundMotifsStr);

	char Gate = GetGate();
	fprintf(f, "%s", m_Query->m_Label.c_str());
	fprintf(f, "\t%s", m_Class.c_str());
	//fprintf(f, "\t%.3f", m_LEFPPM);
	fprintf(f, "\t%.3g", m_DomScore);
	fprintf(f, "\t%.3g", m_ABCScore);
	fprintf(f, "\t%s", FoundMotifsStr.c_str());
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

void ShapeSearcher::ToPmlABC(FILE *f) const
	{
	if (f == 0)
		return;

	if (m_PosA == UINT_MAX || m_PosB == UINT_MAX || m_PosC == UINT_MAX)
		return;

	string SeqA, SeqB, SeqC;
	GetShapeSeq(m_ShapeIndexA, SeqA);
	GetShapeSeq(m_ShapeIndexB, SeqB);
	GetShapeSeq(m_ShapeIndexC, SeqC);

	fprintf(f, "color gray70\n");
	fprintf(f, "select pepseq %s\n", SeqA.c_str());
	fprintf(f, "color tv_blue, sele\n");
	fprintf(f, "select pepseq %s\n", SeqB.c_str());
	fprintf(f, "color tv_green, sele\n");
	fprintf(f, "select pepseq %s\n", SeqC.c_str());
	fprintf(f, "color tv_red, sele\n");
	fprintf(f, "deselect\n");
	}

void ShapeSearcher::ToPml(FILE *f, const string &LoadName) const
	{
	if (f == 0)
		return;

	fprintf(f, "reinitialize;\n");
	if (LoadName != "")
		fprintf(f, "load %s;\n", LoadName.c_str());
	fprintf(f, "color gray70;\n");
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		string Seq, Color;
		GetShapeSeq(i, Seq);
		if (Seq == "" || Seq == ".")
			continue;

		const char *Name = GetShapeName(i);
		GetColor(i, Color);
		if (opt_pmlrevmotifs)
			reverse(Seq.begin(), Seq.end());

		fprintf(f, "select Motif_%s, pepseq %s;\n", Name, Seq.c_str());
		fprintf(f, "color %s, Motif_%s;\n", Color.c_str(), Name);
		}
	fprintf(f, "deselect\n");
	if (opt_pml_savepse)
		{
		const string &SaveName = m_Query->m_Label;
		fprintf(f, "save %s.pse\n", SaveName.c_str());
		}
	}

void ShapeSearcher::GetColor(unsigned MotifIndex, string &Color) const
	{
	string Name = (string) GetShapeName(MotifIndex);
	if (Name == "A")
		Color = "tv_blue";
	else if (Name == "B")
		Color = "tv_green";
	else if (Name == "C")
		Color = "tv_red";
	else if (Name == "F1")
		Color = "deepteal";
	else if (Name == "F2")
		Color = "cyan";
	else if (Name == "D")
		Color = "yellow";
	else if (Name == "E")
		Color = "orange";
	else if (Name == "H")
		Color = "salmon";
	else if (Name == "J")
		Color = "magenta";
	}

void ShapeSearcher::GetShapeIndexes(vector<uint> &Indexes) const
	{
	Indexes.resize(m_ShapeCount);
	for (uint i = 0; i < m_ShapeCount; ++i)
		Indexes[i] = i;
	}

void ShapeSearcher::IncludesBoolsToStr(const vector<bool> &Includes,
  string &Str) const
	{
	asserta(SIZE(Includes) == m_ShapeCount);
	Str.clear();
	for (uint ShapeIndex = 0; ShapeIndex < m_ShapeCount; ++ShapeIndex)
		{
		const string ShapeName = (string) GetShapeName(ShapeIndex);
		if (Includes[ShapeIndex])
			Str += ShapeName;
		else
			{
			uint n = SIZE(ShapeName);
			Str += "-";
			for (uint i = 1; i < n; ++i)
				Str += ".";
			}
		}
	}

uint ShapeSearcher::IncludesStrToBools(const string &What,
  const string &Str, vector<bool> &Includes) const
	{
	Includes.resize(m_ShapeCount);
	if (Str == "*")
		{
		for (uint i = 0; i < m_ShapeCount; ++i)
			Includes[i] = true;
		return m_ShapeCount;
		}

	uint FoundCount = 0;
	for (uint ShapeIndex = 0; ShapeIndex < m_ShapeCount; ++ShapeIndex)
		{
		const string ShapeName = (string) GetShapeName(ShapeIndex);
		size_t n = Str.find(ShapeName);
		bool Found = (n >= 0 && n < SIZE(Str));
		if (Found)
			++FoundCount;
		Includes[ShapeIndex] = Found;
		}
	return FoundCount;
	}

void ShapeSearcher::BoolsToIndexVec1(const vector<bool> &Includes,
  vector<uint> &ShapeIndexes) const
	{
	ShapeIndexes.clear();
	asserta(SIZE(Includes) == m_ShapeCount);
	for (uint i = 0; i < m_ShapeCount; ++i)
		if (Includes[i])
			ShapeIndexes.push_back(i);
	}

void ShapeSearcher::BoolsToIndexVec2(const vector<bool> &Includes,
  vector<uint> &ShapeIndexes) const
	{
	ShapeIndexes.clear();
	asserta(SIZE(Includes) == m_ShapeCount);
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		if (Includes[i])
			ShapeIndexes.push_back(i);
		else
			ShapeIndexes.push_back(UINT_MAX);
		}
	asserta(SIZE(ShapeIndexes) == m_ShapeCount);
	}

#include "myutils.h"
#include "scop40bit.h"
#include "sort.h"
#include "outputfiles.h"

static const uint NRBINS = 32;

static uint ScoreToBin(double Score, double MinScore, double MaxScore)
	{
	if (Score < MinScore)
		return 0;
	if (Score >= MaxScore)
		return NRBINS-1;
	double r = (Score - MinScore)/(MaxScore - MinScore);
	asserta(r >= 0 && r <= 1);
	uint Bin = uint(r*(NRBINS-1));
	asserta(Bin < NRBINS);
	return Bin;
	}

double ConvertEvalueToScore(double Score)
	{
	asserta(Score >= 0);
	if (Score <= 1e-99)
		return 99;
	Score = -log10(Score);
	return Score;
	}

static void ScoresToTsv(FILE *f, const vector<double> &Scores,
  const vector<bool> &TFs)
	{
	const uint N = SIZE(Scores);
	if (N == 0)
		Die("N==0");
	asserta(SIZE(TFs) == N);
	double MinScore = DBL_MAX;
	double MaxScore = DBL_MAX;
	for (uint i = 0; i < N; ++i)
		{
		double Score = Scores[i];
		if (i == 0)
			{
			MinScore = Score;
			MaxScore = Score;
			}
		else
			{
			MinScore = min(Score, MinScore);
			MaxScore = max(Score, MaxScore);
			}
		}
	vector<uint> Counts(NRBINS);
	vector<uint> TCounts(NRBINS);
	double NT = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint Bin = ScoreToBin(Scores[i], MinScore, MaxScore);
		++Counts[Bin];
		if (TFs[i])
			{
			++NT;
			++TCounts[Bin];
			}
		}

	uint Sum = 0;
	double SumFract = 0;
	for (uint Bin = 0; Bin < NRBINS; ++Bin)
		{
		uint n = Counts[Bin];
		uint tn = TCounts[Bin];
		double Fract = double(n)/N;
		double TFract = double(tn)/NT;
		Sum += n;
		SumFract += Fract;
		double BinLo = Bin*(MaxScore - MinScore)/NRBINS;
		double BinHi = (Bin + 1)*(MaxScore - MinScore)/NRBINS;
		double BinMid = (BinLo + BinHi)/2;
		fprintf(g_ftsv, "%u\t%.3g\t%u\t%.4g\t%4.g\n",
		  Bin, BinMid, n, Fract, TFract);
		}
	asserta(Sum == N);
	asserta(feq(SumFract, 1.0));
	}

void cmd_scop40bit_evd()
	{
	asserta(optset_ref); // name of domain
	asserta(optset_tsv);
	const string &DomName = opt_ref;
	const string &FN = g_Arg1;
	SCOP40Bit SB;
	SB.ReadDomInfo();
	SB.ReadHits_Bin(FN);
	const uint DomIdx = SB.GetDomIdx(DomName);
	const uint HitCount = SB.GetHitCount();
	vector<double> Scores;
	vector<double> TScores;
	bool ScoresAreEvalues = SB.m_ScoresAreEvalues;
	double MinScore = 0;
	double MaxScore = 0;
	vector<bool> TFs;
	for (uint i = 0; i < HitCount; ++i)
		{
		uint DomIdx1 = SB.m_DomIdx1s[i];
		uint DomIdx2 = SB.m_DomIdx2s[i];
		if (DomIdx1 == DomIdx || DomIdx2 == DomIdx)
			{
			bool TF = SB.IsT(DomIdx1, DomIdx2);
			double Score = SB.m_Scores[i];
			if (ScoresAreEvalues)
				Score = ConvertEvalueToScore(Score);
			Scores.push_back(Score);
			TFs.push_back(TF);
			}
		}

	sort(Scores.begin(), Scores.end());
	const uint N = SIZE(Scores);
	const double OutlierFract = 0.05;
	const uint Lo = uint(N*OutlierFract);
	const uint Hi = uint(N*(1.0 - OutlierFract));
	vector<double> TrimmedScores;
	vector<bool> TrimmedTFs;
	for (uint i = Lo ; i < Hi; ++i)
		{
		TrimmedScores.push_back(Scores[i]);
		TrimmedTFs.push_back(TFs[i]);
		}
	ScoresToTsv(g_ftsv, TrimmedScores, TrimmedTFs);
	}

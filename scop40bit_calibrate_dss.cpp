#include "myutils.h"
#include "scop40bit.h"
#include "sort.h"
#include "outputfiles.h"

void cmd_scop40bit_calibrate_dss()
	{
	asserta(optset_tsv);
	const string &FN = g_Arg1;
	SCOP40Bit SB;
	SB.ReadDomInfo();
	SB.ReadHits_Bin(FN);
	const uint DomCount = SB.GetDomCount();
	const uint HitCount = SB.GetHitCount();
	vector<float> Scores;
	vector<float> TScores;
	float MinScore = 0;
	float MaxScore = 0;
	vector<bool> TFs;
	vector<float> SumScore2(DomCount);
	vector<uint> NegativeScoreCounts(DomCount);
	vector<float> DomToScoreFirstFP;
	DomToScoreFirstFP.resize(HitCount, 0);
	for (uint i = 0; i < HitCount; ++i)
		{
		uint DomIdx1 = SB.m_DomIdx1s[i];
		uint DomIdx2 = SB.m_DomIdx2s[i];
		float Score = SB.m_Scores[i];
		uint Fam1 = SB.m_DomIdxToFamIdx[DomIdx1];
		uint Fam2 = SB.m_DomIdxToFamIdx[DomIdx2];
		if (Fam1 != Fam2)
			DomToScoreFirstFP[DomIdx1] = max(DomToScoreFirstFP[DomIdx1], Score);
		if (Score < 0)
			{
			++NegativeScoreCounts[DomIdx1];
			SumScore2[DomIdx1] += Score*Score;
			}
		}

	float Sumx = 0;
	float Sumx2 = 0;
	float Sumy = 0;
	float Sumy2 = 0;
	float Sumxy = 0;
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		const string &DomName = SB.m_Doms[DomIdx];
		uint n = NegativeScoreCounts[DomIdx];
		float Sum2 = SumScore2[DomIdx];
		float ScoreFirstFP = DomToScoreFirstFP[DomIdx];
		ScoreFirstFP = powf(ScoreFirstFP, 0.2);
		float Var = 0;
		if (n > 0)
			Var = sqrtf(Sum2/n);
		Sumx += Var;
		Sumx2 += Var*Var;
		Sumy += ScoreFirstFP;
		Sumy2 += ScoreFirstFP*ScoreFirstFP;
		Sumxy += Var*ScoreFirstFP;
		fprintf(g_ftsv, "%u\t%s\t%u\t%.3g\t%.3g\n",
		  DomIdx, DomName.c_str(), n, Var, ScoreFirstFP);
		}

	float top = Sumxy - (Sumx*Sumy)/DomCount;
	float bottomL = Sumx2 - (Sumx*Sumx)/DomCount;
	float bottomR = Sumy2 - (Sumy*Sumy)/DomCount;
	float Correl = top/sqrtf(bottomL*bottomR);
	ProgressLog("Correl = %.3g\n", Correl);
	}

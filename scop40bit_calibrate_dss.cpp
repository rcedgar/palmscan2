#include "myutils.h"
#include "scop40bit.h"
#include "sort.h"
#include "outputfiles.h"

void cmd_scop40bit_calibrate_dss()
	{
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

	vector<float> Pows;
	//Pows.push_back(0.1f);
	//Pows.push_back(0.2f);
	//Pows.push_back(0.4f);
	Pows.push_back(0.5f);
	//Pows.push_back(0.6f);
	//Pows.push_back(0.7f);
	//Pows.push_back(0.8f);
	//Pows.push_back(1.0f);
	//Pows.push_back(1.5f);
	//Pows.push_back(2.0f);
	//Pows.push_back(3.0f);
	vector<float> Vars;
	for (vector<float>::const_iterator iter = Pows.begin();
	  iter != Pows.end(); ++iter)
		{
		float Pow = *iter;

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
			float Var = 0;
			if (n > 0)
				Var = sqrtf(Sum2/n);
			Var = powf(Var, Pow);
			Vars.push_back(Var);
			Sumx += Var;
			Sumx2 += Var*Var;
			Sumy += ScoreFirstFP;
			Sumy2 += ScoreFirstFP*ScoreFirstFP;
			Sumxy += Var*ScoreFirstFP;
			Pf(g_ftsv, "%u\t%s\t%u\t%.3g\t%.3g\n",
			  DomIdx, DomName.c_str(), n, Var, ScoreFirstFP);
			}

		float top = Sumxy - (Sumx*Sumy)/DomCount;
		float bottomL = Sumx2 - (Sumx*Sumx)/DomCount;
		float bottomR = Sumy2 - (Sumy*Sumy)/DomCount;
		float Correl = top/sqrtf(bottomL*bottomR);
		ProgressLog("Pow = %.2f correl = %.4f\n", Pow, Correl);
		}

	if (optset_output)
		{
		for (uint i = 0; i < HitCount; ++i)
			{
			uint DomIdx1 = SB.m_DomIdx1s[i];
			float Score = SB.m_Scores[i];
			float Var = Vars[DomIdx1];
			SB.m_Scores[i] = Score/Var;
			}

		SB.WriteHits_Bin(opt_output);
		}
	}

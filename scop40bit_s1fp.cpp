#include "myutils.h"
#include "scop40bit.h"
#include "sort.h"

void cmd_scop40bit_s1fp()
	{
	const string &FN = g_Arg1;
	SCOP40Bit SB;
	SB.ReadDomInfo();
	SB.ReadHits_Bin(FN);
	uint NT, NF;
	SB.CalcNXs(NT, NF);
	ProgressLog("NT %u (%s)\n", NT, IntToStr(NT));
	ProgressLog("NF %u (%s)\n", NF, IntToStr(NF));

	const uint DomCount = SB.GetDomCount();
	const uint HitCount = SB.GetHitCount();
	vector<float> DomToScoreFirstFP;
	if (SB.m_ScoresAreEvalues)
		DomToScoreFirstFP.resize(HitCount, FLT_MAX);
	else
		DomToScoreFirstFP.resize(HitCount, 0);

	for (uint i = 0; i < HitCount; ++i)
		{
		uint Dom1 = SB.m_DomIdx1s[i];
		uint Dom2 = SB.m_DomIdx2s[i];
		float Score = SB.m_Scores[i];
		uint Fam1 = SB.m_DomIdxToFamIdx[Dom1];
		uint Fam2 = SB.m_DomIdxToFamIdx[Dom2];
		if (Fam1 != Fam2)
			{
			if (SB.m_ScoresAreEvalues)
				DomToScoreFirstFP[Dom1] = min(DomToScoreFirstFP[Dom1], Score);
			else
				DomToScoreFirstFP[Dom1] = max(DomToScoreFirstFP[Dom1], Score);
			}
		}

	uint GoodCount = 0;
	for (uint i = 0; i < HitCount; ++i)
		{
		uint Dom1 = SB.m_DomIdx1s[i];
		uint Dom2 = SB.m_DomIdx2s[i];
		float Score = SB.m_Scores[i];
		uint Fam1 = SB.m_DomIdxToFamIdx[Dom1];
		uint Fam2 = SB.m_DomIdxToFamIdx[Dom2];
		if (Fam1 == Fam2)
			{
			if (SB.m_ScoresAreEvalues)
				{
				if (Score < DomToScoreFirstFP[Dom1])
					++GoodCount;
				}
			else
				{
				if (Score > DomToScoreFirstFP[Dom1])
					++GoodCount;
				}
			}
		}
	ProgressLog("Sens1FP %u (%s) %s\n",
	  GoodCount, IntToStr(GoodCount), g_Arg1.c_str());
	}

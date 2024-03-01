#include "myutils.h"
#include "seqdb.h"
#include "xprof.h"
#include "pdbchain.h"
#include "alpha.h"
#include "quarts.h"

/***
palmscan2.exe \
  -xfeatures_3di d:/int/scop40/out/domains_scop.cal \
  -input d:/a/res/dave_grant/scop40/scop40.3di \
  -log -xfeatures_3di.log 

[A]   Ang_m2_p2 Min=0, LoQ=23.6, Med=43.4, HiQ=70.7, Max=161, Avg=49.1, StdDev=31.4
[C]   Ang_m2_p2 Min=0, LoQ=84.4, Med=107, HiQ=112, Max=177, Avg=97.3, StdDev=22.9
[D]   Ang_m2_p2 Min=0, LoQ=19, Med=38.4, HiQ=63.8, Max=174, Avg=42.6, StdDev=30.5
[E]   Ang_m2_p2 Min=0, LoQ=11.4, Med=20.7, HiQ=41.4, Max=162, Avg=28.7, StdDev=23
...

Avg id 21.6%
***/

void cmd_xfeatures_3di()
	{
	const string &InputFileName = opt_xfeatures_3di;
	vector<PDBChain *> Chains;
	ReadChains(InputFileName, Chains);
	const uint ChainCount = SIZE(Chains);

	XProf::InitScoreTable();

	SeqDB Fa3di;
	Fa3di.FromFasta(opt_input);
	Fa3di.SetLabelToIndex();

	XProf XP;
	vector<vector<vector<double> > > Values(20);
	for (uint Letter = 0; Letter < 20; ++Letter)
		Values[Letter].resize(XFEATS);

	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Pass 1");
		const PDBChain &Chain = *Chains[ChainIndex];
		XP.Init(Chain);
		const string &Label = Chain.m_Label;
		string Seq3di;
		bool Ok = Fa3di.GetSeqByLabel(Label, Seq3di, false);
		if (!Ok)
			continue;
		const uint L = SIZE(Seq3di);
		if (Chain.GetSeqLength() != L)
			continue;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			char c = Seq3di[Pos];
			uint Letter = g_CharToLetterAmino[c];
			if (Letter >= 20)
				continue;
			for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
				{
				uint notused_iValue;
				double Value;
				XP.GetFeature(FeatureIndex, Pos, Value, notused_iValue);
				Values[Letter][FeatureIndex].push_back(Value);
				}
			}
		}

	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		{
		Log("\n");
		for (uint Letter = 0; Letter < 20; ++Letter)
			{
			char c = g_LetterToCharAmino[Letter];
			const vector<double> &v = Values[Letter][FeatureIndex];
			QuartsDouble QD;
			GetQuartsDouble(v, QD);
			Log("[%c]  %10.10s ", c, XProf::GetFeatureName(FeatureIndex));
			QD.LogMe();
			}
		}

	Log("static double g_3di_Medians[20][%u] = {\n", XFEATS);
	Log("//       ");
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		Log("  %8.8s", XProf::GetFeatureName(FeatureIndex));
	Log("\n");
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		char c = g_LetterToCharAmino[Letter];
		Log("/* %c */ { ", c);
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			const vector<double> &v = Values[Letter][FeatureIndex];
			QuartsDouble QD;
			GetQuartsDouble(v, QD);
			Log(" %8.4g", QD.Med);
			if (FeatureIndex+1 != XFEATS)
				Log(",");
			}
		Log("}, // %c\n", c);
		}
	Log("};\n");

	double SumPctId = 0;
	uint NumPctId = 0;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Pass 2");
		const PDBChain &Chain = *Chains[ChainIndex];
		XP.Init(Chain);
		const string &Label = Chain.m_Label;
		string Seq3di;
		bool Ok = Fa3di.GetSeqByLabel(Label, Seq3di, false);
		if (!Ok)
			continue;
		const uint L = SIZE(Seq3di);
		if (Chain.GetSeqLength() != L)
			continue;
		string MySeq;
		string Annot;
		uint Ids = 0;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			char c = Seq3di[Pos];
			uint Letter = g_CharToLetterAmino[c];
			if (Letter >= 20)
				continue;

			vector<double> FeatureValues;
			for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
				{
				uint notused_iValue;
				double Value;
				XP.GetFeature(FeatureIndex, Pos, Value, notused_iValue);
				FeatureValues.push_back(Value);
				}
			uint MyLetter3di = XProf::Get3di(FeatureValues);
			char MyChar3di = g_LetterToCharAmino[MyLetter3di];
			if (MyChar3di == c)
				{
				++Ids;
				Annot += "|";
				}
			else
				Annot += " ";
			MySeq += MyChar3di;
			}
		double PctId = GetPct(Ids, L);
		++NumPctId;
		SumPctId += PctId;
		Log(">%s (%.1f%%)\n", Chain.m_Label.c_str(), PctId);
		Log("%s\n", Seq3di.c_str());
		Log("%s\n", MySeq.c_str());
		}
	ProgressLog("Avg id %.1f%%\n", SumPctId/NumPctId);
	}

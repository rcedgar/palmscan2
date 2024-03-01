#include "myutils.h"
#include "xprof.h"
#include "pdbchain.h"
#include "outputfiles.h"

/***
Convert .cal file to tsv with one line per position =
  XProf feature vector values.

palmscan2 \
  -xfeatures d:/int/scop40/out/domains_scop.cal \
  -tsv d:/int/scop40/out/xfeatures.tsv

# head d:/int/scop40/out/xfeatures.tsv | columns.py
aa  Ang_m2_p2  Ang_m3_p3  ED_p4  ED_m4  NU  ND
 A          0          0   6.07      0   0   0
 Y          0          0   5.92      0   7   7
 I        113          0   6.38      0   0  11
 A        115       34.1   6.33      0   2  12
 K        109         28   6.23   6.07  15   9

***/

void cmd_xfeatures()
	{
	const string &InputFileName = opt_xfeatures;
	vector<PDBChain *> Chains;
	ReadChains(InputFileName, Chains);
	const uint ChainCount = SIZE(Chains);
	asserta(g_ftsv);

	fprintf(g_ftsv, "aa");
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		fprintf(g_ftsv, "\t%s", XProf::GetFeatureName(FeatureIndex));
	fprintf(g_ftsv, "\n");

	XProf XP;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Running");
		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Seq = Chain.m_Seq;
		XP.Init(Chain);
		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			fprintf(g_ftsv, "%c", Seq[Pos]);
			for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
				{
				double Value;
				uint notused_iValue;
				XP.GetFeature(FeatureIndex, Pos,
				  Value, notused_iValue);
				fprintf(g_ftsv, "\t%.3g", Value);
				}
			fprintf(g_ftsv, "\n");
			}
		}
	}

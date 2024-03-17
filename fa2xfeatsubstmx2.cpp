#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
#include "outputfiles.h"

void GetScopDomFromLabel(const string &Label, string &Dom);

/***
Optimize binning of 3D feature by expected score.

Input
=====
(*) cal for chains d:/int/scop40/out/domains_scop.cal
# head /d/int/scop40/out/domains_scop.cal
>d12asa_/d.104.1.1
A       12.501  39.048  28.539
Y       15.552  39.410  26.282

(*) Pair-wise alignments of scop40 domains, created by
  palmscan2 -tm_scop d:/int/scop40/out/domains_scop.cal -output scop40.fa2 -tsv scop40.tsv -log tm.log -threads 10

# head /d/a/res/dave_grant/scop40/scop40.fa2

 Dom     Fam     TM     %id
 vvvvvvv vvvvvvv vvvvvv vvvv
>d1a04a1/a.4.6.2/0.7293/31.2
ERDVNQLTPRERDILKLIAQ-GLPNKMIARRLDITESTVKVHVKHMLKKMKLKSRVEAAVWVHQERIF--
>d1fsea_/a.4.6.2
SKP-L-LTKREREVFELL-VQDKTTKEIASELFISEKTVRNHISNAMQKLGVKGRSQAVVELLRMGELEL

>d1914a1/d.49.1.1/0.6033/21.1
--FQTWEEFSRAAEKLYLADP--MKVRVVLKYRHVD-------GNLCIKVTDDLVCLVYRTDQAQDVKKIEK-FHSQLMR-LMVAKESRNV-
>d1914a2/d.49.1.1
MVLLESEQFLTELTRLFQKCRSSGSVFITLKKY--DEGLEPAENKCLLRATDGKRKISTVV-SSKEVNKFQMAY-SNLLRANMDGLKK---R
***/

void cmd_fa2xfeatsubstmx2()
	{
	asserta(optset_feature);
	const uint FeatureIndex = opt_feature;
	const char *FeatureName = DSS::GetFeatureName(FeatureIndex);
	ProgressLog("Feature=%s\n", DSS::GetFeatureName(FeatureIndex));

	SeqDB Input;
	Input.FromFasta(opt_fa2xfeatsubstmx2, true);

	double MinTM = 0.6;
	double MaxTM = 0.8;
	if (optset_mintm)
		MinTM = opt_mintm;
	if (optset_maxtm)
		MaxTM = opt_maxtm;

	vector<PDBChain *> Chains;
	ReadChains(opt_train_cal, Chains);
	const uint ChainCount = SIZE(Chains);
	map<string, uint> DomToChainIndex;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const string &Label = Chains[ChainIndex]->m_Label;
		vector<string> Fields;
		string Dom;
		GetScopDomFromLabel(Label, Dom);
		DomToChainIndex[Dom] = ChainIndex;
		}

	const uint SeqCount = Input.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;
	DSS QX;
	DSS RX;
	uint AlphaSize = DSS::GetAlphaSize(FeatureIndex);
	LogOdds LO;
	LO.Init(AlphaSize);
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Processing");
		const string &QLabel = Input.GetLabel(2*PairIndex);
		const string &RLabel = Input.GetLabel(2*PairIndex+1);
		vector<string> Fields;
		Split(QLabel, Fields, '/');
		asserta(SIZE(Fields) == 4);
		const string &QDom = Fields[0];
		const string &Fam = Fields[1];
		const string sTM = Fields[2];
		const string sPctId = Fields[3];
		double TM = StrToFloat(sTM);
		if (TM < MinTM || TM > MaxTM)
			continue;

		string RDom;
		GetScopDomFromLabel(RLabel, RDom);
		uint QChainIndex = DomToChainIndex[QDom];
		uint RChainIndex = DomToChainIndex[RDom];
		const PDBChain &QChain = *Chains[QChainIndex];
		const PDBChain &RChain = *Chains[RChainIndex];
		uint QL = QChain.GetSeqLength();
		uint RL = RChain.GetSeqLength();
		QX.Init(QChain);
		RX.Init(RChain);
		const string &QRow = Input.GetSeq(2*PairIndex);
		const string &RRow = Input.GetSeq(2*PairIndex+1);
		const uint ColCount = SIZE(QRow);
		asserta(SIZE(QRow) == SIZE(RRow));

		uint QPos = 0;
		uint RPos = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char q = QRow[Col];
			char r = RRow[Col];
			if (!isgap(q) && !isgap(r))
				{
				uint ValueQ = QX.GetFeature(FeatureIndex, QPos);
				uint ValueR = RX.GetFeature(FeatureIndex, RPos);
				LO.AddPair(ValueQ, ValueR);
				}
			if (!isgap(q))
				++QPos;
			if (!isgap(r))
				++RPos;
			}
		}

	vector<vector<double> > ScoreMx;
	double ExpectedScore = LO.GetLogOddsMx(ScoreMx);
	ProgressLog("Expected score [%2u] %7.3g %s\n",
	  AlphaSize,
	  ExpectedScore,
	  FeatureName);
	FILE *f = CreateStdioFile(opt_output);
	vector<double> Freqs;
	vector<vector<double> > FreqMx;
	LO.GetFreqs(Freqs);
	LO.GetFreqMx(FreqMx);
	LO.MxToSrc(f, string("Score_") + string(FeatureName), ScoreMx);
	LO.MxToSrc(f, string("Freqs_") + string(FeatureName), FreqMx);
	LO.VecToSrc(f, FeatureName, Freqs);
	CloseStdioFile(f);
	}

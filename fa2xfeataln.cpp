#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "xprof.h"
#include "outputfiles.h"

/***
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

Output
======
TODO
***/

void GetScopDomFromLabel(const string &Label, string &Dom)
	{
	vector<string> Fields;
	Split(Label, Fields, '/');
	asserta(SIZE(Fields) >= 2);
	Dom = Fields[0];
	}

void cmd_fa2xfeataln()
	{
	SeqDB Input;
	Input.FromFasta(opt_fa2xfeataln, true);

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

	asserta(optset_mintm && optset_maxtm);
	double MinTM = opt_mintm;
	double MaxTM = opt_maxtm;

	const uint SeqCount = Input.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;
	XProf QX;
	XProf RX;
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
		QX.Init(QChain);
		RX.Init(RChain);
		const string &QRow = Input.GetSeq(2*PairIndex);
		const string &RRow = Input.GetSeq(2*PairIndex+1);
		const uint ColCount = SIZE(QRow);
		asserta(SIZE(QRow) == SIZE(RRow));

		vector<string> QFeatRows(6);
		vector<string> RFeatRows(6);
		uint QPos = 0;
		uint RPos = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char q = QRow[Col];
			char r = RRow[Col];
			vector<uint> qiv;
			vector<uint> riv;
			for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
				{
				if (isgap(q))
					QFeatRows[FeatureIndex].push_back('-');
				else
					{
					double Value;
					uint iValue;
					QX.GetFeature(FeatureIndex, QPos, Value, iValue);
					QFeatRows[FeatureIndex].push_back('a' + iValue);
					qiv.push_back(iValue);
					}

				if (isgap(r))
					RFeatRows[FeatureIndex].push_back('-');
				else
					{
					double Value;
					uint iValue;
					RX.GetFeature(FeatureIndex, RPos, Value, iValue);
					RFeatRows[FeatureIndex].push_back('a' + iValue);
					riv.push_back(iValue);
					}
				}
			if (!isgap(q) && !isgap(r))
				{
				asserta(SIZE(qiv) == XFEATS);
				asserta(SIZE(riv) == XFEATS);
				if (g_ftsv)
					{
					fprintf(g_ftsv, "%c", q);
					for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
						fprintf(g_ftsv, "%d", qiv[FeatureIndex]);
					fprintf(g_ftsv, "\t");
					fprintf(g_ftsv, "%c", r);
					for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
						fprintf(g_ftsv, "%d", riv[FeatureIndex]);
					fprintf(g_ftsv, "\n");
					}
				}

			if (!isgap(q))
				++QPos;
			if (!isgap(r))
				++RPos;
			}

		Log("\n____");
		for (uint i = 0; i < SIZE(QRow); ++i)
			Log("_");
		Log("\n");
		Log(">%s\n", QLabel.c_str());
		Log(">%s\n", RLabel.c_str());
		Log("\n");
		Log("Qa  %s\n", QRow.c_str());
		Log("Ra  %s\n", RRow.c_str());
		Log("\n");
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			Log("Q%u  ", FeatureIndex);
			Log("%s\n", QFeatRows[FeatureIndex].c_str());
			Log("R%u  ", FeatureIndex);
			Log("%s\n", RFeatRows[FeatureIndex].c_str());
			Log("\n");
			}
		}
	}

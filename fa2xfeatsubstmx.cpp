#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
//#include "xprof.h"
#include "dss.h"
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

double ValueVecsToSubstMx(const vector<uint> &Values1,
  const vector<uint> &Values2, uint &K)
	{
	const uint PairCount = SIZE(Values1);
	asserta(SIZE(Values2) == PairCount);

	K = 0;
	for (uint i = 0; i < PairCount; ++i)
		{
		if (Values1[i] != UINT_MAX)
			K = max(K, Values1[i]);
		if (Values2[i] != UINT_MAX)
			K = max(K, Values2[i]);
		}
	asserta(K <= 20);
	K += 1;
	asserta(K <= 20);

	vector<uint> CountVec(K);
	vector<vector<uint> > CountMx(K);
	for (uint i = 0; i < K; ++i)
		CountMx[i].resize(K);

	uint LetterPairCount = 0;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		uint iq = Values1[PairIndex];
		uint ir = Values2[PairIndex];
		if (iq < K && ir < K)
			{
			LetterPairCount += 2;
			CountVec[iq] += 1;
			CountVec[ir] += 1;
			CountMx[iq][ir] += 1;
			CountMx[ir][iq] += 1;
			}
		}

	double SumFreq = 0;
	vector<double> Freqs;
	for (uint i = 0; i < K; ++i)
		{
		char c = g_LetterToCharAmino[i];
		double Freq = double(CountVec[i])/(LetterPairCount);
		SumFreq += Freq;
		Log("%c %.4f\n", c, Freq);
		Freqs.push_back(Freq);
		}
	asserta(feq(SumFreq, 1.0));

	FILE *f = CreateStdioFile(opt_output);
	if (f == 0)
		{
		Warning("-output not set");
		return 0;
		}

	double SumFreq2 = 0;
	double H = 0;
// RH=Relative entropy
//   =mean expected score per residue pair
//   =Kullback–Leibler divergence
// Altschul, SF (1991) "Amino acid substitution matrices from an
//   information theoretic perspective." JMB 219(3): 555-565.
	fprintf(f, "aa\tfreq");
	for (uint i = 0; i < K; ++i)
		{
		char c = g_LetterToCharAmino[i];
		fprintf(f, "\t%c", c);
		}
	fprintf(f, "\n");

	double ES = 0;
	for (uint i = 0; i < K; ++i)
		{
		char c = g_LetterToCharAmino[i];
		double Freq = double(CountVec[i])/(2*LetterPairCount);
		SumFreq += Freq;
		fprintf(f, "%c", c);
		fprintf(f, "\t%.4f", Freq);
		for (uint j = 0; j < K; ++j)
			{
			uint n = CountMx[i][j];
			double ObsFreq = double(n)/double(LetterPairCount);
			SumFreq2 += ObsFreq;
			double ExpFreq = double(Freqs[i]*Freqs[j]);
			double Ratio = ObsFreq/ExpFreq;
			double Score = log(Ratio);
			ES += ObsFreq*Score;
			fprintf(f, "\t%.4f", Score);
			}
		fprintf(f, "\n");
		}
	asserta(feq(SumFreq2, 1.0));
	CloseStdioFile(f);

	Log("static double g_AminoSubstMx[K][K] = {\n");
	Log("//      ");
	for (uint i = 0; i < K; ++i)
		{
		char c = g_LetterToCharAmino[i];
		Log("        %c", c);
		}
	Log("\n");
	for (uint i = 0; i < K; ++i)
		{
		char c = g_LetterToCharAmino[i];
		Log("/* %c */ {", c);
		for (uint j = 0; j < K; ++j)
			{
			uint n = CountMx[i][j];
			double ObsFreq = double(n)/double(LetterPairCount);
			SumFreq2 += ObsFreq;
			double ExpFreq = double(Freqs[i]*Freqs[j]);
			double Ratio = ObsFreq/ExpFreq;
			double Score = log(Ratio);
			Log(" %7.4f", Score);
			if (j != 19)
				Log(",");
			}
		Log("}, // %c\n", c);
		}
	Log("};\n");

	return ES;
	}

void cmd_fa2xfeatsubstmx()
	{
	asserta(optset_feature);
	const uint FeatureIndex = opt_feature;
	ProgressLog("Feature=%s\n", DSS::GetFeatureName(FeatureIndex));

	SeqDB Input;
	Input.FromFasta(opt_fa2xfeatsubstmx, true);

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
	//XProf QX;
	//XProf RX;
	DSS QX;
	DSS RX;
	vector<uint> Values1;
	vector<uint> Values2;
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

		vector<string> QFeatRows(6);
		vector<string> RFeatRows(6);
		uint QPos = 0;
		uint RPos = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char q = QRow[Col];
			char r = RRow[Col];
			uint ValueQ = UINT_MAX;
			uint ValueR = UINT_MAX;
			if (isgap(q))
				QFeatRows[FeatureIndex].push_back('-');
			else
				{
				asserta(QPos < QL);
				ValueQ = QX.GetFeature(FeatureIndex, QPos);
				QFeatRows[FeatureIndex].push_back('a' + ValueQ);
				}

			if (isgap(r))
				RFeatRows[FeatureIndex].push_back('-');
			else
				{
				asserta(RPos < RL);
				ValueR = RX.GetFeature(FeatureIndex, RPos);
				RFeatRows[FeatureIndex].push_back('a' + ValueR);
				}
			if (!isgap(q) && !isgap(r))
				{
				Values1.push_back(ValueQ);
				Values2.push_back(ValueR);
				if (QPos >= 3 && RPos >= 3 && QPos + 3 <= QL && RPos + 3 <= RL)
					{
					if (g_ftsv)
						{
						fprintf(g_ftsv, "%c", q);
						fprintf(g_ftsv, "\t%u", ValueQ);
						fprintf(g_ftsv, "\t");
						fprintf(g_ftsv, "%c", r);
						fprintf(g_ftsv, "\t%u", ValueR);
						fprintf(g_ftsv, "\t%s\t%u", QDom.c_str(), QPos);
						fprintf(g_ftsv, "\t%s\t%u", RDom.c_str(), RPos);
						fprintf(g_ftsv, "\n");
						}
					}
				}

			if (!isgap(q))
				++QPos;
			if (!isgap(r))
				++RPos;
			}
		}
	uint K;
	double ES = ValueVecsToSubstMx(Values1, Values2, K);
	ProgressLog("Expected score %s [%u] %.3g\n",
	  DSS::GetFeatureName(FeatureIndex), K, ES);
	}

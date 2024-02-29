#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"

/***
Input
=====
Pair-wise alignments of scop40 domains, created by
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
Log-odds score matrix
# head /d/a/res/dave_grant/scop40/substmx.0.4.0.6.tsv
aa      freq    A       C       D       E       F       G       H      
A       0.0415  0.5330  -0.1004 -0.1395 -0.0768 -0.0880 -0.0375 -0.1575
C       0.0076  -0.1004 2.9846  -0.6905 -0.6766 -0.2257 -0.4758 0.0235
D       0.0269  -0.1395 -0.6905 0.9609  0.3571  -0.4771 0.0191  0.0140

***/

void cmd_fa2substmx()
	{
	SeqDB Input;
	Input.FromFasta(opt_fa2substmx, true);

	asserta(optset_mintm && optset_maxtm);
	double MinTM = opt_mintm;
	double MaxTM = opt_maxtm;

	vector<uint> CountVec(20);
	vector<vector<uint> > CountMx(20);
	for (uint i = 0; i < 20; ++i)
		CountMx[i].resize(20);

	const uint SeqCount = Input.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Processing");
		const string &QRow = Input.GetSeq(2*PairIndex);
		const string &RRow = Input.GetSeq(2*PairIndex+1);
		const string &QLabel = Input.GetLabel(2*PairIndex);
		const string &RLabel = Input.GetLabel(2*PairIndex+1);
		const uint ColCount = SIZE(QRow);
		asserta(SIZE(QRow) == SIZE(RRow));
		vector<string> Fields;
		Split(QLabel, Fields, '/');
		asserta(SIZE(Fields) == 4);
		const string &Dom = Fields[0];
		const string &Fam = Fields[1];
		const string sTM = Fields[2];
		const string sPctId = Fields[3];
		double TM = StrToFloat(sTM);
		if (TM < MinTM || TM > MaxTM)
			continue;

		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char q = QRow[Col];
			char r = RRow[Col];
			if (isgap(q) || isgap(r))
				continue;
			uint iq = g_CharToLetterAmino[q];
			uint ir = g_CharToLetterAmino[r];
			if (iq >= 20 || ir >= 20)
				continue;
			LetterPairCount += 2;
			CountVec[iq] += 1;
			CountVec[ir] += 1;
			CountMx[iq][ir] += 1;
			CountMx[ir][iq] += 1;
			}
		}

	double SumFreq = 0;
	vector<double> Freqs;
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		double Freq = double(CountVec[i])/(LetterPairCount);
		SumFreq += Freq;
		Log("%c %.4f\n", c, Freq);
		Freqs.push_back(Freq);
		}
	asserta(feq(SumFreq, 1.0));

	FILE *f = CreateStdioFile(opt_output);

	double SumFreq2 = 0;
	double H = 0;
// RH=Relative entropy
//   =mean expected score per residue pair
//   =Kullback–Leibler divergence
// Altschul, SF (1991) "Amino acid substitution matrices from an
//   information theoretic perspective." JMB 219(3): 555-565.
	double RH = 0;
	fprintf(f, "aa\tfreq");
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		fprintf(f, "\t%c", c);
		}
	fprintf(f, "\n");

	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		double Freq = double(CountVec[i])/(2*LetterPairCount);
		SumFreq += Freq;
		fprintf(f, "%c", c);
		fprintf(f, "\t%.4f", Freq);
		for (uint j = 0; j < 20; ++j)
			{
			uint n = CountMx[i][j];
			double ObsFreq = double(n)/double(LetterPairCount);
			SumFreq2 += ObsFreq;
			double ExpFreq = double(Freqs[i]*Freqs[j]);
			double Ratio = ObsFreq/ExpFreq;
			double Score = log(Ratio);
			H -= ObsFreq*log(ObsFreq);
			RH += ObsFreq*Score;
			fprintf(f, "\t%.4f", Score);
			}
		fprintf(f, "\n");
		}
	asserta(feq(SumFreq2, 1.0));
	CloseStdioFile(f);
	ProgressLog("TM %.4f - %.4f, entropy %.3g, relative entropy %.3g\n",
	  MinTM, MaxTM, H, RH);
	}

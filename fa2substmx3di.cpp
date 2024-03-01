#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "outputfiles.h"

/***
palmscan2 \
  -fa2substmx3di ../scop40/scop40.fa2 \
  -input d:/a/res/dave_grant/scop40/scop40.3di \
  -output ../scop40/substmx3di_0.4_0.6.tsv \
  -log substmx3di.log -mintm 0.6 -maxtm 0.8 \
  -aln ../scop40/scop40_3di_0.4_0.6.fa2

Input
=====
(*) Pair-wise alignments of scop40 domains, created by
  palmscan2 \
    -tm_scop d:/int/scop40/out/domains_scop.cal \
	-output scop40.fa2 \
	-tsv scop40.tsv \
	-log tm.log

(*) 3di sequences for scop40 domains created by hacked foldseek 
  source in /d/int/foldseek_src (see also D:\a\doc\notebooks\2024-02-26_build_foldseek.txt)

# head /d/a/res/dave_grant/scop40/scop40.3di

>d1tyeb3/g.16.2.1
DCDQLQVDDDDDQVVSCVRDLQKKAAPDPPADVVDRRIDGPVVVVVSPGDPPRMDRD
>d1d1da1/a.28.3.1
DQCVVQADDPVDDPLRSLVVSLVSLVVDPDDPPCSFQVSLVRCCVDPDPLLVVLLVVDDPVDGGSPVSVVSSCVVPDDPD
>d2dyrb2/f.17.2.1
DDDPPDDDDDDDDDPLVVLVVVLVVVVVVVVVVVVVVVVVVVVVVVPDPDDDPDDDDDVPVVVVVVVVVVVVVVVSVPSV


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
# head /d/a/res/dave_grant/scop40/substmx3di_0.4_0.6.tsv
aa      freq    A       C       D       E       F       G       H      
A       0.0415  0.5330  -0.1004 -0.1395 -0.0768 -0.0880 -0.0375 -0.1575
C       0.0076  -0.1004 2.9846  -0.6905 -0.6766 -0.2257 -0.4758 0.0235
D       0.0269  -0.1395 -0.6905 0.9609  0.3571  -0.4771 0.0191  0.0140

***/

void GetScopDomFromLabel(const string &Label, string &Dom);

void cmd_fa2substmx3di()
	{
	SeqDB Fa2;
	Fa2.FromFasta(opt_fa2substmx3di, true);

	SeqDB Fa3di;
	Fa3di.FromFasta(opt_input);
	const uint SeqCount3di = Fa3di.GetSeqCount();
	map<string, string> DomToSeq3di;
	for (uint i = 0; i < SeqCount3di; ++i)
		{
		const string &Label = Fa3di.GetLabel(i);
		const string &Seq3di = Fa3di.GetSeq(i);
		string Dom;
		GetScopDomFromLabel(Label, Dom);
		DomToSeq3di[Dom] = Seq3di;
		}

	asserta(optset_mintm && optset_maxtm);
	double MinTM = opt_mintm;
	double MaxTM = opt_maxtm;

	vector<uint> CountVec(20);
	vector<vector<uint> > CountMx(20);
	for (uint i = 0; i < 20; ++i)
		CountMx[i].resize(20);

	const uint SeqCount = Fa2.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Processing");
		const string &QRow = Fa2.GetSeq(2*PairIndex);
		const string &RRow = Fa2.GetSeq(2*PairIndex+1);
		const string &QLabel = Fa2.GetLabel(2*PairIndex);
		const string &RLabel = Fa2.GetLabel(2*PairIndex+1);
		const uint ColCount = SIZE(QRow);
		asserta(SIZE(QRow) == SIZE(RRow));
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

		map<string, string>::const_iterator Qiter = DomToSeq3di.find(QDom);
		map<string, string>::const_iterator Riter = DomToSeq3di.find(RDom);
		if (Qiter == DomToSeq3di.end() || Riter == DomToSeq3di.end())
			continue;
		const string &QSeq3di = Qiter->second;
		const string &RSeq3di = Riter->second;
		uint QL = SIZE(QSeq3di);
		uint RL = SIZE(RSeq3di);
		uint QPos = 0;
		uint RPos = 0;
		string QRow3di;
		string RRow3di;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char q = QRow[Col];
			char r = RRow[Col];
			if (!isgap(q) && !isgap(r))
				{
				asserta(QPos < QL);
				asserta(RPos < RL);
				char q3di = QSeq3di[QPos];
				char r3di = RSeq3di[RPos];
				QRow3di += q3di;
				RRow3di += r3di;
				uint iq = g_CharToLetterAmino[q3di];
				uint ir = g_CharToLetterAmino[r3di];
				if (iq >= 20 || ir >= 20)
					continue;
				LetterPairCount += 2;
				CountVec[iq] += 1;
				CountVec[ir] += 1;
				CountMx[iq][ir] += 1;
				CountMx[ir][iq] += 1;
				}
			if (isgap(q))
				QRow3di += '-';
			else
				++QPos;

			if (isgap(r))
				RRow3di += '-';
			else
				++RPos;
			}
		if (g_faln)
			{
			fprintf(g_faln, ">%s\n", QLabel.c_str());
			fprintf(g_faln, "%s\n", QRow3di.c_str());
			fprintf(g_faln, ">%s\n", RLabel.c_str());
			fprintf(g_faln, "%s\n", RRow3di.c_str());
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

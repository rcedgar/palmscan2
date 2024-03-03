#include "myutils.h"
#include "pdbchain.h"
#include "xbinner.h"
#include "xprof.h"
#include "seqdb.h"
#include "alpha.h"
#include "outputfiles.h"

extern double g_XScoreMx[20][20];

static void WriteAln(FILE *f,
	const string &Label1, const string &Label2,
	const string &Row1, const string &Row2,
	const string &SSRow1, const string &SSRow2,
	const string &XRow1, const string &XRow2)
	{
	if (f == 0)
		return;
	fprintf(f, "\n\n");
	fprintf(f, "____________________________________________________________________\n");

	const uint ColCount = SIZE(Row1);
	asserta(SIZE(Row2) == ColCount);
	asserta(SIZE(SSRow1) == ColCount);
	asserta(SIZE(SSRow2) == ColCount);
	asserta(SIZE(XRow1) == ColCount);
	asserta(SIZE(XRow2) == ColCount);

	string Annot;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char x1 = XRow1[Col];
		char x2 = XRow2[Col];
		if (isgap(x1) || isgap(x2))
			Annot += ' ';
		else if (x1 == x2)
			Annot += '|';
		else
			{
			uint letter1 = g_CharToLetterAmino[x1];
			uint letter2 = g_CharToLetterAmino[x2];
			double Score = g_XScoreMx[letter1][letter2];
			if (Score > 0.5)
				Annot += '+';
			else if (Score > 0.2)
				Annot += '.';
			else
				Annot += ' ';
			}
		}

	const unsigned ROWLEN = 80;
	unsigned BlockCount = (ColCount + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned From = BlockIndex*ROWLEN;
		unsigned To = From + ROWLEN;
		if (To >= ColCount)
			To = ColCount;

		for (uint Pos = From; Pos < To; ++Pos) fputc(Row1[Pos], f);  fprintf(f, "  aa >%s", Label1.c_str()); fputc('\n', f);
		for (uint Pos = From; Pos < To; ++Pos) fputc(SSRow1[Pos], f);fprintf(f, "  ss >%s", Label1.c_str()); fputc('\n', f);
		for (uint Pos = From; Pos < To; ++Pos) fputc(XRow1[Pos], f); fprintf(f, "   X >%s", Label1.c_str()); fputc('\n', f);
		for (uint Pos = From; Pos < To; ++Pos) fputc(Annot[Pos], f); fputc('\n', f);
		for (uint Pos = From; Pos < To; ++Pos) fputc(XRow2[Pos], f); fprintf(f, "   X >%s", Label2.c_str()); fputc('\n', f);
		for (uint Pos = From; Pos < To; ++Pos) fputc(SSRow2[Pos], f);fprintf(f, "  ss >%s", Label2.c_str()); fputc('\n', f);
		for (uint Pos = From; Pos < To; ++Pos) fputc(Row2[Pos], f);  fprintf(f, "  aa >%s", Label2.c_str()); fputc('\n', f);
		fprintf(f, "\n");
		}
	}

void cmd_calfa2x()
	{
	const string &CalFN = opt_calfa2x;
	vector<PDBChain *> Chains;
	ReadChains(CalFN, Chains);
	const uint ChainCount = SIZE(Chains);

	SeqDB Alns;
	Alns.FromFasta(opt_input, true);
	const uint AlnsSeqCount = Alns.GetSeqCount();
	asserta(AlnsSeqCount%2 == 0);
	const uint AlnsPairCount = AlnsSeqCount/2;

	XProf::InitScoreTable();
	XBinner::InitCentroids();
	XBinner XB;

	FILE *f = CreateStdioFile(opt_output);
	XProf XP;
	map<string, string> LabelToXSeq;
	map<string, string> LabelToSS;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Converting cal to X");
		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;

		XP.Init(Chain);
		const uint L = Chain.GetSeqLength();
		string XSeq;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			vector<double> Values;
			XP.GetFeatures(Pos, Values);
			uint XLetter = XB.GetLetter(Values);
			char XChar = g_LetterToCharAmino[XLetter];
			XSeq += XChar;
			}
		string SS;
		Chain.GetSS(SS);
		LabelToXSeq[Label] = XSeq;
		LabelToSS[Label] = SS;
		SeqToFasta(f, Label, XSeq);
		}
	CloseStdioFile(f);

	FILE *f2 = CreateStdioFile(opt_output2);
	FILE *f4 = CreateStdioFile(opt_output4);
	uint Lerr = 0;
	uint nok = 0;
	uint notfound = 0;
	for (uint PairIndex = 0; PairIndex < AlnsPairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, AlnsPairCount, "Pairs, %u Lerrs, %u notfound",
		  Lerr, notfound);
		const string &Label1 = Alns.GetLabel(2*PairIndex);
		const string &Label2 = Alns.GetLabel(2*PairIndex+1);

		vector<string> Fields;
		Split(Label1, Fields, '/');
		asserta(SIZE(Fields) == 4);
		string FixedLabel1 = Fields[0] + "/" + Fields[1];

		const map<string, string>::const_iterator iter1 = LabelToXSeq.find(FixedLabel1);
		const map<string, string>::const_iterator iter2 = LabelToXSeq.find(Label2);
		if (iter1 == LabelToXSeq.end() || iter2 == LabelToXSeq.end())
			{
			++notfound;
			continue;
			}

		const string &XSeq1 = iter1->second;
		const string &XSeq2 = iter2->second;

		const map<string, string>::const_iterator iterss1 = LabelToSS.find(FixedLabel1);
		const map<string, string>::const_iterator iterss2 = LabelToSS.find(Label2);
		asserta(iterss1 != LabelToSS.end() && iterss2 != LabelToSS.end());

		const string &SS1 = iterss1->second;
		const string &SS2 = iterss2->second;

		const uint L1 = SIZE(XSeq1);
		const uint L2 = SIZE(XSeq2);

		const string &Row1 = Alns.GetSeq(2*PairIndex);
		const string &Row2 = Alns.GetSeq(2*PairIndex+1);

		string XRow1;
		string XRow2;

		string SSRow1;
		string SSRow2;

		const uint ColCount = SIZE(Row1);
		asserta(SIZE(Row2) == ColCount);
		uint Pos1 = 0;
		uint Pos2 = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char c1 = Row1[Col];
			char c2 = Row2[Col];
			bool gap1 = isgap(c1);
			bool gap2 = isgap(c2);

			if (gap1)
				{
				XRow1 += '-';
				SSRow1 += ' ';
				}
			else
				{
				XRow1 += XSeq1[Pos1];
				SSRow1 += SS1[Pos1];
				++Pos1;
				}

			if (gap2)
				{
				XRow2 += '-';
				SSRow2 += ' ';
				}
			else
				{
				XRow2 += XSeq2[Pos2];
				SSRow2 += SS2[Pos2];
				++Pos2;
				}
			}
		if (Pos1 != L1 || Pos2 != L2)
			{
			++Lerr;
			continue;
			}
		++nok;
		SeqToFasta(f2, Label1, XRow1);
		SeqToFasta(f2, Label2, XRow2);

		SeqToFasta(f4, Label1 + "/aa", Row1);
		SeqToFasta(f4, Label1 + "/ss", SSRow1);
		SeqToFasta(f4, Label1 + "/x", XRow1);
		SeqToFasta(f4, Label2 + "/aa", Row2);
		SeqToFasta(f4, Label2 + "/ss", SSRow2);
		SeqToFasta(f4, Label2 + "/x", XRow2);

		WriteAln(g_faln, Label1, Label2, Row1, Row2, SSRow1, SSRow2, XRow1, XRow2);
		}
	CloseStdioFile(f2);
	CloseStdioFile(f4);
	ProgressLog("%u ok, %u Lerr\n", nok, Lerr);
	}

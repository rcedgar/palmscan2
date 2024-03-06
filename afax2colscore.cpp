#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "outputfiles.h"

extern double g_Blosum62[20][20];
extern double g_XScoreMx[20][20];
static const double GAPSCORE = 0;
static const double MAXGAPFRACT = 0.2;

void MakeRow(const string &InRow, const string &InSeq,
  char GapChar, string &OutRow);

/***
Input is structure MSA & aa sequences (alignment if any ignored)
	-afax2colscore seqs.afax	# X (structure alphabet)
	-input seqs.afa/fa			# aa alphabet
***/
static void GetFreqs(string &Col, vector<double> &Freqs)
	{
	Freqs.clear();
	Freqs.resize(20, 0);
	vector<uint> Counts;
	const uint SeqCount = SIZE(Col);
	uint n = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		char c = Col[i];
		uint Letter = g_CharToLetterAmino[c];
		if (Letter >= 20)
			continue;
		++Counts[Letter];
		++n;
		}
	if (n == 0)
		return;
	for (uint i = 0; i < 20; ++i)
		Freqs[i] = double(Counts[i])/n;
	}

static double GetScore(const double Mx[20][20], const string &Col)
	{
	const uint SeqCount = SIZE(Col);
	uint n = 0;
	double SumScore = 0;
	uint GapCount = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		char ci = Col[i];
		if (isgap(ci))
			{
			++GapCount;
			continue;
			}
		uint Letteri = g_CharToLetterAmino[ci];
		if (Letteri >= 20)
			continue;
		for (uint j = 1; j < i; ++j)
			{
			uint Letterj = g_CharToLetterAmino[Col[j]];
			if (Letterj >= 20)
				continue;
			SumScore += Mx[Letteri][Letterj];
			++n;
			}
		}
	double Score = GAPSCORE*GapCount;
	if (n > 0)
		Score += SumScore/n;
	return Score;
	}

//static double GetScore(const double Mx[20][20], const vector<double> &Freqs)
//	{
//	asserta(SIZE(Freqs) == 20);
//	double SumScore = 0;
//	for (uint i = 0; i < 20; ++i)
//		{
//		double fi = Freqs[i];
//		SumScore += Mx[i][i]*fi*fi;
//		for (uint j = 0; j < i; ++j)
//			{
//			double fj = Freqs[j];
//			SumScore += 2*Mx[i][j]*fi*fj;
//			}
//		}
//	return SumScore;
//	}
//
//static double GetXScore(const vector<double> &Freqs)
//	{
//	return GetScore(g_XScoreMx, Freqs);
//	}
//
//static double GetAAScore(const vector<double> &Freqs)
//	{
//	return GetScore(g_Blosum62, Freqs);
//	}

void cmd_afax2colscore()
	{
	asserta(optset_input);
	const string &MSAXFN = opt_afax2colscore;

	SeqDB Seqs;
	Seqs.FromFasta(opt_input, false);
	Seqs.SetLabelToIndex();

	SeqDB MSAx;
	MSAx.FromFasta(MSAXFN, true);

	SeqDB Input;
	Input.FromFasta(opt_input, false);

	SeqDB MSAX;
	MSAX.FromFasta(MSAXFN, true);
	const uint SeqCount = MSAX.GetSeqCount();
	Input.SetLabelToIndex();

	vector<string> XRows;
	vector<string> AARows;

	const uint ColCount = MSAX.GetColCount();
	for (uint SeqIndexX = 0; SeqIndexX < SeqCount; ++SeqIndexX)
		{
		const string &Label = MSAX.GetLabel(SeqIndexX);
		string AASeq;
		bool Found = Input.GetSeqByLabel(Label, AASeq, false);
		if (!Found)
			{
			Warning("Not found in input >%s", Label.c_str());
			continue;
			}

		const string &XRow = MSAX.GetSeq(SeqIndexX);

		string AARow;
		MakeRow(XRow, AASeq, '-', AARow);

		asserta(SIZE(XRow) == ColCount);
		asserta(SIZE(AARow) == ColCount);

		XRows.push_back(XRow);
		AARows.push_back(AARow);
		}
	ProgressLog("Done.\n");

	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		string XCol;
		string AACol;
		uint GapCount = 0;
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			char x = XRows[SeqIndex][ColIndex];
			char aa = AARows[SeqIndex][ColIndex];
			XCol.push_back(x);
			AACol.push_back(aa);
			if (isgap(aa))
				++GapCount;
			}
		if (double(GapCount)/SeqCount > MAXGAPFRACT)
			continue;
		//vector<double> XFreqs;
		//vector<double> AAFreqs;
		//GetFreqs(XCol, XFreqs);
		//GetFreqs(AACol, AAFreqs);
		double AAScore = GetScore(g_Blosum62, AACol);
		double XScore = GetScore(g_XScoreMx, XCol);
		double Score = max(AAScore, XScore);
		if (g_ftsv != 0)
			fprintf(g_ftsv, "%u\t%.1f\t%.1f\t%.1f\t'%s'\t'%s'\n",
			  ColIndex, Score, XScore, AAScore, XCol.c_str(), AACol.c_str());
		}
	}

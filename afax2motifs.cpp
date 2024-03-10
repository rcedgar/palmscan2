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
	-afax2smm seqs.afax	# X (structure alphabet)
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

static void GetGaplessBlocks(const vector<double> &GapFracts, 
  double MaxFract, uint MinLength,
  vector<uint> &Los, vector<uint> &Lengths)
	{
	Los.clear();
	Lengths.clear();
	const uint ColCount = SIZE(GapFracts);
	uint Lo = UINT_MAX;
	for (uint Col = 0; Col <= ColCount; ++Col)
		{
		double GapFract = (Col == ColCount ? 0 : GapFracts[Col]);
		if (GapFract <= MaxFract)
			{
			if (Lo == UINT_MAX)
				Lo = Col;
			}
		else
			{
			if (Lo != UINT_MAX)
				{
				uint Length = Col - Lo;
				if (Length > MinLength)
					{
					Los.push_back(Lo);
					Lengths.push_back(Length);
					}
				Lo = UINT_MAX;
				}
			}
		}
	}

static void GetGapFracts(const vector<string> &AARows,
  vector<double> &GapFracts)
	{
	GapFracts.clear();
	const uint SeqCount = SIZE(AARows);
	if (SeqCount == 0)
		return;
	const uint ColCount = SIZE(AARows[0]);
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		uint GapCount = 0;
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			asserta(SeqIndex < SIZE(AARows));
			asserta(ColIndex < SIZE(AARows[SeqIndex]));
			char aa = AARows[SeqIndex][ColIndex];
			if (isgap(aa))
				++GapCount;
			}
		GapFracts.push_back(double(GapCount)/SeqCount);
		}
	}

static void MakeRows(const SeqDB &MSAX, const SeqDB &Input,
  vector<string> &Labels,
  vector<string> &XRows,
  vector<string> &AARows)
	{
	Labels.clear();
	XRows.clear();
	AARows.clear();
	asserta(!Input.m_LabelToIndex.empty());
	const uint SeqCount = MSAX.GetSeqCount();
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

		Labels.push_back(Label);
		XRows.push_back(XRow);
		AARows.push_back(AARow);
		}
	}

static uint ColToPos(uint Col, const string &Row)
	{
	uint Pos = 0;
	const uint ColCount = SIZE(Row);
	asserta(Col < ColCount);
	for (uint Col2 = 0; Col2 < ColCount; ++Col2)
		{
		char c = Row[Col2];
		if (!isgap(c))
			{
			if (Col2 >= Col)
				return Pos;
			++Pos;
			}
		}
	asserta(false);
	return UINT_MAX;
	}

static uint GetMotifs(const vector<string> &Labels,
	const vector<string> &AARows,
	const vector<uint> &ColLos, 
	const vector<uint> &ColLengths,
	vector<vector<uint> > &MotifStartsVec,
	vector<vector<string> > &MotifSeqsVec)
	{
	MotifStartsVec.clear();
	MotifSeqsVec.clear();
	const uint SeqCount = SIZE(Labels);
	asserta(SIZE(AARows) == SeqCount);
	asserta(SeqCount > 0);
	const uint MotifCount = SIZE(ColLos);
	asserta(SIZE(ColLengths) == MotifCount);

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		vector<uint> MotifStarts;
		vector<string> MotifSeqs;

		const string &AARow = AARows[SeqIndex];
		const string &Label = Labels[SeqIndex];
		for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
			{
			uint ColLo = ColLos[MotifIndex];
			uint ColLength = ColLengths[MotifIndex];
			uint PosLo = ColToPos(ColLo, AARow);
			string MotifSeq;
			for (uint Col = ColLo; Col < ColLo + ColLength; ++Col)
				{
				char c = AARow[Col];
				MotifSeq += c;
				}
			MotifStarts.push_back(PosLo);
			MotifSeqs.push_back(MotifSeq);
			}
		MotifStartsVec.push_back(MotifStarts);
		MotifSeqsVec.push_back(MotifSeqs);
		}
	return MotifCount;
	}

static void WriteMotifs(FILE *f,
	const vector<string> &Labels,
	vector<vector<uint> > &MotifStartsVec,
	vector<vector<string> > &MotifSeqsVec)
	{
	if (f == 0)
		return;
	const uint SeqCount = SIZE(Labels);
	asserta(SeqCount > 0);
	asserta(SIZE(MotifStartsVec) == SeqCount);
	asserta(SIZE(MotifSeqsVec) == SeqCount);
	const uint MotifCount = SIZE(MotifStartsVec[0]);
	fprintf(f, "Label");
	for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
		fprintf(f, "\tMf%u", MotifIndex+1);
	fprintf(f, "\n");
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = Labels[SeqIndex];
		const vector<uint> &MotifStarts = MotifStartsVec[SeqIndex];
		const vector<string> &MotifSeqs = MotifSeqsVec[SeqIndex];
		asserta(SIZE(MotifStarts) == MotifCount);
		asserta(SIZE(MotifSeqs) == MotifCount);
		fprintf(f, "%s", Label.c_str());
		for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
			fprintf(f, "\t%s", MotifSeqs[MotifIndex].c_str());
		fprintf(f, "\n");
		}
	}

static void WriteMotifsMSAs(const string &FNPrefix,
	const vector<string> &Labels,
	vector<vector<uint> > &MotifStartsVec,
	vector<vector<string> > &MotifSeqsVec)
	{
	if (FNPrefix == "")
		return;

	const uint SeqCount = SIZE(Labels);
	asserta(SeqCount > 0);
	asserta(SIZE(MotifStartsVec) == SeqCount);
	asserta(SIZE(MotifSeqsVec) == SeqCount);
	const uint MotifCount = SIZE(MotifStartsVec[0]);
	for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
		{
		string FastaFN = FNPrefix;
		Psa(FastaFN, "mf%u", MotifIndex + 1);
		FILE *f = CreateStdioFile(FastaFN);
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			const string &Label = Labels[SeqIndex];
			uint Start = MotifStartsVec[SeqIndex][MotifIndex];
			const vector<string> &MotifSeqs = MotifSeqsVec[SeqIndex];
			fprintf(f, ">%s/%u\n", Label.c_str(), Start+1);
			fprintf(f, "%s\n", MotifSeqsVec[SeqIndex][MotifIndex].c_str());
			}
		CloseStdioFile(f);
		}
	}

const char *GetPymolColor(uint MotifIndex)
	{
	static const char *Colors[] =
		{
		"tv_red",
		"tv_green",
		"tv_blue",
		"cyan",
		"tv_orange",
		"magenta",
		"lightpink",
		"white",
		};
	static const size_t ColorCount = sizeof(Colors)/sizeof(Colors[0]);
	return Colors[MotifIndex%ColorCount];
	}

static void WritePMLs(FILE *f,
	const vector<string> &Labels,
	vector<vector<uint> > &MotifStartsVec,
	vector<vector<string> > &MotifSeqsVec)
	{
	if (f == 0)
		return;
	const uint SeqCount = SIZE(Labels);
	asserta(SeqCount > 0);
	asserta(SIZE(MotifStartsVec) == SeqCount);
	asserta(SIZE(MotifSeqsVec) == SeqCount);
	const uint MotifCount = SIZE(MotifStartsVec[0]);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = Labels[SeqIndex];
		fprintf(f, "\n\n%s\n", Label.c_str());
		fprintf(f, "color gray70\n");
		fprintf(f, "bg white\n");
		const vector<uint> &MotifStarts = MotifStartsVec[SeqIndex];
		const vector<string> &MotifSeqs = MotifSeqsVec[SeqIndex];
		asserta(SIZE(MotifStarts) == MotifCount);
		asserta(SIZE(MotifSeqs) == MotifCount);
		for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
			{
			char M = 'A' + MotifIndex;
			const char *Color = GetPymolColor(MotifIndex);
			const char *MotifSeq = MotifSeqs[MotifIndex].c_str();
			fprintf(f, "select Mtf%c, pepseq %s\n", M, MotifSeq);
			fprintf(f, "color %s, Mtf%c\n", Color, M);
			}
		fprintf(f, "deselect\n");
		fprintf(f, "\n");
		}
	}

void cmd_afax2motifs()
	{
	asserta(optset_input);
	const string &MSAXFN = opt_afax2motifs;

	SeqDB Seqs;
	Seqs.FromFasta(opt_input, false);
	Seqs.SetLabelToIndex();

	SeqDB Input;
	Input.FromFasta(opt_input, false);

	SeqDB MSAX;
	MSAX.FromFasta(MSAXFN, true);
	const uint SeqCount = MSAX.GetSeqCount();
	Input.SetLabelToIndex();

	vector<string> Labels;
	vector<string> XRows;
	vector<string> AARows;
	MakeRows(MSAX, Input, Labels, XRows, AARows);

	vector<double> GapFracts;
	GetGapFracts(AARows, GapFracts);

	const double MaxFract = 0.2;
	const uint MinLength = 7;
	vector<uint> ColLos;
	vector<uint> ColLengths;
	GetGaplessBlocks(GapFracts, MaxFract, MinLength, ColLos, ColLengths);
	const uint MotifCount = SIZE(ColLos);
	asserta(SIZE(ColLengths) == MotifCount);

	vector<vector<uint> > MotifStartsVec;
	vector<vector<string> > MotifSeqsVec;
	const uint MotifCount2 = GetMotifs(Labels, AARows, ColLos, ColLengths,
	  MotifStartsVec, MotifSeqsVec);
	asserta(MotifCount2 == MotifCount);

	WriteMotifs(g_ftsv, Labels, MotifStartsVec, MotifSeqsVec);
	WritePMLs(g_fpml, Labels, MotifStartsVec, MotifSeqsVec);
	if (optset_output2)
		WriteMotifsMSAs(opt_output2, Labels, MotifStartsVec, MotifSeqsVec);
	}

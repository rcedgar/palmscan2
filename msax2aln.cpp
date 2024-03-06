#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "outputfiles.h"
#include <map>
#include <set>

static void MakeRow(const string &InRow,
  const string &InSeq, char GapChar, string &OutRow)
	{
	const uint ColCount = SIZE(InRow);
	uint Pos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = InRow[Col];
		if (isgap(c))
			{
			OutRow += GapChar;
			continue;
			}
		asserta(Pos < SIZE(InSeq));
		OutRow += InSeq[Pos++];
		}
	asserta(Pos == SIZE(InSeq));
	}

void cmd_msax2aln()
	{
	const string &MSAXFN = opt_msax2aln; // struct MSA in 3DI or X alphabet
	const string &CalFN = opt_calin;	 // .cal with chains (for aa & ss)

	SeqDB MSAX;
	MSAX.FromFasta(MSAXFN, true);
	set<string> LabelSet;
	const uint SeqCount = MSAX.GetSeqCount();
	for (uint i = 0; i < SeqCount; ++i)
		LabelSet.insert(MSAX.GetLabel(i));
	if (SIZE(LabelSet) != SeqCount)
		Die("Dupe labels");

	vector<PDBChain *> Chains;
	ReadChains(CalFN, Chains);
	const uint ChainCount = SIZE(Chains);
	map<string, string> LabelToSS;
	map<string, string> LabelToSeq;
	set<string> FoundSet;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Converting cal to X");
		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;
		if (LabelSet.find(Label) == LabelSet.end())
			continue;
		FoundSet.insert(Label);

		const string &Seq = Chain.m_Seq;

		string SS;
		Chain.GetSS(SS);
		asserta(SIZE(Seq) == SIZE(SS));
		LabelToSeq[Label] = Seq;
		LabelToSS[Label] = SS;
		}
	if (FoundSet.empty())
		Die("None found in .cal");

	vector<string> XRows;
	vector<string> AARows;
	vector<string> SSRows;
	const uint ColCount = MSAX.GetColCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = MSAX.GetLabel(SeqIndex);
		if (FoundSet.find(Label) == FoundSet.end())
			{
			Warning("Not found in .cal >%s", Label.c_str());
			continue;
			}
		const string &XRow = MSAX.GetSeq(SeqIndex);
		const string &AASeq = LabelToSeq[Label];
		const string &SS = LabelToSS[Label];
		string AARow;
		string SSRow;
		MakeRow(XRow, AASeq, '-', AARow);
		MakeRow(XRow, SS, ' ', SSRow);

		asserta(SIZE(XRow) == ColCount);
		asserta(SIZE(AARow) == ColCount);
		asserta(SIZE(SSRow) == ColCount);

		XRows.push_back(XRow);
		AARows.push_back(AARow);
		SSRows.push_back(SSRow);
		}
	const uint AlignedSeqCount = SIZE(XRows);
	asserta(SIZE(AARows) == AlignedSeqCount);
	asserta(SIZE(SSRows) == AlignedSeqCount);

	const unsigned ROWLEN = 80;

	if (optset_aln)
		{
		FILE *f = g_faln;
		unsigned BlockCount = (ColCount + ROWLEN - 1)/ROWLEN;
		for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
			{
			unsigned From = BlockIndex*ROWLEN;
			unsigned To = From + ROWLEN;
			if (To >= ColCount)
				To = ColCount;

			for (uint SeqIndex = 0; SeqIndex < AlignedSeqCount; ++SeqIndex)
				{
				const string &Label = MSAX.GetLabel(SeqIndex);
				const string &XRow = XRows[SeqIndex];
				const string &AARow = AARows[SeqIndex];
				const string &SSRow = SSRows[SeqIndex];

				for (uint Pos = From; Pos < To; ++Pos) fputc(SSRow[Pos], f); fprintf(f, "  ss >%s", Label.c_str()); fputc('\n', f);
				for (uint Pos = From; Pos < To; ++Pos) fputc(AARow[Pos], f); fprintf(f, "  aa >%s", Label.c_str()); fputc('\n', f);
				for (uint Pos = From; Pos < To; ++Pos) fputc(XRow[Pos], f);  fprintf(f, "   X >%s", Label.c_str()); fputc('\n', f);
				fprintf(f, "\n");
				}
			for (uint Pos = From; Pos < To; ++Pos) fputc('#', f); fputc('\n', f); fputc('\n', f);
			}
		}

	if (optset_output2)
		{
		FILE *f2 = CreateStdioFile(opt_output2);
		for (uint SeqIndex = 0; SeqIndex < AlignedSeqCount; ++SeqIndex)
			{
			const string &Label = MSAX.GetLabel(SeqIndex);
			const string &AARow = AARows[SeqIndex];
			SeqToFasta(f2, Label, AARow);
			}
		}
	}

#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "outputfiles.h"
#include <map>
#include <set>

void MakeRow(const string &InRow, const string &InSeq,
  char GapChar, string &OutRow)
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

void cmd_afax2afa()
	{
	const string &MSAXFN = opt_afax2afa; // struct MSA in 3Di or X alphabet
	asserta(optset_output);
	FILE *fOut = CreateStdioFile(opt_output);

	SeqDB Input;
	Input.FromFasta(opt_input, false);

	SeqDB MSAX;
	MSAX.FromFasta(MSAXFN, true);
	const uint SeqCount = MSAX.GetSeqCount();
	Input.SetLabelToIndex();

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
		SeqToFasta(fOut, Label, AARow);
		}
	CloseStdioFile(fOut);
	}

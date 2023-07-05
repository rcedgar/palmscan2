#include "myutils.h"
#include "shapesearcher.h"
#include "rdrpsearcher.h"
#include "motifsettings.h"

void GetTrainingMotifs(const string &FileName,
  const vector<PDBChain *> &Chains, vector<string> &ChainLabels,
  vector<string> &MotifNames, vector<uint> &MotifLengths,
  vector<vector<string> > &MotifSeqsVec)
	{
	ChainLabels.clear();
	MotifNames.clear();
	MotifLengths.clear();
	MotifSeqsVec.clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	vector<string> Fields;
	Split(Line, Fields, '\t');
	asserta(Fields[0] == "Label");
	asserta(SIZE(Fields) > 1);
	const uint MotifCount = SIZE(Fields) - 1;
	uint IxA = UINT_MAX;
	uint IxB = UINT_MAX;
	uint IxC = UINT_MAX;
	for (uint i = 0; i < MotifCount; ++i)
		{
		const string &Name = Fields[i+1];
		if (Name == "A")
			IxA = i;
		else if (Name == "B")
			IxB = i;
		else if (Name == "C")
			IxC = i;
		MotifLengths.push_back(UINT_MAX);
		MotifNames.push_back(Name);
		}

	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == MotifCount + 1);

		const string &ChainLabel = Fields[0];
		ChainLabels.push_back(ChainLabel);

		vector<string> MotifSeqs;
		for (uint i = 0; i < MotifCount; ++i)
			MotifSeqs.push_back(Fields[i+1]);
		MotifSeqsVec.push_back(MotifSeqs);

		if (IxA != UINT_MAX)
			CheckA(MotifSeqs[IxA]);
		if (IxB != UINT_MAX)
			CheckB(MotifSeqs[IxB]);
		if (IxC != UINT_MAX)
			CheckC(MotifSeqs[IxC]);

		for (uint i = 0; i < MotifCount; ++i)
			{
			const string &MotifSeq = MotifSeqs[i];
			if (MotifSeq == "" || MotifSeq == ".")
				continue;
			uint L = SIZE(MotifSeq);
			if (MotifLengths[i] == UINT_MAX)
				MotifLengths[i] = L;
			else
				asserta(L == MotifLengths[i]);
			}
		}
	CloseStdioFile(f);
	WriteMotifSettings(g_fLog);
	ClearMotifGoodCounts();
	}

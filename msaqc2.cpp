#include "myutils.h"
#include "alpha.h"
#include "seqdb.h"
#include "rdrpmodel.h"
#include "msaqc2.h"
#include "outputfiles.h"

static const double MIN_MODEL_SCORE = 10.0;

void MSAQC2::GetColMaps(const string &UnalignedSeq, const string &AlignedSeq,
  vector<uint> &UnalignedPosToAlignedCol,
  vector<uint> &AlignedColToUnalginedPos)
	{
	UnalignedPosToAlignedCol.clear();
	AlignedColToUnalginedPos.clear();

	const uint L = SIZE(UnalignedSeq);
	const uint m_ColCount = SIZE(AlignedSeq);
	asserta(L <= m_ColCount);

	UnalignedPosToAlignedCol.resize(L, UINT_MAX);
	AlignedColToUnalginedPos.resize(m_ColCount, UINT_MAX);

	uint Pos = 0;
	for (uint Col = 0; Col < m_ColCount; ++Col)
		{
		char c = AlignedSeq[Col];
		if (isgap(c))
			continue;
		asserta(toupper(UnalignedSeq[Pos]) == toupper(AlignedSeq[Col]));
		asserta(UnalignedPosToAlignedCol[Pos] == UINT_MAX);
		asserta(AlignedColToUnalginedPos[Col] == UINT_MAX);

		UnalignedPosToAlignedCol[Pos] = Col;
		AlignedColToUnalginedPos[Col] = Pos;
		++Pos;
		}
	asserta(Pos == L);
	}

uint MSAQC2::GetMotifIndex(char c) const
	{
	if (c == 'A')
		return 0;
	else if (c == 'B')
		return 1;
	else if (c == 'C')
		return 2;
	asserta(false);
	return UINT_MAX;
	}

uint MSAQC2::GetTopIx(const vector<uint> &v, uint MinCount) const
	{
	asserta(SIZE(v) > 0);
	uint Max = 0;
	uint Ix = UINT_MAX;
	for (uint i = 0; i < SIZE(v); ++i)
		{
		uint n = v[i];
		if (n > Max)
			{
			Max = n;
			Ix = i;
			}
		}
	return Ix;
	}

float MSAQC2::GetGapFract(uint ColIndex) const
	{
	uint GapCount = 0;
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		asserta(SeqIndex < SIZE(m_AlignedSeqs));
		const string &Seq = m_AlignedSeqs[SeqIndex];
		asserta(ColIndex < SIZE(Seq));
		char c = Seq[ColIndex];
		if (isgap(c))
			++GapCount;
		}
	float Fract = float(GapCount)/float(m_ColCount);
	return Fract;
	}

void MSAQC2::Init(const SeqDB &MSA)
	{
	m_Mod.Clear();

	asserta(MSA.IsAligned());
	m_Ls.clear();
	m_Ls.push_back(LA);
	m_Ls.push_back(LB);
	m_Ls.push_back(LC);

	m_MSA = &MSA;
	
	m_ColCount = MSA.GetColCount();
	m_SeqCount = MSA.GetSeqCount();
	
	m_UnalignedSeqs.clear();
	m_AlignedSeqs.clear();

	m_UnalignedSeqs.resize(m_SeqCount);
	m_AlignedSeqs.resize(m_SeqCount);
	m_ColIndexToGapFract.resize(m_ColCount, FLT_MAX);
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, m_SeqCount, "Strip gaps");
		m_AlignedSeqs[SeqIndex] = MSA.GetSeq(SeqIndex);
		MSA.GetSeq_StripGaps(SeqIndex, m_UnalignedSeqs[SeqIndex]);
		}

	for (uint ColIndex = 0; ColIndex < m_ColCount; ++ColIndex)
		{
		ProgressStep(ColIndex, m_ColCount, "Gap fracts");
		float GapFract = GetGapFract(ColIndex);
		m_ColIndexToGapFract[ColIndex] = GapFract;
		}

	m_UnalignedPosToAlignedColVec.clear();
	m_AlignedColToUnalignedPosVec.clear();

	m_UnalignedPosToAlignedColVec.resize(m_SeqCount);
	m_AlignedColToUnalignedPosVec.resize(m_SeqCount);
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, m_SeqCount, "Col maps");
		const string &UnalignedSeq = m_UnalignedSeqs[SeqIndex];
		const string &AlignedSeq = m_AlignedSeqs[SeqIndex];
		GetColMaps(UnalignedSeq, AlignedSeq,
		  m_UnalignedPosToAlignedColVec[SeqIndex],
		  m_AlignedColToUnalignedPosVec[SeqIndex]);
		}

	GetRdrpModel(m_Mod);
	m_RS.Init(m_Mod);

	m_SeqIndexToPosVec.clear();
	m_SeqIndexToColVec.clear();
	m_SeqIndexToPosVec.resize(m_SeqCount);
	m_SeqIndexToColVec.resize(m_SeqCount);
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		m_SeqIndexToPosVec[SeqIndex].resize(NPOS, UINT_MAX);
		m_SeqIndexToColVec[SeqIndex].resize(NPOS, UINT_MAX);
		}

	vector<vector<uint> > PosToColToCountVec;
	PosToColToCountVec.resize(NPOS);
	for (uint i = 0; i < NPOS; ++i)
		PosToColToCountVec[i].resize(m_ColCount, 0);

	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, m_SeqCount, "PSSM search");

		const string Label = MSA.GetLabel(SeqIndex);

		asserta(SeqIndex < SIZE(m_UnalignedSeqs));
		const string &UnalignedSeq = m_UnalignedSeqs[SeqIndex];

		asserta(SeqIndex < SIZE(m_UnalignedPosToAlignedColVec));
		const vector<uint> &UnalignedPosToAlignedCol =
		  m_UnalignedPosToAlignedColVec[SeqIndex];

		m_RS.Search(Label, UnalignedSeq);
		bool IsHit = (m_RS.m_TopPalmHit.m_Score >= MIN_MODEL_SCORE);
		if (!IsHit)
			continue;

		uint APos = m_RS.GetMotifPos(0);
		uint BPos = m_RS.GetMotifPos(1);
		uint CPos = m_RS.GetMotifPos(2);

		vector<uint> PosVec;

		//m_Mod.GetSuperMotifPositions(PosVec);
		//if (PosVec.empty())
		//	continue;

		PosVec.push_back(APos+3);
		PosVec.push_back(BPos+1);
		PosVec.push_back(CPos+3);
		asserta(SIZE(PosVec) == NPOS);
		for (uint i = 0; i < NPOS; ++i)
			{
			uint Pos = PosVec[i];
			m_SeqIndexToPosVec[SeqIndex][i] = Pos;
			if (Pos == UINT_MAX)
				continue;
			asserta(Pos < SIZE(UnalignedSeq));
			uint Col = UnalignedPosToAlignedCol[Pos];
			m_SeqIndexToColVec[SeqIndex][i] = Col;
			asserta(Col < m_ColCount);
			PosToColToCountVec[i][Col] += 1;
			}
		}

	const uint Majority = m_SeqCount/2 + 1;
	
	m_ConsensusColVec.clear();
	m_ConsensusColVec.resize(NPOS);
	for (uint i = 0; i < NPOS; ++i)
		{
		uint Col = GetTopIx(PosToColToCountVec[i], Majority);
		m_ConsensusColVec[i] = Col;
		}
	}

bool MSAQC2::CheckIncreasingOrder(const vector<uint> &v) const
	{
	const uint N = SIZE(v);
	for (uint i = 1; i < N; ++i)
		if (v[i] <= v[i-1])
			return false;
	return true;
	}

void MSAQC2::GetCorrectlyAlignedLetterCount(uint SeqIndex,
  uint &CorrectCount, uint &LetterCount) const
	{
	CorrectCount = 0;
	LetterCount = 0;
	asserta(SeqIndex < SIZE(m_UnalignedPosToAlignedColVec));
	const vector<uint> &PosToCol = m_UnalignedPosToAlignedColVec[SeqIndex];

	asserta(SeqIndex < SIZE(m_SeqIndexToPosVec));
	const vector<uint> &PosVec = m_SeqIndexToPosVec[SeqIndex];

	asserta(SIZE(m_ConsensusColVec) == NPOS);
	for (uint i = 0; i < NPOS; ++i)
		{
		uint MotifCol = m_ConsensusColVec[i];
		uint Pos = PosVec[i];
		if (Pos == UINT_MAX)
			continue;
		++LetterCount;
		asserta(Pos < SIZE(PosToCol));
		uint Col = PosToCol[Pos];
		if (Col == MotifCol)
			++CorrectCount;
		}
	}

static uint GetLoCol(const vector<uint> &ColVec)
	{
	uint Lo = UINT_MAX;
	for (uint i = 0; i < SIZE(ColVec); ++i)
		{
		uint Col = ColVec[i];
		if (Col == UINT_MAX)
			continue;
		if (Col < Lo)
			Lo = Col;
		}
	return Lo;
	}

static uint GetHiCol(const vector<uint> &ColVec)
	{
	uint Hi = UINT_MAX;
	for (uint i = 0; i < SIZE(ColVec); ++i)
		{
		uint Col = ColVec[i];
		if (Col == UINT_MAX)
			continue;
		if (Hi == UINT_MAX || Col > Hi)
			Hi = Col;
		}
	return Hi;
	}

uint MSAQC2::GetCorrectSeqCount() const
	{
	uint CorrectSeqCount = 0;
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		uint CorrectCount;
		uint LetterCount;
		GetCorrectlyAlignedLetterCount(SeqIndex, CorrectCount, LetterCount);
		asserta(CorrectCount <= LetterCount);
		if (CorrectCount == LetterCount)
			++CorrectSeqCount;
		}
	return CorrectSeqCount;
	}

void MSAQC2::GetCorrectLetterCount(uint &CorrectCount, uint &LetterCount) const
	{
	CorrectCount = 0;
	LetterCount = 0;
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		uint CorrectCount1;
		uint LetterCount1;
		GetCorrectlyAlignedLetterCount(SeqIndex, CorrectCount1, LetterCount1);
		asserta(CorrectCount1 <= LetterCount1);
		CorrectCount += CorrectCount1;
		LetterCount += LetterCount1;
		}
	}

void MSAQC2::WriteAlnA(FILE *f) const
	{
	if (f == 0)
		return;
	if (m_ConsensusColVec.empty())
		return;
	asserta(SIZE(m_ConsensusColVec) == NPOS);
	uint Col1 = m_ConsensusColVec[0];
	fprintf(f, "\n");
	fprintf(f, "=================================================\n");
	fprintf(f, "Motif A\n");
	fprintf(f, "=================================================\n");
	fprintf(f, "Consensus column for motif A catalytic D =");
	if (Col1 == UINT_MAX)
		fprintf(f, " *\n");
	else
		 fprintf(f, " %u\n", Col1+1);

	if (Col1 == UINT_MAX)
		return;

	const uint MARGIN = 20;
	uint LoCol = (Col1 < 20 ? 0 : Col1 - 20);
	uint HiCol = Col1 + 20;
	if (HiCol >= m_ColCount)
		HiCol = m_ColCount - 1;
	asserta(LoCol <= HiCol);
	fprintf(f, "\n");
	uint A1BadCount = 0;
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		asserta(SeqIndex < SIZE(m_AlignedSeqs));
		uint SeqCol1 = m_SeqIndexToColVec[SeqIndex][0];
		bool A1Bad = false;
		if (SeqCol1 != UINT_MAX && SeqCol1 != Col1)
			{
			A1Bad = true;
			++A1BadCount;
			}

		const string &Seq = m_AlignedSeqs[SeqIndex];
		for (uint Col = LoCol; Col <= HiCol; ++Col)
			{
			asserta(Col < SIZE(Seq));
			char c = Seq[Col];
			asserta(SeqIndex < SIZE(m_SeqIndexToColVec));
			asserta(SIZE(m_SeqIndexToColVec[SeqIndex]) == NPOS);

			if (Col == Col1)
				fprintf(f, " |");
			fprintf(f, "%c", c);
			if (Col == Col1)
				fprintf(f, "| ");
			}
		const char *Label = m_MSA->GetLabel(SeqIndex).c_str();
		if (A1Bad)
			fprintf(f, "  ** disagree");
		fprintf(f, "  >%s", Label);
		fprintf(f, "\n");
		}
	fprintf(f, "\n");
	fprintf(f, "Total disagreements motif A %u (%.1f%%)\n",
	  A1BadCount, GetPct(A1BadCount, m_SeqCount));
	}

void MSAQC2::WriteAlnB(FILE *f) const
	{
	if (f == 0)
		return;
	if (m_ConsensusColVec.empty())
		return;
	asserta(SIZE(m_ConsensusColVec) == NPOS);
	uint Col1 = m_ConsensusColVec[1];
	fprintf(f, "\n");
	fprintf(f, "=================================================\n");
	fprintf(f, "Motif B\n");
	fprintf(f, "=================================================\n");
	fprintf(f, "Consensus column for motif B catalytic G =");
	if (Col1 == UINT_MAX)
		fprintf(f, " *\n");
	else
		 fprintf(f, " %u\n", Col1+1);

	if (Col1 == UINT_MAX)
		return;

	const uint MARGIN = 20;
	uint LoCol = (Col1 < 20 ? 0 : Col1 - 20);
	uint HiCol = Col1 + 20;
	if (HiCol >= m_ColCount)
		HiCol = m_ColCount - 1;
	asserta(LoCol <= HiCol);
	fprintf(f, "\n");
	uint B1BadCount = 0;
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		asserta(SeqIndex < SIZE(m_AlignedSeqs));
		uint SeqCol1 = m_SeqIndexToColVec[SeqIndex][1];
		bool B1Bad = false;
		if (SeqCol1 != UINT_MAX && SeqCol1 != Col1)
			{
			B1Bad = true;
			++B1BadCount;
			}

		const string &Seq = m_AlignedSeqs[SeqIndex];
		for (uint Col = LoCol; Col <= HiCol; ++Col)
			{
			asserta(Col < SIZE(Seq));
			char c = Seq[Col];
			asserta(SeqIndex < SIZE(m_SeqIndexToColVec));
			asserta(SIZE(m_SeqIndexToColVec[SeqIndex]) == NPOS);

			if (Col == Col1)
				fprintf(f, " |");
			fprintf(f, "%c", c);
			if (Col == Col1)
				fprintf(f, "| ");
			}
		const char *Label = m_MSA->GetLabel(SeqIndex).c_str();
		if (B1Bad)
			fprintf(f, "  ** disagree **");
		fprintf(f, "  >%s", Label);
		fprintf(f, "\n");
		}
	fprintf(f, "\n");
	fprintf(f, "Total disagreements motif B %u (%.1f%%)\n",
	  B1BadCount, GetPct(B1BadCount, m_SeqCount));
	}

void MSAQC2::WriteAlnC(FILE *f) const
	{
	if (f == 0)
		return;
	if (m_ConsensusColVec.empty())
		return;
	asserta(SIZE(m_ConsensusColVec) == NPOS);
	uint Col1 = m_ConsensusColVec[2];
	fprintf(f, "\n");
	fprintf(f, "=================================================\n");
	fprintf(f, "Motif C\n");
	fprintf(f, "=================================================\n");
	fprintf(f, "Consensus column for motif C catalytic D =");
	if (Col1 == UINT_MAX)
		fprintf(f, " *\n");
	else
		 fprintf(f, " %u\n", Col1+1);

	if (Col1 == UINT_MAX)
		return;

	const uint MARGIN = 20;
	uint LoCol = (Col1 < 20 ? 0 : Col1 - 20);
	uint HiCol = Col1 + 20;
	if (HiCol >= m_ColCount)
		HiCol = m_ColCount - 1;
	asserta(LoCol <= HiCol);
	fprintf(f, "\n");
	uint C1BadCount = 0;
	uint C2BadCount = 0;
	uint C3BadCount = 0;
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		asserta(SeqIndex < SIZE(m_AlignedSeqs));
		uint SeqCol1 = m_SeqIndexToColVec[SeqIndex][2];
		bool C1Bad = false;
		if (SeqCol1 != UINT_MAX && SeqCol1 != Col1)
			{
			C1Bad = true;
			++C1BadCount;
			}

		const string &Seq = m_AlignedSeqs[SeqIndex];
		for (uint Col = LoCol; Col <= HiCol; ++Col)
			{
			asserta(Col < SIZE(Seq));
			char c = Seq[Col];
			asserta(SeqIndex < SIZE(m_SeqIndexToColVec));
			asserta(SIZE(m_SeqIndexToColVec[SeqIndex]) == NPOS);

			if (Col == Col1)
				fprintf(f, " |");
			fprintf(f, "%c", c);
			if (Col == Col1)
				fprintf(f, "| ");
			}
		const char *Label = m_MSA->GetLabel(SeqIndex).c_str();
		if (C1Bad)
			fprintf(f, " ** disagree **");
		fprintf(f, "  >%s", Label);
		fprintf(f, "\n");
		}
	fprintf(f, "\n");
	fprintf(f, "Total disagreements motif C %u (%.1f%%)\n",
	  C1BadCount, GetPct(C1BadCount, m_SeqCount));
	}

void MSAQC2::WriteFev(FILE *f, const string &Name) const
	{
	if (f == 0)
		return;

	vector<bool> IsMs;
	uint MatchStateCount = GetColIsMatchStateVec(IsMs);

	asserta(SIZE(m_ConsensusColVec) == 3);
	uint ColA = m_ConsensusColVec[0];
	uint ColB = m_ConsensusColVec[1];
	uint ColC = m_ConsensusColVec[2];

	const uint ColCount = SIZE(IsMs);
	uint MatchStateA = UINT_MAX;
	uint MatchStateB = UINT_MAX;
	uint MatchStateC = UINT_MAX;
	uint M = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (!IsMs[Col])
			continue;

		if (Col == ColA)
			{
			asserta(MatchStateA == UINT_MAX);
			MatchStateA = M;
			}
		else if (Col == ColB)
			{
			asserta(MatchStateB == UINT_MAX);
			MatchStateB = M;
			}
		else if (Col == ColC)
			{
			asserta(MatchStateC == UINT_MAX);
			MatchStateC = M;
			}
		++M;
		}

	fprintf(g_ffev, "%s", Name.c_str());
	fprintf(g_ffev, "\tM=%u", MatchStateCount);
	if (MatchStateA != UINT_MAX)
		fprintf(g_ffev, "\tMA=%u", MatchStateA);

	if (MatchStateB != UINT_MAX)
		fprintf(g_ffev, "\tMB=%u", MatchStateB);

	if (MatchStateC != UINT_MAX)
		fprintf(g_ffev, "\tMC=%u", MatchStateC);

	fprintf(g_ffev, "\n");
	}

void cmd_msaqc2()
	{
	const string &InputFileName = opt_msaqc2;
	FILE *fOut = CreateStdioFile(opt_output);
	MSAQC2 QC;
	SeqDB MSA1;
	MSA1.FromFasta(InputFileName, true);
	asserta(MSA1.IsAligned());

	SeqDB MSA;
	const uint SeqCount1 = MSA1.GetSeqCount();
	uint ConsensusCount = 0;
	for (uint i = 0; i < SeqCount1; ++i)
		{
		const string Label = string(MSA1.GetLabel(i));
		if (Label == "CONSENSUS")
			{
			++ConsensusCount;
			continue;
			}
		const string &Seq = MSA1.GetSeq(i);
		MSA.AddSeq(Label, Seq);
		}
	if (ConsensusCount > 0)
		ProgressLog("%u CONSENSUS seqs removed\n", ConsensusCount);

	QC.Init(MSA);

	uint CorrectSeqCount = QC.GetCorrectSeqCount();
	uint CorrectLetterCount;
	uint TotalLetterCount;
	QC.GetCorrectLetterCount(CorrectLetterCount, TotalLetterCount);
	if (TotalLetterCount == 0)
		{
		fprintf(fOut, "None found\n");
		CloseStdioFile(fOut);
		ProgressLog("None found");
		return;
		}
	asserta(TotalLetterCount > 0);
	asserta(CorrectLetterCount <= TotalLetterCount);
	uint WrongLetterCount = TotalLetterCount - CorrectLetterCount;
	double WrongLetterPct = GetPct(WrongLetterCount, TotalLetterCount);
	uint SeqCount = QC.m_SeqCount;
	uint DisagreeSeqCount = SeqCount - CorrectSeqCount;
	double Pct = GetPct(DisagreeSeqCount, SeqCount);
	double AgreePct = GetPct(CorrectSeqCount, SeqCount);
	if (fOut != 0)
		{
		QC.WriteAlnA(fOut);
		QC.WriteAlnB(fOut);
		QC.WriteAlnC(fOut);

		fprintf(fOut, "@QC2 %s:", InputFileName.c_str());
		fprintf(fOut, "  %u/%u seqs (%.2f%%)",
		  DisagreeSeqCount, SeqCount, Pct);
		fprintf(fOut, ", %u/%u letters disagree (%.2f%%)\n",
		  WrongLetterCount, TotalLetterCount, WrongLetterPct);

		CloseStdioFile(fOut);
		fOut = 0;
		}

	QC.WriteFev(g_ffev, InputFileName);

	ProgressLog("QC2 %s", InputFileName.c_str());
	ProgressLog("  %u/%u seqs (%.2f%%)",
	  DisagreeSeqCount, SeqCount, Pct);
	ProgressLog(", %u/%u letters disagree (%.2f%%)\n",
	  WrongLetterCount, TotalLetterCount, WrongLetterPct);
	}

bool MSAQC2::GetColIsMatchState(uint ColIndex) const
	{
	const uint SeqCount = m_MSA->GetSeqCount();
	uint GapCount = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Seq = m_MSA->GetSeq(i);
		asserta(ColIndex < SIZE(Seq));
		char c = Seq[ColIndex];
		if (isgap(c))
			++GapCount;
		}
	double Yes = double(GapCount)/double(SeqCount) < 0.5;
	return Yes;
	}

uint MSAQC2::GetColIsMatchStateVec(vector<bool> &v) const
	{
	v.clear();
	uint n = 0;
	const uint ColCount = m_MSA->GetColCount();
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		bool Yes = GetColIsMatchState(ColIndex);
		if (Yes)
			++n;
		v.push_back(Yes);
		}
	return n;
	}

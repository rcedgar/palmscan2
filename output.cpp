#include "myutils.h"
#include "rdrpsearcher.h"
#include "outputfiles.h"

static omp_lock_t g_OutputLock;

void RdRpSearcher::InitOutput()
	{
	omp_init_lock(&g_OutputLock);
	}

void RdRpSearcher::WriteOutput() const
	{
#pragma omp critical
	{
	WriteReport(g_freport_pssms);
	WriteTsv(g_ftsv);
	}
	}

void RdRpSearcher::WriteReport(FILE *f) const
	{
	if (f == 0)
		return;

	if (m_TopPalmHit.m_Score <= 0)
		return;

	uint GroupIndex = m_TopPalmHit.m_GroupIndex;
	bool Permuted = m_TopPalmHit.m_Permuted;
	float Score = m_TopPalmHit.m_Score;

	string GroupName;
	GetGroupName(GroupIndex, GroupName);

	vector<string> Rows;
	GetAlnRows(Rows);

	fprintf(f, "\n");
	fprintf(f, ">%s\n", m_QueryLabel.c_str());
	for (uint i = 0; i < SIZE(Rows); ++i)
		fprintf(f, "%s\n", Rows[i].c_str());

	fprintf(f, "Score %.1f", Score);
	if (Permuted)
		fprintf(f, " permuted");
	fprintf(f, ", %s", GroupName.c_str());

	float SecondScore = m_SecondPalmHit.m_Score;
	if (SecondScore > 0)
		{
		string SecondGroupName;
		uint SecondGroupIndex = m_SecondPalmHit.m_GroupIndex;
		GetGroupName(SecondGroupIndex, SecondGroupName);
		float Diff = Score - SecondScore;
		fprintf(f, " (%s +%.1f)", SecondGroupName.c_str(), Diff);
		}
	fprintf(f, "\n");
	}

void RdRpSearcher::WriteTsv(FILE *f) const
	{
	if (f == 0)
		return;

	static bool HdrDone = false;
	if (!HdrDone)
		{
		HdrDone = true;

		fprintf(f, "Label");
		fprintf(f, "\tScore");
		fprintf(f, "\tGroup");
		fprintf(f, "\tGroup2");
		fprintf(f, "\tDiff2");
		fprintf(f, "\tABC");
		fprintf(f, "\tQL");
		fprintf(f, "\tLo");
		fprintf(f, "\tHi");
		fprintf(f, "\tPPL");
		fprintf(f, "\tSuff");
		fprintf(f, "\tPosA");
		fprintf(f, "\tSeqA");
		fprintf(f, "\tPosB");
		fprintf(f, "\tSeqB");
		fprintf(f, "\tPosC");
		fprintf(f, "\tSeqC");
		fprintf(f, "\n");
		}

	string QueryLabel = m_QueryLabel;
	float Score = 0;
	string GroupName = ".";
	string SecondGroupName = ".";
	float Diff2 = 0;
	const char *ABC = ".";
	uint QL = SIZE(m_QuerySeq);
	uint Lo = 0;
	uint Hi = 0;
	uint PPL = 0;
	uint Suff = 0;

	if (m_TopPalmHit.m_Score > 0)
		{
		Score = m_TopPalmHit.m_Score;
		uint GroupIndex = m_TopPalmHit.m_GroupIndex;
		GetGroupName(GroupIndex, GroupName);
		if (m_TopPalmHit.m_Permuted)
			ABC = "CAB";
		else
			ABC = "ABC";
		GetSpan(Lo, Hi);
		++Lo;
		++Hi;
		assert(Hi <= QL);
		PPL = Hi - Lo + 1;
		Suff = QL - Hi;
		}

	if (m_SecondPalmHit.m_Score > 0)
		{
		float SecondScore = m_SecondPalmHit.m_Score;
		uint SecondGroupIndex = m_SecondPalmHit.m_GroupIndex;
		GetGroupName(SecondGroupIndex, SecondGroupName);
		Diff2 = Score - SecondScore;
		}

	uint PosA = GetMotifPos(0);
	uint PosB = GetMotifPos(1);
	uint PosC = GetMotifPos(2);

	string SeqA = ".";
	string SeqB = ".";
	string SeqC = ".";

	if (PosA == UINT_MAX)
		PosA = 0;
	else
		{
		GetMotifSeq(0, SeqA);
		++PosA;
		}
	if (PosB == UINT_MAX)
		PosB = 0;
	else
		{
		GetMotifSeq(1, SeqB);
		++PosB;
		}
	if (PosC == UINT_MAX)
		PosC = 0;
	else
		{
		GetMotifSeq(2, SeqC);
		++PosC;
		}

	fprintf(f, "%s", QueryLabel.c_str());
	fprintf(f, "\t%.1f", Score);
	fprintf(f, "\t%s", GroupName.c_str());
	fprintf(f, "\t%s", SecondGroupName.c_str());
	fprintf(f, "\t%+.1f", Diff2);
	fprintf(f, "\t%s", ABC);
	fprintf(f, "\t%u", QL);
	fprintf(f, "\t%u", Lo);
	fprintf(f, "\t%u", Hi);
	fprintf(f, "\t%u", PPL);
	fprintf(f, "\t%u", Suff);
	fprintf(f, "\t%u", PosA);
	fprintf(f, "\t%s", SeqA.c_str());
	fprintf(f, "\t%u", PosB);
	fprintf(f, "\t%s", SeqB.c_str());
	fprintf(f, "\t%u", PosC);
	fprintf(f, "\t%s", SeqC.c_str());
	fprintf(f, "\n");
	}

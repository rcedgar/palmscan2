#include "myutils.h"
#include "rdrpsearcher.h"

static omp_lock_t g_OutputLock;

static FILE *g_fRep = 0;
static FILE *g_fTsv = 0;
static FILE *g_fTri = 0;

void RdRpSearcher::InitOutput()
	{
	omp_init_lock(&g_OutputLock);
	g_fRep = CreateStdioFile(opt_report);
	g_fTsv = CreateStdioFile(opt_tsvout);
	g_fTri = CreateStdioFile(opt_triout);
	}

void RdRpSearcher::CloseOutput()
	{
	CloseStdioFile(g_fRep);
	CloseStdioFile(g_fTsv);
	CloseStdioFile(g_fTri);
	}

void RdRpSearcher::WriteOutput() const
	{
	omp_set_lock(&g_OutputLock);
	WriteReport(g_fRep);
	WriteTsv(g_fTsv);
	WriteTri(g_fTri);
	omp_unset_lock(&g_OutputLock);
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
		Suff = QL - Hi;
		}

	if (m_SecondPalmHit.m_Score > 0)
		{
		float SecondScore = m_SecondPalmHit.m_Score;
		uint SecondGroupIndex = m_SecondPalmHit.m_GroupIndex;
		GetGroupName(SecondGroupIndex, SecondGroupName);
		Diff2 = Score - SecondScore;
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
	fprintf(f, "\t%u", Suff);
	fprintf(f, "\n");
	}

void RdRpSearcher::WriteTri(FILE *f) const
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
		fprintf(f, "\taD");
		fprintf(f, "\tbG");
		fprintf(f, "\tcD");
		fprintf(f, "\tPosaD");
		fprintf(f, "\tPosbG");
		fprintf(f, "\tPoscD");
		fprintf(f, "\n");
		}

	string QueryLabel = m_QueryLabel;
	float Score = 0;
	string GroupName = ".";
	char aD = '.';
	char bG = '.';
	char cD = '.';
	uint PosaD = 0;
	uint PosbG = 0;
	uint PoscD = 0;

	if (m_TopPalmHit.m_Score > 0)
		{
		Score = m_TopPalmHit.m_Score;
		uint GroupIndex = m_TopPalmHit.m_GroupIndex;
		GetGroupName(GroupIndex, GroupName);
		}

	GetTri(aD, bG, cD, PosaD, PosbG, PoscD);

	fprintf(f, "%s", QueryLabel.c_str());
	fprintf(f, "\t%.1f", Score);
	fprintf(f, "\t%s", GroupName.c_str());
	fprintf(f, "\t%c", aD);
	fprintf(f, "\t%c", bG);
	fprintf(f, "\t%c", cD);
	fprintf(f, "\t%u", PosaD);
	fprintf(f, "\t%u", PosbG);
	fprintf(f, "\t%u", PoscD);
	fprintf(f, "\n");
	}

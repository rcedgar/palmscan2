#include "myutils.h"
#include "ppsp.h"
#include "abcxyz.h"

uint PPSP::GetSeqPos(uint i, uint APos, uint BPos, uint CPos)
	{
	if (i >= CIX)
		return CPos + i - (CIX);
	if (i >= AL)
		return BPos + i - AL;
	return APos + i;
	}

char PPSP::GetMotifChar(uint Ix) const
	{
	if (Ix < AL)
		return 'A';
	if (Ix < AL + BL)
		return 'B';
	if (Ix < PPSPL)
		return 'C';
	asserta(false);
	return '!';
	}

void PPSP::GetDistMx(const PDBChain &Q,
  uint APos, uint BPos, uint CPos,
  vector<vector<double> > &DistMx)
	{
	asserta(APos != UINT_MAX);
	asserta(BPos != UINT_MAX);
	asserta(CPos != UINT_MAX);

	DistMx.clear();
	const uint N = CIX + CL;
	DistMx.resize(N);
	for (uint i = 0; i < N; ++i)
		DistMx[i].resize(N, DBL_MAX);

	for (uint i = 0; i < N; ++i)
		{
		uint SeqPosi = GetSeqPos(i, APos, BPos, CPos);
		for (uint j = 0; j < N; ++j)
			{
			uint SeqPosj = GetSeqPos(j, APos, BPos, CPos);
			double d = Q.GetDist(SeqPosi, SeqPosj);
			DistMx[i][j] = d;
			}
		}
	}

void PPSP::ToFile(const string &FileName) const
	{
	if (FileName == "")
		return;
	asserta(SIZE(m_Means) == PPSPL);
	asserta(SIZE(m_StdDevs) == PPSPL);
	FILE *f = CreateStdioFile(FileName);
	asserta(f != 0);
	fprintf(f, "PPSP\n");
	for (uint i = 1; i < PPSPL; ++i)
		{
		fprintf(f, "mean\t%u", i);
		for (uint j = 0; j <= i; ++j)
			{
			double Mean = m_Means[i][j];
			if (i == j)
				asserta(Mean == 0);
			asserta(m_Means[j][i] == m_Means[i][j]);
			fprintf(f, "\t%.3g", Mean);
			}
		fprintf(f, "\n");
		}

	for (uint i = 1; i < PPSPL; ++i)
		{
		fprintf(f, "stddev\t%u", i);
		for (uint j = 0; j <= i; ++j)
			{
			double Mean = m_Means[i][j];
			double StdDev = m_StdDevs[i][j];
			if (i == j)
				asserta(StdDev == 0);
			asserta(m_StdDevs[j][i] == m_StdDevs[i][j]);
			fprintf(f, "\t%.3g", StdDev);
			}
		fprintf(f, "\n");
		}

	CloseStdioFile(f);
	}

void PPSP::FromFile(const string &FileName)
	{
	asserta(FileName != "");
	Clear();
	string Line;
	vector<string> Fields;
	FILE *f = OpenStdioFile(FileName);
	bool Ok = ReadLineStdioFile(f, Line);
	if (!Ok)
		Die("Premature EOF in PPSP file %s", FileName.c_str());
	if (Line != "PPSP")
		Die("Not PSPP file %s", FileName.c_str());

	for (uint i = 1; i < PPSPL; ++i)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			Die("Premature EOF in PPSP file %s", FileName.c_str());
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == i+3);
		asserta(Fields[0] == "mean");
		asserta(StrToUint(Fields[1]) == i);
		for (uint j = 0; j <= i; ++j)
			{
			double Mean = StrToFloat(Fields[j+2]);
			if (i == j)
				asserta(Mean == 0);
			m_Means[i][j] = Mean;
			m_Means[j][i] = m_Means[i][j];
			}
		}

	for (uint i = 1; i < PPSPL; ++i)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			Die("Premature EOF in PPSP file %s", FileName.c_str());
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == i+3);
		asserta(Fields[0] == "stddev");
		asserta(StrToUint(Fields[1]) == i);
		for (uint j = 0; j <= i; ++j)
			{
			double StdDev = StrToFloat(Fields[j+2]);
			if (i == j)
				asserta(StdDev == 0);
			m_StdDevs[i][j] = StdDev;
			m_StdDevs[j][i] = m_StdDevs[i][j];
			}
		}
	}

double PPSP::GetScoreA(const PDBChain &Chain, uint PosA) const
	{
	asserta(Chain.m_Seq[PosA+3] == 'D');
	double Score = GetScore(Chain, PosA, AIX, AL);
	return Score;
	}

double PPSP::GetScoreB(const PDBChain &Chain, uint PosB) const
	{
	asserta(Chain.m_Seq[PosB+1] == 'G');
	double Score = GetScore(Chain, PosB, BIX, BL);
	return Score;
	}

double PPSP::GetScoreC(const PDBChain &Chain, uint PosC) const
	{
	asserta(Chain.m_Seq[PosC+3] == 'D');
	double Score = GetScore(Chain, PosC, CIX, CL);
	return Score;
	}

double PPSP::GetScore(const PDBChain &Chain, uint SeqPos,
  uint Ix, uint L, bool Trace) const
	{
	double Score = GetScore2(Chain, SeqPos, SeqPos, Ix, Ix, L, L, Trace);
	return Score;
	}

double PPSP::GetScore2(const PDBChain &Chain,
  uint SeqPos1, uint SeqPos2,
  uint Ix1, uint Ix2,
  uint L1, uint L2, bool Trace) const
	{
	bool Diag = (Ix1 == Ix2);
	double XS = (Diag ? 1.5 : 2);
	double Sum = 0;
	uint n = 0;
	if (Trace)
		Log("\n");
	for (uint i = 0; i < L1; ++i)
		{
		uint jhi = (Diag ? i : L2);
		for (uint j = 0; j < jhi; ++j)
			{
			double Observed_d = Chain.GetDist(SeqPos1+i, SeqPos2+j);
			double Mu = m_Means[Ix1+i][Ix2+j];
			double Sigma = m_StdDevs[Ix1+i][Ix2+j];
			double y = GetNormal(Mu, XS*Sigma, Observed_d);
			double Max = GetNormal(Mu, XS*Sigma, Mu);
			double Ratio = y/Max;
			Sum += Ratio;
			++n;

			if (Trace)
				{
				char MotifChari = GetMotifChar(Ix1);
				char MotifCharj = GetMotifChar(Ix2);

				char ci = Chain.m_Seq[SeqPos1+i];
				char cj = Chain.m_Seq[SeqPos2+j];

				Log("%c[%3u]%c", MotifChari, SeqPos1+i, ci);
				Log("  %c[%3u]%c", MotifCharj, SeqPos2+j, cj);
				Log("  d %6.2f", Observed_d);
				Log("  mu %6.2f", Mu);
				Log("  sd %6.2f", Sigma);
				Log("  r %6.4f", Ratio);
				Log("\n");
				}
			}
		}
	asserta(n > 0);
	double Score = Sum/n;
	return Score;
	}

double PPSP::GetScoreAB(const PDBChain &Chain, uint PosA, uint PosB) const
	{
	double ScoreAB = GetScore2(Chain, PosA, PosB, AIX, BIX, AL, BL);
	return ScoreAB;
	}

double PPSP::GetScoreBC(const PDBChain &Chain, uint PosB, uint PosC) const
	{
	double ScoreBC = GetScore2(Chain, PosB, PosC, BIX, CIX, BL, CL);
	return ScoreBC;
	}

double PPSP::GetScoreAC(const PDBChain &Chain, uint PosA, uint PosC) const
	{
	double ScoreAC = GetScore2(Chain, PosA, PosC, AIX, CIX, AL, CL);
	return ScoreAC;
	}

double PPSP::GetScore3(const PDBChain &Chain,
  uint PosA, uint PosB, uint PosC) const
	{
	double Score = 0;

	Score += GetScoreA(Chain, PosA);
	Score += GetScoreB(Chain, PosB);
	Score += GetScoreC(Chain, PosC);

	Score += GetScoreAB(Chain, PosA, PosB);
	Score += GetScoreBC(Chain, PosB, PosC);
	Score += GetScoreAC(Chain, PosA, PosC);

	Score /= 6;
	return Score;
	}

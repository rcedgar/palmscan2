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

bool PPSP::GetDistMx(const PDBChain &Q,
  uint APos, uint BPos, uint CPos,
  vector<vector<double> > &DistMx)
	{
	DistMx.clear();
	if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
		return false;

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
	return true;
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

double PPSP::GetScore(const PDBChain &Chain, uint SeqPos, uint Ix, uint L) const
	{
	double Sum = 0;
	uint n = 0;
	for (uint i = 1; i < L; ++i)
		{
		for (uint j = 0; j < i; ++j)
			{
			double Observed_d = Chain.GetDist(SeqPos+i, SeqPos+j);
			double Mu = m_Means[Ix+i][Ix+j];
			double Sigma = m_StdDevs[Ix+i][Ix+j];
			double y = GetNormal(Mu, Sigma, Observed_d);
			double Max = GetNormal(Mu, Sigma, Mu);
			double Ratio = y/Max;
			Sum += Ratio;
			++n;
			}
		}
	asserta(n > 0);
	double Score = Sum/n;
	return Score;
	}

double PPSP::GetScore2(const PDBChain &Chain,
  uint SeqPos1, uint SeqPos2,
  uint Ix1, uint Ix2,
  uint L1, uint L2) const
	{
	double Sum = 0;
	uint n = 0;
	for (uint i = 0; i < L1; ++i)
		{
		for (uint j = 0; j < L2; ++j)
			{
			double Observed_d = Chain.GetDist(SeqPos1+i, SeqPos2+j);
			double Mu = m_Means[Ix1+i][Ix2+j];
			double Sigma = m_StdDevs[Ix1+i][Ix2+j];
			double y = GetNormal(Mu, Sigma, Observed_d);
			double Max = GetNormal(Mu, Sigma, Mu);
			double Ratio = y/Max;
			Sum += Ratio;
			++n;
			}
		}
	asserta(n > 0 && n == L1*L2);
	double Score = Sum/n;
	return Score;
	}

double PPSP::GetScore3(const PDBChain &Chain,
  uint PosA, uint PosB, uint PosC) const
	{
	double ScoreA = GetScoreA(Chain, PosA);
	double ScoreB = GetScoreB(Chain, PosB);
	double ScoreC = GetScoreC(Chain, PosC);

	double ScoreAB = GetScore2(Chain, PosA, PosB, AIX, BIX, AL, BL);
	double ScoreBC = GetScore2(Chain, PosB, PosC, BIX, CIX, BL, CL);
	double ScoreAC = GetScore2(Chain, PosA, PosC, AIX, CIX, AL, CL);

	double Score = 0;

	Score += ScoreA;
	Score += ScoreB;
	Score += ScoreC;

	Score += ScoreAB;
	Score += ScoreBC;
	Score += ScoreAC;

	Score /= 6;

	return Score;
	}

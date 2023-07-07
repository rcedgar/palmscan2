#include "myutils.h"
#include "shapesearcher.h"
#include "motifsettings.h"
#include "sort.h"

#define TRACE	0

void ShapeSearcher::Init(const Shapes &S)
	{
	m_Shapes = &S;
	m_ShapeCount = S.GetShapeCount();
	m_ShapeIndexA = UINT_MAX;
	m_ShapeIndexB = UINT_MAX;
	m_ShapeIndexC = UINT_MAX;
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		if (S.m_Names[i] == "A")
			m_ShapeIndexA = i;
		else if (S.m_Names[i] == "B")
			m_ShapeIndexB = i;
		else if (S.m_Names[i] == "C")
			m_ShapeIndexC = i;
		}
	if (optset_minscore_pp)
		m_MinScoreABC = opt_minscore_pp;
	}

void ShapeSearcher::SetQuery(const PDBChain &Query)
	{
	ClearSearch();
	m_Query = &Query;

	const uint ShapeCount = GetShapeCount();
	m_ShapePosVec.resize(ShapeCount);
	m_ShapeScores.resize(ShapeCount);
	for (uint i = 0; i < ShapeCount; ++i)
		{
		m_ShapePosVec[i] = UINT_MAX;
		m_ShapeScores[i] = 0;
		}
	}

void ShapeSearcher::GetFoundMotifsStr(string &Str) const
	{
	Str.resize(0);
	if (SIZE(m_ShapePosVec) != m_ShapeCount)
		{
		Str = "*";
		return;
		}
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		bool Found = (m_ShapePosVec[i] != UINT_MAX);
		const string Name = (string) GetShapeName(i);
		if (Found)
			Str += Name;
		else
			{
			Str += "-";
			for (uint j = 1; j < SIZE(Name); ++j)
				Str += ".";
			}
		}

	if (m_Permuted)
		{
		size_t n = Str.find("ABC");
		if (n == string::npos)
			Str += "-permuted";
		else
			{
			Str[n] = 'C';
			Str[n+1] = 'A';
			Str[n+2] = 'B';
			}
		}
	}

char ShapeSearcher::GetGate() const
	{
	if (m_Query == 0)
		return '.';
	if (m_ShapeIndexA == UINT_MAX)
		return '.';
	asserta(m_ShapeIndexA < SIZE(m_ShapePosVec));
	if (m_PosA == UINT_MAX)
		return '.';
	const string &Seq = m_Query->m_Seq;
	uint Pos = m_PosA + g_OffAd + 5;
	if (Pos >= SIZE(Seq))
		return '.';
	char Gate = Seq[Pos];
	return Gate;
	}

void ShapeSearcher::GetSubSeq(uint Pos, uint n, string &Seq) const
	{
	m_Query->GetSubSeq(Pos, n, Seq);
	}

void ShapeSearcher::GetA(string &A) const
	{
	GetShapeSeq(m_ShapeIndexA, A);
	}

void ShapeSearcher::GetB(string &B) const
	{
	GetShapeSeq(m_ShapeIndexB, B);
	}

void ShapeSearcher::GetC(string &C) const
	{
	GetShapeSeq(m_ShapeIndexC, C);
	}

void ShapeSearcher::GetShapeSeq(uint ShapeIndex, string &Seq) const
	{
	Seq = ".";
	if (ShapeIndex == UINT_MAX)
		return;
	uint ShapeLength = GetShapeLength(ShapeIndex);
	asserta(ShapeIndex < SIZE(m_ShapePosVec));
	uint Pos = m_ShapePosVec[ShapeIndex];
	m_Query->GetSubSeq(Pos, ShapeLength, Seq);
	}

// dist = start_j - start_i
void ShapeSearcher::GetDistRange(uint ShapeIndexi, uint ShapeIndexj,
  uint &MinDist, uint &MaxDist) const
	{
	asserta(ShapeIndexj > ShapeIndexi);

	MinDist = 0;
	MaxDist = 0;
	for (uint ShapeIndex = ShapeIndexi; ShapeIndex < ShapeIndexj;
	  ++ShapeIndex)
		{
		MinDist += m_Shapes->m_MinNeighborDists[ShapeIndex];
		MaxDist += m_Shapes->m_MaxNeighborDists[ShapeIndex];
		}
	}

double ShapeSearcher::GetScoreShapes(const vector<uint> &ShapeIndexes,
  const vector<uint> &PosVec) const
	{
	double Sum = 0;
	const uint N = SIZE(ShapeIndexes);
	if (N == 0)
		return 0;
	uint Pairs = 0;
	uint FoundCount = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint ShapeIndexi = ShapeIndexes[i];
		uint Posi = PosVec[i];
		if (Posi != UINT_MAX)
			++FoundCount;
		for (uint j = 0; j <= i; ++j)
			{
			uint ShapeIndexj = ShapeIndexes[j];
			uint Posj = PosVec[j];
			Sum += GetScoreShapePair(ShapeIndexi, ShapeIndexj, Posi, Posj);
			++Pairs;
			}
		}
	double Avg = Sum/Pairs;
//	double Score = Avg*double(FoundCount)/N;
	double Score = Avg;
	return Score;
	}

double ShapeSearcher::GetScoreShapePair(uint ShapeIndex1, uint ShapeIndex2,
  uint Pos1, uint Pos2) const
	{
	if (ShapeIndex1 == UINT_MAX || ShapeIndex2 == UINT_MAX)
		return 0;
	if (Pos1 == UINT_MAX || Pos2 == UINT_MAX)
		return 0;
	double Sum = 0;
	uint L1 = GetShapeLength(ShapeIndex1);
	uint L2 = GetShapeLength(ShapeIndex2);
	uint QL = GetQL();
	if (Pos1 + L1 > QL)
		return 0;
	if (Pos2 + L2 > QL)
		return 0;
	uint n = 0;
	for (uint Offset1 = 0; Offset1 < L1; ++Offset1)
		{
		uint StartOffset2 = (ShapeIndex1 == ShapeIndex2 ? Offset1 + 1 : 0);
		for (uint Offset2 = StartOffset2; Offset2 < L2; ++Offset2)
			{
			Sum += GetScoreResiduePair(ShapeIndex1, ShapeIndex2,
			  Pos1, Pos2, Offset1, Offset2);
			++n;
			}
		}
	assert(n > 0);
	double Score = Sum/n;
	return Score;
	}

double ShapeSearcher::GetScoreResiduePair(uint ShapeIndex1, uint ShapeIndex2,
  uint Pos1, uint Pos2, uint Offset1, uint Offset2) const
	{
	bool Diag = (ShapeIndex1 == ShapeIndex2);
	double XS = (Diag ? 1.5 : 2);
	const PDBChain &Query = *m_Query;
	double Observed_d = Query.GetDist(Pos1+Offset1, Pos2+Offset2);
	Observed_d = fabs(Observed_d);
	double Mu = GetMeanDist3(ShapeIndex1, ShapeIndex2, Offset1, Offset2);
	Mu = fabs(Mu);
	double Sigma = GetStdDev3(ShapeIndex1, ShapeIndex2, Offset1, Offset2);
	double y = GetNormal(Mu, XS*Sigma, Observed_d);
	double Max = GetNormal(Mu, XS*Sigma, Mu);
	double Score = y/Max;
	return Score;
	}

double ShapeSearcher::GetSelfScore(uint ShapeIndex, uint Pos) const
	{
	double Score = GetScoreShapePair(ShapeIndex, ShapeIndex, Pos, Pos);
	return Score;
	}

double ShapeSearcher::GetScore(uint ShapeIndex, uint Pos,
  const vector<uint> &PosVec) const
	{
	const uint ShapeCount = GetShapeCount();
	asserta(SIZE(PosVec) == ShapeCount);
	double Sum = 0;
	uint n = 0;
	for (uint ShapeIndex2 = 0; ShapeIndex2 < ShapeCount; ++ShapeIndex2)
		{
		uint Pos2 = PosVec[ShapeIndex2];
		if (Pos2 == UINT_MAX)
			continue;
		++n;
		Sum += GetScoreShapePair(ShapeIndex, ShapeIndex2, Pos, Pos2);
		}
	asserta(n > 0);
	double Score = Sum/n;
	return Score;
	}

void ShapeSearcher::SearchShapeSelfTop(uint ShapeIndex,
  double MinScore, uint MaxHits, vector<uint> &HitPosVec) const
	{
	HitPosVec.clear();
	uint QL = GetQL();
	uint Length = GetShapeLength(ShapeIndex);
	if (QL < Length + 10)
		return;

	vector<uint> TmpHitPosVec;
	vector<double> TmpHitScores;
	SearchShapeSelf(ShapeIndex, m_MinSelfScore, 0, QL-Length, 0, UINT_MAX,
	  TmpHitPosVec, TmpHitScores);
	const uint N = SIZE(TmpHitPosVec);
	asserta(SIZE(TmpHitScores) == N);
	if (N == 0)
		return;
	vector<uint> Order(N);
	QuickSortOrderDesc(TmpHitScores.data(), N, Order.data());
	uint M = min(N, MaxHits);
	for (uint k = 0; k < M; ++k)
		{
		uint i = Order[k];
		HitPosVec.push_back(TmpHitPosVec[i]);
		}
	}

double ShapeSearcher::SearchShapeSelf(uint ShapeIndex, double MinScore,
  uint Lo, uint Hi, char Letter, uint LetterOffset,
  vector<uint> &HitPosVec, vector<double> &HitScores) const
	{
	if (Lo == UINT_MAX || Hi == UINT_MAX)
		return 0;
	HitPosVec.clear();
	HitScores.clear();
	uint L = GetShapeLength(ShapeIndex);
	uint QL = GetQL();
	if (Lo + L > QL)
		return 0;
	double TopScore = 0;
	uint Top = Hi;
	if (Top + L > QL)
		{
		if (QL < L)
			return 0;
		Top = QL - L;
		}
	for (uint Pos = Lo; Pos <= Top; ++Pos)
		{
		if (Letter != 0 && LetterOffset != UINT_MAX &&
		  m_Query->m_Seq[Pos+LetterOffset] != Letter)
			continue;
		double Score = GetSelfScore(ShapeIndex, Pos);
		if (Score < MinScore)
			continue;
		if (Score > TopScore)
			TopScore = Score;
		uint HitCount = SIZE(HitPosVec);
		if (HitCount > 0)
			{
			uint LastHit = HitPosVec[HitCount-1];
			if (Pos - LastHit < L)
				{
				double LastScore = HitScores[HitCount-1];
				if (Score > LastScore)
					{
					HitPosVec[HitCount-1] = Pos;
					HitScores[HitCount-1] = Score;
					}
				continue;
				}
			}
		HitPosVec.push_back(Pos);
		HitScores.push_back(Score);
		}
	return TopScore;
	}

void ShapeSearcher::SearchShapeTopHit(uint ShapeIndex,
  const vector<uint> &PosVec, double MinScore, uint Lo, uint Hi,
  char Letter, uint LetterOffset, uint &Pos, double &Score) const
	{
	Pos = UINT_MAX;
	Score = 0;
	if (Lo == UINT_MAX)
		return;
	asserta(Hi != UINT_MAX);

	Pos = UINT_MAX;
	Score = 0;

	vector<uint> HitPosVec;
	vector<double> Scores;
	SearchShape(ShapeIndex, PosVec, MinScore, Lo, Hi,
	  Letter, LetterOffset, HitPosVec, Scores);

	const uint N = SIZE(HitPosVec);
	for (uint i = 0; i < N; ++i)
		{
		if (Scores[i] > Score)
			{
			Pos = HitPosVec[i];
			Score = Scores[i];
			}
		}
	}

void ShapeSearcher::SearchShape(uint ShapeIndex, const vector<uint> &PosVec,
  double MinScore, uint Lo, uint Hi, char Letter, uint LetterOffset,
  vector<uint> &HitPosVec, vector<double> &HitScores) const
	{
	HitPosVec.clear();
	HitScores.clear();
	if (Lo == UINT_MAX)
		return;
	asserta(Hi != UINT_MAX);
	uint L = GetShapeLength(ShapeIndex);
	uint QL = GetQL();
	//asserta(Hi + L <= QL);
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		if (Letter != 0 && LetterOffset != UINT_MAX &&
		  m_Query->m_Seq[Pos+LetterOffset] != Letter)
			continue;
		double Score = GetScore(ShapeIndex, Pos, PosVec);
		if (Score < MinScore)
			continue;
		uint HitCount = SIZE(HitPosVec);
		if (HitCount > 0)
			{
			uint LastHit = HitPosVec[HitCount-1];
			if (Pos - LastHit < L)
				{
				double LastScore = HitScores[HitCount-1];
				if (Score > LastScore)
					{
					HitPosVec[HitCount-1] = Pos;
					HitScores[HitCount-1] = Score;
					}
				continue;
				}
			}
		HitPosVec.push_back(Pos);
		HitScores.push_back(Score);
		}
	}

void ShapeSearcher::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	static bool HdrDone = false;
#pragma omp critical
	{
	if (!HdrDone)
		{
		HdrDone = true;
		fprintf(f, "Label");
		fprintf(f, "\tClass");
		fprintf(f, "\tPalm_score");
		fprintf(f, "\tPP_score");
		fprintf(f, "\tMotifs");
		fprintf(f, "\tGate");
		for (uint i = 0; i < m_ShapeCount; ++i)
			{
			const char *ShapeName = GetShapeName(i);
			fprintf(f, "\t%s_pos", ShapeName);
			fprintf(f, "\t%s_seq", ShapeName);
			fprintf(f, "\t%s_score", ShapeName);
			}
		fprintf(f, "\n");
		}

	char Gate = GetGate();
	string FoundMotifsStr;
	GetFoundMotifsStr(FoundMotifsStr);

	string Class = ".";
	if (m_ScoreABC >= m_MinScoreABC)
		{
		Class = "PP+";

		switch (Gate)
			{
		case 'D':
		case 'S':
		case 'T':
		case 'G':
		case 'C':
			Class = "RdRp+";
			break;

		case 'Y':
		case 'F':
			Class = "RdRp-";
			break;
			}
		}

	fprintf(f, "%s", m_Query->m_Label.c_str());
	fprintf(f, "\t%s", Class.c_str());
	fprintf(f, "\t%.3g", m_PalmScore);
	fprintf(f, "\t%.3g", m_ScoreABC);
	fprintf(f, "\t%s", FoundMotifsStr.c_str());
	fprintf(f, "\t%c", Gate);
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		double Score = m_ShapeScores[i];
		string Seq;
		GetShapeSeq(i, Seq);
		uint Pos = m_ShapePosVec[i];
		if (Pos == UINT_MAX)
			fprintf(f, "\t.");
		else
			fprintf(f, "\t%u", Pos+1);
		fprintf(f, "\t%s", Seq.c_str());
		fprintf(f, "\t%.3g", Score);
		}
	fprintf(f, "\n");
	}
	}

void ShapeSearcher::ToPmlABC(FILE *f) const
	{
	if (f == 0)
		return;

	if (m_PosA == UINT_MAX || m_PosB == UINT_MAX || m_PosC == UINT_MAX)
		return;

	string SeqA, SeqB, SeqC;
	GetShapeSeq(m_ShapeIndexA, SeqA);
	GetShapeSeq(m_ShapeIndexB, SeqB);
	GetShapeSeq(m_ShapeIndexC, SeqC);

	fprintf(f, "color gray70\n");
	fprintf(f, "select pepseq %s\n", SeqA.c_str());
	fprintf(f, "color tv_blue, sele\n");
	fprintf(f, "select pepseq %s\n", SeqB.c_str());
	fprintf(f, "color tv_green, sele\n");
	fprintf(f, "select pepseq %s\n", SeqC.c_str());
	fprintf(f, "color tv_red, sele\n");
	fprintf(f, "deselect\n");
	}

void ShapeSearcher::ToPml(FILE *f, const string &LoadName) const
	{
	if (f == 0)
		return;

	fprintf(f, "reinitialize;\n");
	if (LoadName != "")
		fprintf(f, "load %s;\n", LoadName.c_str());
	fprintf(f, "color gray70;\n");
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		string Seq, Color;
		GetShapeSeq(i, Seq);
		if (Seq == "" || Seq == ".")
			continue;

		const char *Name = GetShapeName(i);
		GetColor(i, Color);

		fprintf(f, "select Motif_%s, pepseq %s;\n", Name, Seq.c_str());
		fprintf(f, "color %s, Motif_%s;\n", Color.c_str(), Name);
		}
	fprintf(f, "deselect\n");
	if (opt_pml_save)
		{
		const string &SaveName = m_Query->m_Label;
		fprintf(f, "save %s.pse\n", SaveName.c_str());
		}
	}

void ShapeSearcher::GetColor(unsigned MotifIndex, string &Color) const
	{
	string Name = (string) GetShapeName(MotifIndex);
	if (Name == "A")
		Color = "tv_blue";
	else if (Name == "B")
		Color = "tv_green";
	else if (Name == "C")
		Color = "tv_red";
	else if (Name == "F1")
		Color = "deepteal";
	else if (Name == "F2")
		Color = "cyan";
	else if (Name == "D")
		Color = "yellow";
	else if (Name == "E")
		Color = "orange";
	else if (Name == "H")
		Color = "salmon";
	else if (Name == "J")
		Color = "magenta";
	}

void ShapeSearcher::GetShapeIndexes(vector<uint> &Indexes) const
	{
	Indexes.resize(m_ShapeCount);
	for (uint i = 0; i < m_ShapeCount; ++i)
		Indexes[i] = i;
	}

void ShapeSearcher::GetIncludes(vector<bool> &Includes) const
	{
	Includes.resize(m_ShapeCount);
	for (uint i = 0; i < m_ShapeCount; ++i)
		Includes[i] = true;
	}

#pragma once

#include "shapes.h"
#include <map>

class ShapeSearcher
	{
public:
	const PDBChain *m_Query = 0;
	uint m_PosA = UINT_MAX;
	uint m_PosB = UINT_MAX;
	uint m_PosC = UINT_MAX;
	const Shapes *m_Shapes = 0;
	uint m_ShapeCount = UINT_MAX;
	double m_Sigmas = 2.5;
	uint m_ShapeIndexA = UINT_MAX;
	uint m_ShapeIndexB = UINT_MAX;
	uint m_ShapeIndexC = UINT_MAX;
	double m_MinScoreABC = 0.55;
	double m_MinScoreShapePalm = 0.5;
	double m_MinSelfScore = 0.7;
	uint m_MaxTopHitCount = 8;
	uint m_MaxTopHitCountABC = 8;

	double m_ScoreABC = 0;
	double m_PalmScore = 0;

// m_IncludeShapes[i] is true/false to include i'th shape.
	vector<bool> m_IncludeShapes;

	vector<vector<uint> > m_SelfTopHits;

// Search results
	double m_Score = 0;
	vector<uint> m_ShapePosVec;
	vector<double> m_ShapeScores;

public:
	void ClearSearch()
		{
		m_Query = 0;
		m_PosA = UINT_MAX;
		m_PosB = UINT_MAX;
		m_PosC = UINT_MAX;
		m_Score = 0;
		m_ScoreABC = 0;
		m_PalmScore = 0;
		m_SelfTopHits.clear();
		m_IncludeShapes.clear();
		m_ShapePosVec.clear();
		m_ShapeScores.clear();
		}

	void Init(const Shapes &S);

	void ToTsv(FILE *f) const;

	uint GetQL() const
		{
		return m_Query->GetSeqLength();
		}

	const char *GetShapeName(uint ShapeIndex) const
		{
		asserta(m_Shapes != 0);
		asserta(ShapeIndex < SIZE(m_Shapes->m_Names));
		return m_Shapes->m_Names[ShapeIndex].c_str();
		}

	uint GetShapeLength(uint ShapeIndex) const
		{
		return m_Shapes->m_Lengths[ShapeIndex];
		}

	double GetMeanDist3(uint ShapeIndex1, uint ShapeIndex2,
	  uint Offset1, uint Offset2) const
		{
		return m_Shapes->GetMeanDist3(ShapeIndex1, ShapeIndex2, Offset1, Offset2);
		}

	double GetStdDev3(uint ShapeIndex1, uint ShapeIndex2,
	  uint Offset1, uint Offset2) const
		{
		return m_Shapes->GetStdDev3(ShapeIndex1, ShapeIndex2, Offset1, Offset2);
		}

	uint GetShapeCount() const { return m_ShapeCount; }
	void GetShapeIndexes(vector<uint> &Indexes) const;
	void GetIncludes(vector<bool> &Includes) const;
	void GetFirstIndexes(vector<uint> &Indexes) const;
	bool GetNextIndexes(vector<uint> &Indexes) const;

	void SetQuery(const PDBChain &Query);

	void GetDistRange(uint ShapeIndex, uint ShapeIndex2, 
	  uint &MinDist, uint &MaxDist) const;

	void SearchPalm(const PDBChain &Q);

	void SearchShape(uint ShapeIndex, const vector<uint> &PosVec,
	  double MinScore, uint Lo, uint Hi, char Letter, uint LetterOffset,
	  vector<uint> &HitPosVec, vector<double> &HitScores) const;

	void SearchShapeTopHit(uint ShapeIndex, const vector<uint> &PosVec,
	  double MinScore, uint Lo, uint Hi, char Letter, uint LetterOffset,
	  uint &Pos, double &Score) const;

	void SearchShapeSelf(uint ShapeIndex, double MinScore,
	  uint Lo, uint Hi, char Letter, uint LetterOffset,
	  vector<uint> &HitPosVec, vector<double> &HitScores) const;

	void SearchShapeSelfTop(uint ShapeIndex, double MinScore,
	  uint MaxHits, vector<uint> &HitPosVec) const;

	double GetSelfScore(uint ShapeIndex, uint Pos) const;

	double GetScore(uint ShapeIndex, uint Pos,
	  const vector<uint> &PosVec) const;

	double GetScoreShapePair(uint ShapeIndex1, uint ShapeIndex2,
	  uint Pos1, uint Pos2) const;

	double GetScoreShapes(const vector<uint> &ShapeIndexes,
	  const vector<uint> &PosVec) const;

	double GetScoreResiduePair(uint ShapeIndex1, uint ShapeIndex2,
	  uint Pos1, uint Pos2, uint Offset1, uint Offset2) const;

	double SearchABC(bool DoTrace = false);
	void Search(const vector<bool> &IncludeShapes);

	void TestABC1(const PDBChain &Chain,
	  const vector<string> &MotifSeqs,
	  double MinPredScore);

	void GetShapeSeq(uint ShapeIndex, string &Seq) const;
	void GetA(string &Seq) const;
	void GetB(string &Seq) const;
	void GetC(string &Seq) const;

	void GetSubSeq(uint Pos, uint n, string &Seq) const;

	void LogShape(const string &Msg, uint ShapeIndex, uint Pos) const;
	void GetColor(uint MotifIndex, string &Color) const;

	void ToPmlABC(FILE *f) const;
	void ToPml(FILE *f, const string &LoadName) const;

	void LogShapes(const vector<uint> &ShapeIndexes,
	  const vector<string> &ShapeSeqs) const;

	char GetGate() const;
	void GetFoundMotifsStr(string &Str) const;

public:
	static void TestABC(const Shapes &S,
	  const vector<PDBChain *> &Chains,
	  vector<vector<string> > &MotifSeqsVec,
	  double MinPredScore);

	static void SearchABCX(const Shapes &S,
	  const vector<PDBChain *> &Chains, double MinScoreX,
	  vector<string> &SeqXs, vector<double> &Scores);

	static bool GetBeforeABC(const vector<PDBChain *> &Chains,
	  vector<vector<string> > &ABCs, vector<string> &Xs);

	static bool JoinABCX1(
	  const vector<PDBChain *> &Chains,
	  const vector<string> &ChainLabels,
	  const map<string, vector<string> > &LabelToABC,
	  const map<string, string> &LabelToX,
	  vector<string> &Xs,
	  vector<vector<string> > &ABCs,
	  vector<vector<string> > &ABCXs,
	  vector<string> &Names,
	  vector<uint> &Lengths);

	static void JoinABCX2(
	  const vector<vector<string> > &ABCs,
	  const vector<string> &Xs,
	  bool BeforeABC,
	  vector<vector<string> > &ABCXs,
	  vector<string> &Names,
	  vector<uint> &Lengths);

	static void TrainABCX(const vector<PDBChain *> &Chains,
	  const vector<vector<string> > &ABCs,
	  const vector<string> &SeqsX1s,
	  bool BeforeABC,
	  Shapes &S1,
	  Shapes &S2,
	  vector<string> &SeqX2s,
	  vector<string> &SeqX3s,
	  vector<double> &ScoreX2s,
	  vector<double> &ScoreX3s);

	static void LogTrainStats(double MinScoreX,
	  const vector<PDBChain *> &Chains,
	  vector<string> &SeqX2s,
	  vector<string> &SeqX3s,
	  vector<double> &ScoreX2s,
	  vector<double> &ScoreX3s);
	};

double GetNormal(double Mu, double Sigma, double x);

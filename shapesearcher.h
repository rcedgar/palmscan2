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
	uint m_ShapeIndexA = UINT_MAX;
	uint m_ShapeIndexB = UINT_MAX;
	uint m_ShapeIndexC = UINT_MAX;
	uint m_MaxTopHitCount = 8;
	uint m_MaxTopHitCountABC = 8;
	vector<vector<uint> > m_SelfTopHits;
	double m_Sigmas = 2.5;

// m_*Shapes[i] is true/false to include i'th shape.
	vector<bool> m_SearchShapes;
	vector<bool> m_RequireShapes;
	bool m_SearchABCOnly = false;

// For IsHit()
	double m_MinSelfScoreABC = DBL_MAX;
	double m_MinSelfScoreNonABC = DBL_MAX;
	double m_MinABCScore = DBL_MAX;
	double m_MinPalmScore = DBL_MAX;

// For E-value
	double m_MeanFinalScore = 0.500537;
	double m_StdDevFinalScore = 0.0536741;
	double m_Log10DBSize = 6;

// Search results
	double m_ABCScore = 0;
	double m_DomScore = 0;
	bool m_Permuted = false;
	double m_FinalScore = 0;
	//double m_LEFPPM = DBL_MAX;

// Search results
	double m_Score = 0;
	vector<uint> m_ShapePosVec;
	vector<double> m_ShapeScores;

	string m_Class;

public:
	void ClearSearch()
		{
		m_Query = 0;
		m_PosA = UINT_MAX;
		m_PosB = UINT_MAX;
		m_PosC = UINT_MAX;
		m_Score = 0;
		m_ABCScore = 0;
		m_DomScore = 0;
		m_SelfTopHits.clear();
		m_ShapePosVec.clear();
		m_ShapeScores.clear();
		m_Permuted = false;
		m_Class = "";
		m_FinalScore = 0;
		//m_LEFPPM = DBL_MAX;
		}

	void Init(const Shapes &S);
	void LogParams() const;
	void SetParamOpts();
	void SetShapeIndexesABC();

	void ToTsv(FILE *f) const;

	bool IsHit() const;
	void SetClass();

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
	uint IncludesStrToBools(const string &What, const string &Str,
	  vector<bool> &Includes) const;
	void IncludesBoolsToStr(const vector<bool> &Includes,
	  string &Str) const;
	void BoolsToIndexVec1(const vector<bool> &Includes,
	  vector<uint> &ShapeIndexes) const;
	void BoolsToIndexVec2(const vector<bool> &Includes,
	  vector<uint> &ShapeIndexes) const;

	void SetQuery(const PDBChain &Query);

	void GetDistRange(uint ShapeIndex, uint ShapeIndex2, 
	  uint &MinDist, uint &MaxDist) const;

	void SearchDom(const PDBChain &Q);

	void SearchShape(uint ShapeIndex, const vector<uint> &PosVec,
	  double MinScore, uint Lo, uint Hi, char Letter, uint LetterOffset,
	  vector<uint> &HitPosVec, vector<double> &HitScores) const;

	void SearchShapeTopHit(uint ShapeIndex, const vector<uint> &PosVec,
	  double MinScore, uint Lo, uint Hi, char Letter, uint LetterOffset,
	  uint &Pos, double &Score) const;

	double SearchShapeSelf(uint ShapeIndex, double MinScore,
	  uint Lo, uint Hi, char Letter, uint LetterOffset,
	  vector<uint> &HitPosVec, vector<double> &HitScores) const;

	void SearchShapeSelfTop(uint ShapeIndex, double MinScore,
	  uint MaxHits, vector<uint> &HitPosVec) const;

	double GetSelfScore(uint ShapeIndex, uint Pos) const;

	double GetScore(uint ShapeIndex, uint Pos,
	  const vector<uint> &PosVec) const;

	double GetScoreShapePair(uint ShapeIndex1, uint ShapeIndex2,
	  uint Pos1, uint Pos2) const;

	double GetScoreShapes(const vector<uint> &PosVec) const;

	double GetScoreResiduePair(uint ShapeIndex1, uint ShapeIndex2,
	  uint Pos1, uint Pos2, uint Offset1, uint Offset2) const;

	void LogPairwiseScores(const vector<uint> &PosVec) const;

	void SearchABC(bool DoTrace = false);

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

	void TestABC1(const PDBChain &Chain, const vector<string> &MotifSeqs,
	  double MinPredScore);

	void CalibrateAdd(bool Hit) const;
	void SetFinalScore();
	//void SetLEFPPM();

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

	static void LogStats();
	static void StatsToFev(FILE *f);
	static void CalibrateWrite();
	};

double GetNormal(double Mu, double Sigma, double x);

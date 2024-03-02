#include "myutils.h"
#include "xprof.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "sort.h"
#include "xcluster.h"
#include "alpha.h"
#include <map>
#include <set>

void LogCountsMx(const vector<vector<uint> > &CountMx);
void LogScoreMx(const vector<vector<double> > &ScoreMx);
void GetLogOddsMx(const vector<double> &Freqs,
  const vector<vector<double> > &FreqMx,
  vector<vector<double> > &ScoreMx,
  double &H, double &RH);

static const double MINFRACT = 0.0;
static const uint ALPHASIZE = 20;
static const uint REVISITN = 0;
static const uint REVISITK = 32;
static const uint ITERS = 100;
static const double MINSCOREBASE = 1.0;
static const double MAXSCOREDELTA = 0.25;

/***
Input is created by xfeatures

palmscan2 \
  -xcluster d:/int/scop40/out/xfeatures.tsv
  -log xcluster.log

# head -n5 /d/int/dave_grant/scop40/aligned_xfeatvecs_0.6_0.8.tsv | columns.py
Q  Q_Ang_m2_p2  Q_Ang_m3_p3  Q_ED_p4  Q_ED_m4  Q_NU  Q_ND  R  R_Ang_m2_p2  R_Ang_m3_p3  R_ED_p4  R_ED_m4  R_NU  R_ND    QDom  QPos    RDom  RPos
N         83.1         62.8     9.57     11.7     0    14  L         13.9           22     11.6        0     1    19  d1a04a1    4  d1fsea_    3
L         72.9         48.9      5.6     8.06    19     4  L         18.5         10.6     5.74     12.9    21     5  d1a04a1    6  d1fsea_    4
T         34.4         74.9     6.27     9.21     2    16  T         34.9         83.4     6.27       13     3    13  d1a04a1    7  d1fsea_    5
P          125          107     6.25     9.57     0    15  K          125          109      6.2     13.1     1    13  d1a04a1    8  d1fsea_    6

***/

void XCluster::ReadFeatureTsv(const string &FileName)
	{
	m_Aminos.clear();
	m_Aminos1.clear();
	m_Aminos2.clear();
	m_FeatureValuesVec.clear();
	m_FeatureValuesVec1.clear();
	m_FeatureValuesVec2.clear();

	FILE *f = OpenStdioFile(FileName);
	string HdrLine;
	bool Ok = ReadLineStdioFile(f, HdrLine);
	vector<string> HdrFields;
	Split(HdrLine, HdrFields, '\t');
	asserta(SIZE(HdrFields) == 2*XFEATS + 6);
	asserta(HdrFields[0] == "Q");
	asserta(HdrFields[XFEATS+1] == "R");
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		{
		asserta(HdrFields[FeatureIndex+1] == 
		  "Q_" +  (string) XProf::GetFeatureName(FeatureIndex));
		asserta(HdrFields[XFEATS+2+FeatureIndex] == 
		  "R_" + (string) XProf::GetFeatureName(FeatureIndex));
		}

	set<pair<string, uint> > DoneSet;

	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2*XFEATS + 6);

		{
		asserta(SIZE(Fields[0]) == 1);
		char Amino = Fields[0][0];
		vector<double> Values;
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			double Value = StrToFloat(Fields[FeatureIndex+1]);
			Values.push_back(Value);
			}

		m_Aminos1.push_back(Amino);
		m_FeatureValuesVec1.push_back(Values);

		string Dom = Fields[2*XFEATS+2];
		uint Coord = StrToUint(Fields[2*XFEATS+3]);
		pair<string, uint> Pair(Dom, Coord);
		if (DoneSet.find(Pair) == DoneSet.end())
			{
			m_Aminos.push_back(Amino);
			m_FeatureValuesVec.push_back(Values);
			m_Doms.push_back(Dom);
			m_Coords.push_back(Coord);
			DoneSet.insert(Pair);
			}
		}

		{
		asserta(SIZE(Fields[XFEATS+1]) == 1);
		char Amino = Fields[XFEATS+1][0];
		vector<double> Values;
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			double Value = StrToFloat(Fields[XFEATS+2+FeatureIndex]);
			Values.push_back(Value);
			}

		m_Aminos2.push_back(Amino);
		m_FeatureValuesVec2.push_back(Values);

		string Dom = Fields[2*XFEATS+4];
		uint Coord = StrToUint(Fields[2*XFEATS+5]);
		pair<string, uint> Pair(Dom, Coord);
		if (DoneSet.find(Pair) == DoneSet.end())
			{
			m_Aminos.push_back(Amino);
			m_FeatureValuesVec.push_back(Values);
			m_Doms.push_back(Dom);
			m_Coords.push_back(Coord);
			DoneSet.insert(Pair);
			}
		}
		}
	const uint N = SIZE(m_FeatureValuesVec);
	const uint M = SIZE(m_FeatureValuesVec1);
	asserta(SIZE(m_Aminos) == N);
	asserta(SIZE(m_Doms) == N);
	asserta(SIZE(m_Coords) == N);

	asserta(SIZE(m_Aminos1) == M);
	asserta(SIZE(m_Aminos2) == M);
	asserta(SIZE(m_FeatureValuesVec2) == M);

	ProgressLog("%s unique feature vecs\n", IntToStr(N));
	ProgressLog("%s total feature vecs\n", IntToStr(M));
	CloseStdioFile(f);
	}

void XCluster::GetFreqs(const vector<uint> &RefIdxs)
	{
	asserta(SIZE(RefIdxs) == 20);

	m_CountVec.clear();
	m_CountMx.clear();
	m_Freqs.clear();
	m_FreqMx.clear();

	const uint PosPairCount = SIZE(m_Aminos1);
	assert(SIZE(m_Aminos2) == PosPairCount);
	assert(SIZE(m_FeatureValuesVec1) == PosPairCount);
	assert(SIZE(m_FeatureValuesVec2) == PosPairCount);

	uint LetterPairCount = 0;
	m_CountVec.resize(20, 0);
	m_CountMx.resize(20);
	for (uint i = 0; i < 20; ++i)
		m_CountMx[i].resize(20, 1);

	uint Counter = 0;
#pragma omp parallel for
	for (int iPairIndex = 0; iPairIndex < (int) PosPairCount; ++iPairIndex)
		{
		uint PairIndex = uint(iPairIndex);
		char AminoQ = m_Aminos1[PairIndex];
		char AminoR = m_Aminos2[PairIndex];

		const vector<double> &ValuesQ = m_FeatureValuesVec1[PairIndex];
		const vector<double> &ValuesR = m_FeatureValuesVec2[PairIndex];

		uint iq = GetMy3DLetter(AminoQ, ValuesQ, RefIdxs);
		uint ir = GetMy3DLetter(AminoR, ValuesR, RefIdxs);
		if (iq >= 20 || ir >= 20)
			continue;
#pragma omp critical
		{
		ProgressStep(Counter++, PosPairCount, "freqs");
		LetterPairCount += 2;
		m_CountVec[iq] += 1;
		m_CountVec[ir] += 1;
		m_CountMx[iq][ir] += 1;
		m_CountMx[ir][iq] += 1;
		}
		}

	double SumFreq = 0;
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		double Freq = double(m_CountVec[i])/(LetterPairCount);
		SumFreq += Freq;
		m_Freqs.push_back(Freq);
		//Log("Freqs[%c] = %8.6f\n", c, Freq);
		}
	//LogCountsMx(m_CountMx);

	assert(feq(SumFreq, 1.0));

	double SumFreq2 = 0;
	m_FreqMx.resize(20);
	for (uint i = 0; i < 20; ++i)
		{
		for (uint j = 0; j < 20; ++j)
			{
			uint n = m_CountMx[i][j];
			double Freq2 = double(n)/double(LetterPairCount);
			m_FreqMx[i].push_back(Freq2);
			SumFreq2 += Freq2;
			}
		}
	assert(feq(SumFreq2, 1.0));
	}

void XCluster::VecToTsv(FILE *f, uint Idx) const
	{
	if (f == 0)
		return;
	asserta(Idx < SIZE(m_Aminos));
	char aa = m_Aminos[Idx];
	asserta(Idx < SIZE(m_FeatureValuesVec));
	const vector<double> &v = m_FeatureValuesVec[Idx];
	const string &Dom = m_Doms[Idx];
	uint Coord = m_Coords[Idx];
	asserta(SIZE(v) == XFEATS);
	fprintf(f, "\t%u\t%c", Idx, aa);
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		fprintf(f, "\t%.3g", v[FeatureIndex]);
	fprintf(f, "\t%s\t%u", Dom.c_str(), Coord);
	}

void XCluster::HitToTsv(double Score, uint Idx, uint CentroidIdx) const
	{
	if (g_ftsv == 0)
		return;
	fprintf(g_ftsv, "H\t%.1f", Score);
	VecToTsv(g_ftsv, Idx);
	VecToTsv(g_ftsv, CentroidIdx);
	fprintf(g_ftsv, "\n");
	}

void XCluster::CentroidToTsv(uint CentroidIdx) const
	{
	if (g_ftsv == 0)
		return;
	fprintf(g_ftsv, "S\t*");
	VecToTsv(g_ftsv, CentroidIdx);
	fprintf(g_ftsv, "\n");
	}

double XCluster::GetScore(uint Idx1, uint Idx2) const
	{
	char aa1 = m_Aminos[Idx1];
	char aa2 = m_Aminos[Idx2];
	const vector<double> &Values1 = m_FeatureValuesVec[Idx1];
	const vector<double> &Values2 = m_FeatureValuesVec[Idx2];
	double Score = XProf::GetScore2(aa1, aa2, Values1, Values2);
	return Score;
	}

void XCluster::TopsToTsv(FILE *f, const vector<uint> &RefIdxs) const
	{
	if (f == 0)
		return;
	asserta(SIZE(RefIdxs) == 20);
	const uint N = SIZE(m_FeatureValuesVec);
	for (uint Idx = 0; Idx < N; ++Idx)
		{
		ProgressStep(Idx, N, "TopsToTsv");
		double Score;
		uint Letter;
		uint TopIdx = GetTopIdx(Idx, RefIdxs, &Score, &Letter);
		asserta(Letter < 20);
		fprintf(f, "%c\t%.3g", g_LetterToCharAmino[Letter], Score);
		VecToTsv(f, Idx);
		fprintf(f, "\n");
		}
	}

uint XCluster::GetMy3DLetter(char Amino, const vector<double> &Values,
  const vector<uint> &RefIdxs) const
	{
	asserta(SIZE(RefIdxs) == 20);
	double BestScore = -999;
	uint BestIndex = UINT_MAX;
	const uint N = SIZE(RefIdxs);
	asserta(N > 0);
	for (uint i = 0; i < N; ++i)
		{
		uint RefIdx = RefIdxs[i];
		char RefAmino = m_Aminos[RefIdx];
		const vector<double> &RefValues = m_FeatureValuesVec[RefIdx];
		double Score =
		  XProf::GetScore2(Amino, RefAmino, Values, RefValues);
		if (Score >= BestScore)
			{
			BestScore = Score;
			BestIndex = i;
			}
		}
	asserta(BestIndex != UINT_MAX);
	return BestIndex;
	}

uint XCluster::GetTopIdx(uint QueryIdx, const vector<uint> &RefIdxs,
  double *ptrScore, uint *ptrIndex) const
	{
	double BestScore = -999;
	uint BestRefIdx = UINT_MAX;
	uint BestRefIndex = UINT_MAX;
	const uint N = SIZE(RefIdxs);
	asserta(N > 0);
	for (uint i = 0; i < N; ++i)
		{
		uint RefIdx = RefIdxs[i];
		double Score = GetScore(QueryIdx, RefIdx);
		if (Score >= BestScore)
			{
			BestScore = Score;
			BestRefIdx = RefIdx;
			BestRefIndex = i;
			}
		}
	if (ptrScore != 0)
		*ptrScore = BestScore;
	if (ptrIndex != 0)
		*ptrIndex = BestRefIndex;
	asserta(BestRefIdx != UINT_MAX);
	return BestRefIdx;
	}

void XCluster::Cluster(
	const vector<uint> &InputIdxs,
	double MinScore,
	vector<uint> &CentroidIdxs,
	map<uint, uint> &IdxToCentroidIdx,
	map<uint, vector<uint> > &CentroidIdxToMemberIdxs)
	{
	CentroidIdxs.clear();
	IdxToCentroidIdx.clear();
	CentroidIdxToMemberIdxs.clear();
	map<uint, uint> CentroidIdxToCentroidIndex;

	const uint N = SIZE(m_Aminos);
	asserta(SIZE(m_FeatureValuesVec) == N);
	vector<uint> Bins(XFEATS);
	for (uint kkk = 0; kkk < N; ++kkk)
		{
		uint Idx = InputIdxs[kkk];
		ProgressStep(kkk, N, "Clustering minscore %.4f, %u centroids",
			MinScore, SIZE(CentroidIdxs));
		bool Hit = false;
		vector<uint> Empty;
		CentroidIdxToMemberIdxs[Idx] = Empty;
		double BestScore = -999;
		uint BestCentroidIdx = UINT_MAX;
		for (uint j = 0; j < SIZE(CentroidIdxs); ++j)
			{
			uint CentroidIdx = CentroidIdxs[j];
			double Score = GetScore(Idx, CentroidIdx);
			if (Score >= MinScore)
				{
				Hit = true;
				if (Score >= BestScore)
					{
					BestScore = Score;
					BestCentroidIdx = CentroidIdx;
					}
				if (!opt_besthit)
					break;
				}
			}

		if (Hit)
			{
			IdxToCentroidIdx[Idx] = BestCentroidIdx;
			CentroidIdxToMemberIdxs[BestCentroidIdx].push_back(Idx);
			//HitToTsv(BestScore, Idx, BestCentroidIdx);
			}
		else
			{
			uint CentroidIndex = SIZE(CentroidIdxs);
			CentroidIdxToCentroidIndex[Idx] = CentroidIndex;
			CentroidIdxs.push_back(Idx);
			IdxToCentroidIdx[Idx] = Idx;
			CentroidIdxToMemberIdxs[Idx].push_back(Idx);
			//CentroidToTsv(Idx);
			}
		}
	asserta(SIZE(IdxToCentroidIdx) == N);
	asserta(SIZE(CentroidIdxToMemberIdxs) == N);
	}

void XCluster::ScoreMxToDistMx(
	const vector<vector<double> > &ScoreMx,
	vector<vector<double> > &DistMx)
	{
	DistMx.clear();
	const uint N = SIZE(ScoreMx);
	vector<double> Scores;
	for (uint i = 0; i < N; ++i)
		{
		asserta(SIZE(ScoreMx[i]) == N);
		for (uint j = 0; j < i; ++j)
			Scores.push_back(ScoreMx[i][j]);
		}

	const uint M = SIZE(Scores);
	sort(Scores.begin(), Scores.end());
	double MinScore = Scores[0];
	double LoScore = Scores[M/8];
	double MedScore = Scores[M/2];
	double HiScore = Scores[7*M/8];
	double MaxScore = Scores[M-1];
	ProgressLog("Scores min %.3g, lo %.3g, med %.3g, hi %.3g, max %.3g\n",
	  MinScore, LoScore, MedScore, HiScore, MaxScore);

	DistMx.resize(N);
	for (uint i = 0; i < N; ++i)
		DistMx[i].resize(N);

	for (uint i = 0; i < N; ++i)
		{
		DistMx[i][i] = 0;
		for (uint j = 0; j < i; ++j)
			{
			double Score = ScoreMx[i][j];
			asserta(feq(ScoreMx[j][i], Score));
			if (Score < LoScore)
				Score = LoScore;
			if (Score > HiScore)
				Score = HiScore;
			double Dist = (HiScore - Score)/(HiScore - LoScore);
			asserta(Dist >= 0 && Dist <= 1);
			DistMx[i][j] = Dist;
			DistMx[j][i] = Dist;
			}
		}
	}

void XCluster::GetScoreMx(const vector<uint> &Idxs,
  vector<vector<double> > &DistMx) const
	{
	const uint N = SIZE(Idxs);
	DistMx.clear();
	if (N == 0)
		return;
	DistMx.resize(N);
	for (uint i = 0; i < N; ++i)
		DistMx[i].resize(N, DBL_MAX);

	for (uint i = 0; i < N; ++i)
		{
		uint Idxi = Idxs[i];
		for (uint j = 0; j <= i; ++j)
			{
			uint Idxj = Idxs[j];
			double Score = GetScore(Idxi, Idxj);
			Log("%5d  %5d  %.3g\n", i, j, Score);
			DistMx[i][j] = Score;
			DistMx[j][i] = Score;
			}
		}
	}

void XCluster::GetClusterSizes(
	const map<uint, vector<uint> > &CentroidIdxToMemberIdxs,
	const vector<uint> &CentroidIdxs,
	vector<uint> &Order,
	vector<uint> &Sizes) const
	{
	Sizes.clear();
	const uint ClusterCount = SIZE(CentroidIdxs);
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		uint CentroidIdx = CentroidIdxs[ClusterIndex];
		asserta(CentroidIdx < SIZE(CentroidIdxToMemberIdxs));
		map<uint, vector<uint> >::const_iterator iter =
		  CentroidIdxToMemberIdxs.find(CentroidIdx);
		asserta(iter != CentroidIdxToMemberIdxs.end());
		//const vector<uint> &MemberIdxs = CentroidIdxToMemberIdxs[CentroidIdx];
		const vector<uint> &MemberIdxs = iter->second;
		Sizes.push_back(SIZE(MemberIdxs));
		}
	Order.resize(ClusterCount);
	QuickSortOrderDesc(Sizes.data(), ClusterCount, Order.data());
	}

void XCluster::ClustersToTsv(
	const vector<uint> &CentroidIdxs,
	const map<uint, vector<uint> > &CentroidIdxToMemberIdxs) const
	{
	const uint ClusterCount = SIZE(CentroidIdxs);

	vector<uint> Order;
	vector<uint> Sizes;
	GetClusterSizes(CentroidIdxToMemberIdxs, CentroidIdxs,
	  Order, Sizes);
	for (uint k = 0; k < ClusterCount; ++k)
		{
		uint ClusterIndex = Order[k];
		uint CentroidIdx = CentroidIdxs[ClusterIndex];
		uint Size = Sizes[ClusterIndex];
		if (g_ftsv != 0)
			{
			fprintf(g_ftsv, "C\t%u", Size);
			VecToTsv(g_ftsv, CentroidIdx);
			fprintf(g_ftsv, "\n");
			}
		}
	}

void XCluster::MxToTsv(FILE *f, 
  const vector<vector<double> > &Mx) const
	{
	if (f == 0)
		return;
	const uint N = SIZE(Mx);
	for (uint i = 0; i < N; ++i)
		{
		asserta(SIZE(Mx[i]) == N);
		for (uint j = 0; j < N; ++j)
			fprintf(f, "%u\t%u\t%.3g\n", i, j, Mx[i][j]);
		}
	}

uint XCluster::GetFirstRandomHit(uint Letter, const vector<uint> &RefIdxs) const
	{
	const uint N = SIZE(m_FeatureValuesVec);
	for (uint i = 0; i < 1024; ++i)
		{
		uint QueryIdx = randu32()%N;

		uint Index;
		double Score;
		GetTopIdx(QueryIdx, RefIdxs, &Score, &Index);
		if (Index == Letter)
			return QueryIdx;
		}
	return randu32()%N;
	}

void XCluster::LogDCs(const vector<uint> AlphaIdxs,
	double ExpScore) const
	{
	asserta(SIZE(AlphaIdxs) == 20);
	for (uint i = 0; i < 20; ++i)
		{
		uint CentroidIdx = AlphaIdxs[i];
		char aa = m_Aminos[CentroidIdx];
		const vector<double> &v = m_FeatureValuesVec[CentroidIdx];
		Log("DefineCentroid(%u, '%c'", i, aa);
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			Log(", %.8g", v[FeatureIndex]);
		Log("); // %.4f\n", ExpScore);
		}
	}

void cmd_xcluster()
	{
	const string &InputFileName = opt_xcluster;
	XCluster XC;
	XC.ReadFeatureTsv(InputFileName);
	const uint N = SIZE(XC.m_Aminos);
	asserta(!optset_minscore);
	//double OptMinScore = opt_minscore;

	XProf::InitScoreTable();

	double BestExpectedScore = 0;
	vector<double> ExpScores;
	vector<vector<uint> > SelectedIdxVec;
	for (uint Iter = 0; Iter < ITERS; ++Iter)
		{
		vector<uint> InputIdxs;
		for (uint i = 0; i < N; ++i)
			InputIdxs.push_back(i);
		random_shuffle(InputIdxs.begin(), InputIdxs.end());

		vector<uint> CentroidIdxs;
		map<uint, uint> IdxToCentroidIdx;
		map<uint, vector<uint> > CentroidIdxToMemberIdxs;
		double ScoreDeltaFract = (randu32()%100)/100.0;
		double ScoreDelta = ScoreDeltaFract*MAXSCOREDELTA;
		if (randu32()%2 == 0)
			ScoreDelta = -ScoreDelta;
		double MinScore = MINSCOREBASE + ScoreDelta;
		XC.Cluster(InputIdxs, MinScore, CentroidIdxs, IdxToCentroidIdx,
		  CentroidIdxToMemberIdxs);
		const uint ClusterCount = SIZE(CentroidIdxs);
		ProgressLog("%u clusters, minscore %.4f\n", ClusterCount, MinScore);
		if (ClusterCount < 20)
			{
			Warning("%u clusters < 20", ClusterCount);
			continue;
			}

		vector<uint> Order;
		vector<uint> Sizes;
		XC.GetClusterSizes(CentroidIdxToMemberIdxs, CentroidIdxs,
		  Order, Sizes);

		vector<uint> SelectedCentroidIdxs;
		double SumFract = 0;
		for (uint k = 0; k < ALPHASIZE; ++k)
			{
			uint ClusterIndex = Order[k];
			uint CentroidIdx = CentroidIdxs[ClusterIndex];
			uint Size = SIZE(CentroidIdxToMemberIdxs[CentroidIdx]);
			double Fract = double(Size)/N;
			SumFract += Fract;
			if (Fract < MINFRACT)
				continue;
			SelectedCentroidIdxs.push_back(CentroidIdx);
			}

		XC.GetFreqs(SelectedCentroidIdxs);
		double H, RH;
		vector<vector<double> > ScoreMx;
		GetLogOddsMx(XC.m_Freqs, XC.m_FreqMx, ScoreMx, H, RH);
		if (RH > BestExpectedScore)
			{
			XC.LogDCs(SelectedCentroidIdxs, RH);
			BestExpectedScore = RH;
			}

		//LogScoreMx(ScoreMx);
		ProgressLog("Iter %u/%u cltscore %.3g, expected score %.3f (max %.3f)\n",
		  Iter, ITERS, MinScore, RH, BestExpectedScore);

		ExpScores.push_back(RH);
		asserta(SIZE(SelectedCentroidIdxs) == 20);
		SelectedIdxVec.push_back(SelectedCentroidIdxs);
		}

	vector<uint> Order(ITERS);
	QuickSortOrderDesc(ExpScores.data(), ITERS, Order.data());

	for (uint RevisitIndex = 0; RevisitIndex < min(REVISITN, ITERS); ++RevisitIndex)
		{
		uint kkk = Order[RevisitIndex];
		asserta(kkk < SIZE(ExpScores));
		asserta(kkk < SIZE(SelectedIdxVec));
		double OrigExpectedScore = ExpScores[kkk];
		const vector<uint> SelectedCentroidIdxs = SelectedIdxVec[kkk];
		asserta(SIZE(SelectedCentroidIdxs) == 20);
		Progress("Revisit %u exp score %.4f\n", RevisitIndex, OrigExpectedScore);

		for (uint k = 0; k < REVISITK; ++k)
			{
			asserta(SIZE(SelectedCentroidIdxs) == 20);
			vector<uint> VariantIdxs;
			for (uint i = 0; i < 20; ++i)
				{
				uint VariantIdx = XC.GetFirstRandomHit(i, SelectedCentroidIdxs);
				VariantIdxs.push_back(VariantIdx);
				}

			XC.GetFreqs(VariantIdxs);
			double H, RH;
			vector<vector<double> > ScoreMx;
			GetLogOddsMx(XC.m_Freqs, XC.m_FreqMx, ScoreMx, H, RH);
			BestExpectedScore = max(RH, BestExpectedScore);

			//LogScoreMx(ScoreMx);
			ProgressLog("Variant %u/%u, expected score %.3f (orig %.3f, max %.3f)\n",
			  RevisitIndex, k, RH, OrigExpectedScore, BestExpectedScore);

			ExpScores.push_back(RH);
			SelectedIdxVec.push_back(VariantIdxs);
			asserta(SIZE(SelectedCentroidIdxs) == 20);
			}
		}

	uint top = 0;
	double topscore = 0;
	for (uint i = 0; i < SIZE(ExpScores); ++i)
		{
		if (ExpScores[i] > topscore)
			{
			topscore = ExpScores[i];
			top = i;
			}
		}

	const vector<uint> SelectedCentroidIdxs = SelectedIdxVec[top];

	XC.GetFreqs(SelectedCentroidIdxs);
	double H, RH;
	vector<vector<double> > ScoreMx;
	GetLogOddsMx(XC.m_Freqs, XC.m_FreqMx, ScoreMx, H, RH);
	BestExpectedScore = max(RH, BestExpectedScore);
	XC.LogDCs(SelectedCentroidIdxs, topscore);

	LogScoreMx(ScoreMx);
	ProgressLog("Final expected score %.3g (max %.3g)\n",
		RH, BestExpectedScore);
	}

#include "myutils.h"
#include "xprof.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "sort.h"
#include "xcluster.h"
#include "alpha.h"
#include <map>
#include <set>

static double MINFRACT = 0.0;
static uint TOPN = 20;

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
	m_FeatureValuesVec.clear();

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
		string Dom = Fields[2*XFEATS+2];
		uint Coord = StrToUint(Fields[2*XFEATS+3]);
		pair<string, uint> Pair(Dom, Coord);
		if (DoneSet.find(Pair) == DoneSet.end())
			{
			for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
				{
				double Value = StrToFloat(Fields[FeatureIndex+1]);
				Values.push_back(Value);
				}
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
		string Dom = Fields[2*XFEATS+4];
		uint Coord = StrToUint(Fields[2*XFEATS+5]);
		pair<string, uint> Pair(Dom, Coord);
		if (DoneSet.find(Pair) == DoneSet.end())
			{
			for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
				{
				double Value = StrToFloat(Fields[XFEATS+2+FeatureIndex]);
				Values.push_back(Value);
				}
			m_Aminos.push_back(Amino);
			m_FeatureValuesVec.push_back(Values);
			m_Doms.push_back(Dom);
			m_Coords.push_back(Coord);
			DoneSet.insert(Pair);
			}
		}
		}
	const uint N = SIZE(m_FeatureValuesVec);
	asserta(SIZE(m_Aminos) == N);
	asserta(SIZE(m_Doms) == N);
	asserta(SIZE(m_Coords) == N);
	ProgressLog("%s feature vecs\n", IntToStr(N));
	CloseStdioFile(f);
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
	vector<uint> &IdxToCentroidIdx,
	vector<vector<uint> > &CentroidIdxToMemberIdxs)
	{
	CentroidIdxs.clear();
	IdxToCentroidIdx.clear();
	CentroidIdxToMemberIdxs.clear();
	map<uint, uint> CentroidIdxToCentroidIndex;

	const uint N = SIZE(InputIdxs);
	vector<uint> Bins(XFEATS);
	for (uint i = 0; i < N; ++i)
		{
		uint Idx = InputIdxs[i];
		ProgressStep(i, N, "Clustering minscore %.1f, %u centroids",
			MinScore, SIZE(CentroidIdxs));
		bool Hit = false;
		vector<uint> Empty;
		CentroidIdxToMemberIdxs.push_back(Empty);
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
			IdxToCentroidIdx.push_back(BestCentroidIdx);
			CentroidIdxToMemberIdxs[BestCentroidIdx].push_back(Idx);
			HitToTsv(BestScore, Idx, BestCentroidIdx);
			}
		else
			{
			uint CentroidIndex = SIZE(CentroidIdxs);
			CentroidIdxToCentroidIndex[Idx] = CentroidIndex;
			CentroidIdxs.push_back(Idx);
			IdxToCentroidIdx.push_back(Idx);
			CentroidIdxToMemberIdxs[Idx].push_back(Idx);
			CentroidToTsv(Idx);
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
	const vector<vector<uint> > &CentroidIdxToMemberIdxs,
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
		const vector<uint> &MemberIdxs = CentroidIdxToMemberIdxs[CentroidIdx];
		Sizes.push_back(SIZE(MemberIdxs));
		}
	Order.resize(ClusterCount);
	QuickSortOrderDesc(Sizes.data(), ClusterCount, Order.data());
	}

void XCluster::ClustersToTsv(
	const vector<uint> &CentroidIdxs,
	const vector<vector<uint> > &CentroidIdxToMemberIdxs) const
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

void cmd_xcluster()
	{
	const string &InputFileName = opt_xcluster;
	XCluster XC;
	XC.ReadFeatureTsv(InputFileName);

	asserta(optset_minscore);
	double MinScore = opt_minscore;

	XProf::InitScoreTable();

	vector<uint> CentroidIdxs;
	vector<uint> IdxToCentroidIdx;
	vector<vector<uint> > CentroidIdxToMemberIdxs;
	vector<uint> InputIdxs;
	const uint N = SIZE(XC.m_Aminos);
	for (uint Idx = 0; Idx < N; ++Idx)
		InputIdxs.push_back(Idx);

	XC.Cluster(InputIdxs, MinScore,
	  CentroidIdxs, IdxToCentroidIdx, CentroidIdxToMemberIdxs);
	const uint ClusterCount = SIZE(CentroidIdxs);
	ProgressLog("%u clusters\n", ClusterCount);

	XC.ClustersToTsv(CentroidIdxs, CentroidIdxToMemberIdxs);

	vector<uint> Order;
	vector<uint> Sizes;
	XC.GetClusterSizes(CentroidIdxToMemberIdxs, CentroidIdxs,
	  Order, Sizes);

	vector<uint> SelectedCentroidIdxs;
	double SumFract = 0;
	for (uint k = 0; k < min(ClusterCount, TOPN); ++k)
		{
		uint ClusterIndex = Order[k];
		uint CentroidIdx = CentroidIdxs[ClusterIndex];
		uint Size = SIZE(CentroidIdxToMemberIdxs[CentroidIdx]);
		double Fract = double(Size)/N;
		SumFract += Fract;
		if (Fract < MINFRACT)
			continue;
		Log("[%3u]  %5u  %7.5f\n", k, Size, Fract);
		SelectedCentroidIdxs.push_back(CentroidIdx);
		}
	Log("Total fract %.5f\n", SumFract);

	const uint SelectedClusterCount = SIZE(SelectedCentroidIdxs);
	ProgressLog("%u / %u clusters selected\n",
		SelectedClusterCount, ClusterCount);

	vector<vector<double> > ScoreMx;
	XC.GetScoreMx(SelectedCentroidIdxs, ScoreMx);

	vector<vector<double> > DistMx;
	XC.ScoreMxToDistMx(ScoreMx, DistMx);
	asserta(SIZE(DistMx) == SelectedClusterCount);
	asserta(optset_output);
	FILE *f = CreateStdioFile(opt_output);
	fprintf(f, "%u\n", SelectedClusterCount);
	for (uint i = 0; i < SelectedClusterCount; ++i)
		{
		uint CentroidIdx = SelectedCentroidIdxs[i];
		uint Size = SIZE(CentroidIdxToMemberIdxs[CentroidIdx]);
		asserta(Size >= 1);
		double Fract = double(Size)/N;
		fprintf(f, "%u\t%u\t%.5f", i, Size, Fract);
		XC.VecToTsv(f, CentroidIdx);
		fprintf(f, "\n");
		}
	XC.MxToTsv(f, DistMx);
	CloseStdioFile(f);

	for (uint i = 0; i < SelectedClusterCount; ++i)
		{
		uint CentroidIdx = SelectedCentroidIdxs[i];
		uint Size = SIZE(CentroidIdxToMemberIdxs[CentroidIdx]);
		char aa = XC.m_Aminos[CentroidIdx];
		const vector<double> &v = XC.m_FeatureValuesVec[CentroidIdx];
		Log("DefineCentroid(%u, '%c'", i, aa);
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			Log(", %.3g", v[FeatureIndex]);
		Log(");\n");
		}

	FILE *f2 = CreateStdioFile(opt_output2);
	XC.TopsToTsv(f2, SelectedCentroidIdxs);
	CloseStdioFile(f2);
	}

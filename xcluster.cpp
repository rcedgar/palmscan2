#include "myutils.h"
#include "xprof.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "sort.h"
#include "xcluster.h"
#include "xbinner.h"
#include "alpha.h"
#include <map>
#include <set>

void LogCountsMx(const vector<vector<uint> > &CountMx);
void LogScoreMx(const vector<vector<double> > &ScoreMx);

static const double MINFRACT = 0.0;
static const uint ALPHASIZE = 20;
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
	m_X2.FromTsv(FileName);
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

	const uint N = SIZE(m_X2.m_Aminos);
	asserta(SIZE(m_X2.m_FeatureValuesVec) == N);
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

double XCluster::GetScore(uint Idx1, uint Idx2) const
	{
	char AminoChar1 = m_X2.m_Aminos[Idx1];
	char AminoChar2 = m_X2.m_Aminos[Idx2];
	uint AminoLetter1 = g_CharToLetterAmino[AminoChar1];
	uint AminoLetter2 = g_CharToLetterAmino[AminoChar2];
	const vector<double> &Values1 = m_X2.m_FeatureValuesVec[Idx1];
	const vector<double> &Values2 = m_X2.m_FeatureValuesVec[Idx2];
	double Score = 
	  XProf::GetScore_Letters2(AminoLetter1, AminoLetter2, Values1, Values2);
	return Score;
	}

void XCluster::DefineBinnerLetters(const vector<uint> AlphaIdxs) const
	{
	asserta(SIZE(AlphaIdxs) == 20);
	for (uint i = 0; i < 20; ++i)
		{
		uint CentroidIdx = AlphaIdxs[i];
		char aa = m_X2.m_Aminos[CentroidIdx];
		const vector<double> &Values = m_X2.m_FeatureValuesVec[CentroidIdx];
		XBinner::SetCentroid(i, aa, Values);
		}
	}

void XCluster::LogDCs(const vector<uint> AlphaIdxs,
	double ExpScore) const
	{
	asserta(SIZE(AlphaIdxs) == 20);
	for (uint i = 0; i < 20; ++i)
		{
		uint CentroidIdx = AlphaIdxs[i];
		char aa = m_X2.m_Aminos[CentroidIdx];
		const vector<double> &v = m_X2.m_FeatureValuesVec[CentroidIdx];
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
	const X2Data &X2 = XC.m_X2;

	const uint N = SIZE(XC.m_X2.m_Aminos);
	asserta(!optset_minscore);

	XProf::InitScoreTable();
	XBinnerC XBC;
	XBC.m_ptrX2 = &X2;

	double BestExpScore = 0;
	vector<double> ExpScores;
	vector<vector<uint> > SelectedIdxVec;
	for (uint Iter = 0; Iter < ITERS; ++Iter)
		{
		vector<uint> InputIdxs;
		for (uint i = 0; i < N; ++i)
			InputIdxs.push_back(i);
		Shuffle(InputIdxs);

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

		XBC.m_Idxs = SelectedCentroidIdxs;
		vector<double> Freqs;
		vector<vector<double> > FreqMx;
		XBC.GetFreqs(X2, Freqs, FreqMx);

		XBinner XB;

		vector<vector<double> > ScoreMx;
		double ExpScore =
		  XBinner::GetLogOddsMx(Freqs, FreqMx, ScoreMx);
		if (ExpScore > BestExpScore)
			{
			BestExpScore = ExpScore;
			XC.LogDCs(SelectedCentroidIdxs, ExpScore);
			}

		ProgressLog("Iter %u/%u radius %.3g, expected score %.3f (max %.3f)\n",
		  Iter, ITERS, MinScore, ExpScore, BestExpScore);

		ExpScores.push_back(ExpScore);
		asserta(SIZE(SelectedCentroidIdxs) == 20);
		SelectedIdxVec.push_back(SelectedCentroidIdxs);
		}
	}

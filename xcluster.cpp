#include "myutils.h"
#include "xprof.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "sort.h"
#include <set>

static const uint SUBSET1 = 1024;
static const uint SUBSET2 = 128;
static double MINSCORE1 = 2.5;

/***
Input is created by xfeatures

palmscan2 \
  -xcluster d:/int/scop40/out/xfeatures.tsv
  -log xcluster.log

# head d:/int/scop40/out/xfeatures.tsv | columns.py
aa  Ang_m2_p2  Ang_m3_p3  ED_p4  ED_m4  NU  ND
 A          0          0   6.07      0   0   0
 Y          0          0   5.92      0   7   7
 I        113          0   6.38      0   0  11
 A        115       34.1   6.33      0   2  12
 K        109         28   6.23   6.07  15   9

***/

void ReadXFeaturesTsv(const string &FileName,
  vector<char> &Aminos, vector<vector<double> > &FeatureValuesVec)
	{
	Aminos.clear();
	FeatureValuesVec.clear();
	FILE *f = OpenStdioFile(FileName);
	string HdrLine;
	bool Ok = ReadLineStdioFile(f, HdrLine);
	vector<string> HdrFields;
	Split(HdrLine, HdrFields, '\t');
	asserta(SIZE(HdrFields) == XFEATS + 1);
	asserta(HdrFields[0] == "aa");
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		asserta(HdrFields[FeatureIndex+1] == 
		  (string) XProf::GetFeatureName(FeatureIndex));

	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == XFEATS+1);
		asserta(SIZE(Fields[0]) == 1);
		char Amino = Fields[0][0];
		vector<double> Values;
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			double Value = StrToFloat(Fields[FeatureIndex+1]);
			Values.push_back(Value);
			}
		Aminos.push_back(Amino);
		FeatureValuesVec.push_back(Values);
		}
	}

uint SelectFinal(const set<uint> &Final20,
  vector<uint> &SubsetIdxs2,
  const vector<vector<double> > &ScoreMx2)
	{
	asserta(SIZE(SubsetIdxs2) == SUBSET2);
	asserta(SIZE(ScoreMx2) == SUBSET2);
	random_shuffle(SubsetIdxs2.begin(), SubsetIdxs2.end());
	if (Final20.empty())
		return SubsetIdxs2[0];

	uint BestIdx = UINT_MAX;
	double LowestSumScore = DBL_MAX;
	for (uint i = 0; i < SUBSET2; ++i)
		{
		uint Idxi = SubsetIdxs2[i];
		if (Final20.find(Idxi) != Final20.end())
			continue;
		double SumScore = 0;
		for (uint j = 0; j < SUBSET2; ++j)
			{
			uint Idxj = SubsetIdxs2[j];
			if (Final20.find(Idxj) == Final20.end())
				continue;
			double Score = ScoreMx2[i][j];
			SumScore += Score;
			}
		if (SumScore < LowestSumScore)
			{
			LowestSumScore = SumScore;
			BestIdx = Idxi;
			}
		}
	asserta(BestIdx != UINT_MAX);
	return BestIdx;
	}

void cmd_xcluster()
	{
	const string &InputFileName = opt_xcluster;
	vector<char> Aminos;
	vector<vector<double> > FeatureValuesVec;
	ReadXFeaturesTsv(InputFileName, Aminos, FeatureValuesVec);
	const uint N = SIZE(FeatureValuesVec);
	asserta(SIZE(Aminos) == N);
	ProgressLog("%s feature vecs\n", IntToStr(N));

	XProf::InitScoreTable();

	vector<uint> SubsetIdxs1;
	SubsetIdxs1.reserve(N);
	for (uint i = 0; i < N; ++i)
		SubsetIdxs1.push_back(i);
	random_shuffle(SubsetIdxs1.begin(), SubsetIdxs1.end());
	SubsetIdxs1.resize(SUBSET1);

	vector<uint> Counts(SUBSET1);
	vector<uint> Bins(XFEATS);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Pass 1");
		char aai = Aminos[i];
		const vector<double> Valuesi = FeatureValuesVec[i];
		for (uint j = 0; j < SUBSET1; ++j)
			{
			uint SubsetIdx = SubsetIdxs1[j];
			char aaj = Aminos[SubsetIdx];
			const vector<double> Valuesj = FeatureValuesVec[SubsetIdx];
			for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
				{
				double Valuei = Valuesi[FeatureIndex];
				double Valuej = Valuesj[FeatureIndex];
				double Diff = XProf::GetDiff(FeatureIndex, Valuei, Valuej);
				uint Bin = XProf::GetFeatureBin(FeatureIndex, Diff);
				Bins[FeatureIndex] = Bin;
				}
			double Score = XProf::GetScore(aai, aaj, Bins);
			if (Score >= MINSCORE1)
				Counts[j] += 1;
			}
		}

	vector<uint> Order(SUBSET1);
	QuickSortOrderDesc(Counts.data(), SUBSET1, Order.data());

	vector<uint> SubsetIdxs2;
	vector<uint> SubsetSizes2;
	for (uint i = 0; i < SUBSET2; ++i)
		{
		uint k = Order[i];
		uint SubsetIdx = SubsetIdxs1[k];
		uint Count = Counts[k];
		SubsetIdxs2.push_back(SubsetIdx);
		SubsetSizes2.push_back(Count);
		Log("[%5u]  %7u  %7.1f%%\n", i, Count, GetPct(Count, N));
		}

	vector<vector<double> > ScoreMx2(SUBSET2);
	for (uint i = 0; i < SUBSET2; ++i)
		ScoreMx2[i].resize(SUBSET2, DBL_MAX);

	for (uint i = 0; i < SUBSET2; ++i)
		{
		uint Idxi = SubsetIdxs2[i];
		char aai = Aminos[Idxi];
		const vector<double> Valuesi = FeatureValuesVec[Idxi];
		for (uint j = 0; j < i; ++j)
			{
			uint Idxj = SubsetIdxs2[j];
			char aaj = Aminos[Idxj];
			const vector<double> Valuesj = FeatureValuesVec[Idxj];
			for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
				{
				double Valuei = Valuesi[FeatureIndex];
				double Valuej = Valuesj[FeatureIndex];
				double Diff = XProf::GetDiff(FeatureIndex, Valuei, Valuej);
				uint Bin = XProf::GetFeatureBin(FeatureIndex, Diff);
				Bins[FeatureIndex] = Bin;
				}
			double Score = XProf::GetScore(aai, aaj, Bins);
			ScoreMx2[i][j] = Score;
			ScoreMx2[j][i] = Score;
			}
		}

	set<uint> Final20;
	for (uint i = 0; i < 20; ++i)
		{
		uint Idx = SelectFinal(Final20, SubsetIdxs2, ScoreMx2);
		Final20.insert(Idx);
		}

	vector<uint> Final20Vec;
	for (set<uint>::const_iterator iter = Final20.begin();
	  iter != Final20.end(); ++iter)
		Final20Vec.push_back(*iter);

	for (uint i = 0; i < 20; ++i)
		{
		uint Idxi = Final20Vec[i];
		char aai = Aminos[Idxi];
		const vector<double> Valuesi = FeatureValuesVec[Idxi];
		double SumScore = 0;
		Log("[%2u]", i);
		for (uint j = 0; j < 20; ++j)
			{
			uint Idxj = Final20Vec[j];
			char aaj = Aminos[Idxj];
			const vector<double> Valuesj = FeatureValuesVec[Idxj];
			for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
				{
				double Valuei = Valuesi[FeatureIndex];
				double Valuej = Valuesj[FeatureIndex];
				double Diff = XProf::GetDiff(FeatureIndex, Valuei, Valuej);
				uint Bin = XProf::GetFeatureBin(FeatureIndex, Diff);
				Bins[FeatureIndex] = Bin;
				}
			double Score = XProf::GetScore(aai, aaj, Bins);
			Log("  %7.1f", Score);
			SumScore += Score;
			}
		Log("  : %7.1f\n", SumScore/20);
		}
	}

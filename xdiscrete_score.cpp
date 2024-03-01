#include "myutils.h"
#include "xprof.h"

/***
Tests XProf::GetScore() on pairs of aligned structure profile positions from SCOP40.

palmscan2 \
  -xdiscrete_score d:/int/dave_grant/scop40/aligned_xfeatvecs_0.6_0.8.tsv \
  -log xdiscrete_score.log 

Output in log:

Mean true score 2.68, rand -5.44
***/

void cmd_xdiscrete_score()
	{
	FILE *f = OpenStdioFile(opt_xdiscrete_score);
	string HdrLine;
	vector<string> HdrFields;
	bool Ok = ReadLineStdioFile(f, HdrLine);
	asserta(Ok);
	Split(HdrLine, HdrFields, '\t');

	XProf::InitScoreTable();

	string Line;
	vector<string> Fields;
	vector<char> AminoQs;
	vector<char> AminoRs;
	vector<vector<uint> > BinsVec;
	double SumTrueScore = 0;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields[0]) == 1);
		asserta(SIZE(Fields[XFEATS+1]) == 1);
		char AminoQ = Fields[0][0];
		char AminoR = Fields[XFEATS+1][0];
		vector<uint> Bins;
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			double ValueQ = StrToFloat(Fields[FeatureIndex+1]);
			double ValueR = StrToFloat(Fields[XFEATS+FeatureIndex+2]);
			double Diff = XProf::GetDiff(FeatureIndex, ValueQ, ValueR);
			uint Bin = XProf::GetFeatureBin(FeatureIndex, Diff);
			Bins.push_back(Bin);
			}
		double Score = XProf::GetScore(AminoQ, AminoR, Bins);
		AminoQs.push_back(AminoQ);
		AminoRs.push_back(AminoR);
		BinsVec.push_back(Bins);
		SumTrueScore += Score;
		}
	const uint N = SIZE(AminoQs);
	asserta(SIZE(AminoRs) == N);
	asserta(SIZE(BinsVec) == N);

	double SumRandScore = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint idxq = randu32()%N;
		uint idxr = randu32()%N;
		char AminoQ = AminoQs[idxq];
		char AminoR = AminoRs[idxr];
		const vector<uint> &BinsQ = BinsVec[idxq];
		const vector<uint> &BinsR = BinsVec[idxr];
		vector<uint> Bins;
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			double ValueQ = StrToFloat(Fields[FeatureIndex+1]);
			double ValueR = StrToFloat(Fields[XFEATS+FeatureIndex+2]);
			double Diff = XProf::GetDiff(FeatureIndex, ValueQ, ValueR);
			uint Bin = XProf::GetFeatureBin(FeatureIndex, Diff);
			Bins.push_back(Bin);
			}
		double Score = XProf::GetScore(AminoQ, AminoR, Bins);
		SumRandScore += Score;
		}
	ProgressLog("Mean true score %.3g, rand %.3g\n",
	  SumTrueScore/N, SumRandScore/N);
	}

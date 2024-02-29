#include "myutils.h"
#include "xprof.h"

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
			uint BinQ = XProf::GetFeatureBin(FeatureIndex, ValueQ);
			uint BinR = XProf::GetFeatureBin(FeatureIndex, ValueR);
			Bins.push_back(BinQ);
			}
		double Score = XProf::GetScore(AminoQ, AminoR, Bins);
		Log("%c %c %.3g\n", AminoQ, AminoR, Score);
		}
	}

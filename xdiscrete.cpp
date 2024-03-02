#include "myutils.h"
#include "xprof.h"

/***
Discretize absolute differences between feature values in aligned position pairs.
Bin diffs into XBINS bins with equal numbers of aligned pairs.

Input
======
tsv with features of aligned pairs of residues
d:/int/dave_grant/scop40/aligned_xfeatvecs_0.6_0.8.tsv
Generated by 
  palmscan2 \
    -fa2xfeataln ../scop40/scop40.fa2 \
	-train_cal d:/int/scop40/out/domains_scop.cal \
	-mintm 0.6 \
	-maxtm 0.8 \
	-tsv d:/int/dave_grant/scop40/aligned_xfeatvecs_0.6_0.8.tsv

Ang_m2_p2, Ang_m3_p3 range from 0..179 degrees

# head d:/int/dave_grant/scop40/aligned_xfeatvecs_0.6_0.8.tsv
Q  Q_Ang_m2_p2  Q_Ang_m3_p3  Q_ED_p4  Q_ED_m4  Q_NU  Q_ND  R  R_Ang_m2_p2  R_Ang_m3_p3  R_ED_p4  R_ED_m4  R_NU  R_ND
E            0            0     11.7        0     0     0  S            0            0     12.9        0     0     0
R            0            0     9.04        0    14     4  K            0            0       13        0     5     8
D         17.5            0     8.06        0     2    14  P         23.1            0     13.1        0    12     9
N         83.1         62.8     9.57     11.7     0    14  L         13.9           22     11.6        0     1    19

***/

void cmd_xdiscrete()
	{
	asserta(optset_lo && optset_hi);
	uint ColQ = opt_lo;
	uint ColR = opt_hi;
	FILE *f = OpenStdioFile(opt_xdiscrete);
	string HdrLine;
	vector<string> HdrFields;
	bool Ok = ReadLineStdioFile(f, HdrLine);
	asserta(Ok);
	Split(HdrLine, HdrFields, '\t');
	const string &HdrQ = HdrFields[ColQ];
	const string &HdrR = HdrFields[ColR];
	asserta(StartsWith(HdrQ, "Q_"));
	asserta(StartsWith(HdrR, "R_"));
	const string &FeatureName = HdrQ.substr(2);

	bool Angle = false;
	bool NX = false;
	if (FeatureName == "Ang_m2_p2")
		Angle = true;
	else if (FeatureName == "Ang_m3_p3")
		Angle = true;
	else if (FeatureName == "ED_m4")
		;
	else if (FeatureName == "NU")
		NX = true;
	else if (FeatureName == "ND")
		NX = true;

	ProgressLog("feature %s, angle %c, NX %c\n",
	  FeatureName.c_str(), tof(Angle), tof(NX));
	asserta(HdrR.substr(2) == FeatureName);

	vector<double> ValueQs;
	vector<double> ValueRs;
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		const double ValueQ = StrToFloat(Fields[ColQ]);
		const double ValueR = StrToFloat(Fields[ColR]);
		ValueQs.push_back(ValueQ);
		ValueRs.push_back(ValueR);
		}
	const uint N = SIZE(ValueQs);
	ProgressLog("%s values\n", IntToStr(N));
	vector<double> Diffs;
	for (uint i = 0; i < N; ++i)
		{
		if (Angle)
			{
			double d1 = fabs(ValueQs[i] - ValueRs[i]);
			double d2 = fabs(fabs(ValueQs[i] - ValueRs[i]) - 180);
			double Diff = min(d1, d2);
			asserta(Diff <= 90);
			Diffs.push_back(Diff);
			}
		else
			{
			double Diff = fabs(ValueQs[i] - ValueRs[i]);
			Diffs.push_back(Diff);
			}
		}
	sort(Diffs.begin(), Diffs.end());
	vector<double> BinLos;
	if (NX)
		{
		for (uint i = 0; i < XBINS; ++i)
			BinLos.push_back(i);
		}
	else
		{
		for (uint i = 0; i < XBINS; ++i)
			{
			uint k = (N*i)/(XBINS+1);
			double dlo = Diffs[k];
			BinLos.push_back(dlo);
			ProgressLog("[%2u]  %.3g\n", i, dlo);
			}
		BinLos.push_back(DBL_MAX);
		}
	ProgressLog(" Max  %.3g\n", Diffs[N-1]);

	vector<uint> TrueCounts(10);
	for (uint i = 0; i < N; ++i)
		{
		double ValueQ = ValueQs[i];
		double ValueR = ValueRs[i];
		double Diff = 0;
		if (Angle)
			{
			double d1 = fabs(ValueQ - ValueR);
			double d2 = fabs(fabs(ValueQ - ValueR) - 180);
			Diff = min(d1, d2);
			asserta(Diff <= 90);
			}
		else
			Diff = fabs(ValueQ - ValueR);

		for (uint Bin = 0; Bin < XBINS; ++Bin)
			{
			if (Diff < BinLos[Bin+1])
				{
				TrueCounts[Bin] += 1;
				break;
				}
			}
		}

	vector<uint> RandCounts(10);
	for (uint i = 0; i < N; ++i)
		{
		uint idxq = randu32()%N;
		uint idxr = randu32()%N;
		double ValueQ = ValueQs[idxq];
		double ValueR = ValueRs[idxr];
		double Diff = 0;
		if (Angle)
			{
			double d1 = fabs(ValueQ - ValueR);
			double d2 = fabs(fabs(ValueQ - ValueR) - 180);
			Diff = min(d1, d2);
			asserta(Diff <= 90);
			}
		else
			Diff = fabs(ValueQ - ValueR);

		for (uint Bin = 0; Bin < XBINS; ++Bin)
			{
			if (Diff < BinLos[Bin+1])
				{
				RandCounts[Bin] += 1;
				break;
				}
			}
		}

	for (uint Bin = 0; Bin < XBINS; ++Bin)
		{
		uint nt = TrueCounts[Bin];
		uint nr = RandCounts[Bin];
		double TrueFreq = double(nt)/N;
		double RandFreq = double(nr)/N;
		double Score = TrueFreq < 1e-6 ? 0 : log(TrueFreq/RandFreq);
		ProgressLog("%10.10s [%2u]  %10.3g  %10u  %10u  %8.6f  %8.6f  %8.6f\n",
		  FeatureName.c_str(), BinLos[Bin], Bin, nt, nr, TrueFreq, RandFreq, Score);
		}

	for (uint Bin = 0; Bin < XBINS; ++Bin)
		{
		uint nt = TrueCounts[Bin];
		uint nr = RandCounts[Bin];
		double TrueFreq = double(nt)/N;
		double RandFreq = double(nr)/N;
		double Score = TrueFreq < 1e-6 ? 0 : log(TrueFreq/RandFreq);
		ProgressLog("	FeatureScoreBin(\"%s\", %u, %.3g, %.3g);\n",
		  FeatureName.c_str(),
		  Bin,
		  BinLos[Bin],
		  Score);
		}
	}
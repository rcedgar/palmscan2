#include "myutils.h"
#include "scop40bit.h"
#include "sort.h"

void cmd_scop40bit_roc()
	{
	asserta(optset_output);
	const string &FN = g_Arg1;
	SCOP40Bit SB;
	SB.ReadDomInfo();
	SB.ReadHits_Bin(FN);
	uint NT, NF;
	SB.CalcNXs(NT, NF);
	ProgressLog("NT %s\n", IntToStr(NT));
	ProgressLog("NF %s\n", IntToStr(NF));

	vector<double> ScoreSteps;
	vector<uint> TPCounts;
	vector<uint> FPCounts;
	SB.GetROCSteps(ScoreSteps, TPCounts, FPCounts);

	const uint StepCount = SIZE(ScoreSteps);
	asserta(SIZE(TPCounts) == StepCount);
	asserta(SIZE(FPCounts) == StepCount);
	FILE *fOut = CreateStdioFile(opt_output);
	double FPRAtTPR90 = DBL_MAX;
	double LastFPR = 0;
	for (uint StepIdx = 0; StepIdx < StepCount; ++StepIdx)
		{
		uint NTP = TPCounts[StepIdx];
		uint NFP = FPCounts[StepIdx];
		double TPR = double(NTP)/NT;
		double FPR = double(NFP)/NF;
		if (TPR > 0.9 && FPRAtTPR90 == DBL_MAX)
			FPRAtTPR90 = LastFPR;
		Log("%.8g\t%u\t%u\t%.6f\t%.6f\n", ScoreSteps[StepIdx],
		  NTP, NFP, TPR, FPR);
		fprintf(fOut, "%.8g\t%.6f\t%.6f\n",
		  ScoreSteps[StepIdx], TPR, FPR);
		LastFPR = FPR;
		}
	CloseStdioFile(fOut);
	ProgressLog("FPR %.3g at TPR=0.9\n", FPRAtTPR90);
	}

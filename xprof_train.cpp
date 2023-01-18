#include "myutils.h"
#include "pdbchain.h"
#include "xtrainer.h"
#include "xprof.h"
#include "outputfiles.h"

static char ScoreChar(double Score, double MaxScore)
	{
	asserta(Score <= MaxScore);
	if (Score <= 0)
		return ' ';

	if (Score >= MaxScore*0.75)
		return '*';
	if (Score >= MaxScore*0.5)
		return '+';
	if (Score >= MaxScore*0.25)
		return '.';
	if (Score >= MaxScore*0.1)
		return '_';
	return ' ';
	}

void cmd_xprof_train()
	{
	vector<PDBChain *> Chains;
	ReadChains(opt_xprof_train, Chains);

	XTrainer XT;
	XT.SetChains(Chains);

	const double MINTM = 0.6;
	const double MAXTM = 0.9;

	FILE *f = OpenStdioFile(opt_ref);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 7);
		double TM = StrToFloat(Fields[0]);
		if (TM < MINTM || TM > MAXTM)
			continue;
		const string &Label1 = Fields[1];
		const string &Label2 = Fields[2];
		const string &Row1 = Fields[5];
		const string &Row2 = Fields[6];

		XT.SetAln(Label1, Label2, Row1, Row2);
		}

	const uint FeatureCount = XProf::GetFeatureCount();
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		ProgressStep(FeatureIndex, FeatureCount, "Feature value dist");
		XT.SetValues(FeatureIndex);

		Log("\n");
		Log("Feature %s\n", XProf::GetFeatureName(FeatureIndex));
		XT.LogValueHist();
		}

	XT.Train();

	Log("\n\n");
	for (uint i = 0; i < FeatureCount; ++i)
		{
		ProgressStep(i, FeatureCount, "Correl");
		const char *Namei = XProf::GetFeatureName(i);
		for (uint j = i+1; j < FeatureCount; ++j)
			{
			const char *Namej = XProf::GetFeatureName(j);
			double r = XT.CorrelFeatures(i, j);
			Log("Correl  %10.10s  %10.10s  %6.3f\n", Namei, Namej, r);
			}
		}

	XT.LogFinal();

	vector<vector<double> > Mx;
	XT.GetScoreMx(0, 1, Mx);

	const uint L0 = XT.m_Chains[0]->GetSeqLength();
	const uint L1 = XT.m_Chains[1]->GetSeqLength();
	double MaxScore = 0;
	for (uint i = 0; i < L0; ++i)
		for (uint j = 0; j < L1; ++j)
			MaxScore = max(MaxScore, Mx[i][j]);

	if (g_ftsv_dot != 0)
		{
		FILE *f = g_ftsv_dot;
		fprintf(f, "%u\t%u\n", L0, L1);
		for (uint i = 0; i < L0; ++i)
			{
			for (uint j = 0; j < L1; ++j)
				{
				double Score = Mx[i][j];
				fprintf(f, "%.3f\n", Score/MaxScore);
				}
			}
		}
	}

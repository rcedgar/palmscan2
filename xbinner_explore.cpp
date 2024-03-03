#include "myutils.h"
#include "xprof.h"
#include "alpha.h"
#include "xbinner.h"
#include "x2data.h"

void cmd_xbinner_explore()
	{
	XBinner::InitCentroids();
	XProf::InitScoreTable();

	X2Data X2;
	X2.FromTsv(opt_xbinner_explore);

	vector<double> Freqs;
	vector<vector<double> > FreqMx;
	XBinner XB;
	XB.GetFreqs(X2, Freqs, FreqMx);

	vector<vector<double> > ScoreMx;
	double ExpScore = XBinner::GetLogOddsMx(Freqs, FreqMx, ScoreMx);

	ProgressLog("Expected score %.3f\n", ExpScore);
	}

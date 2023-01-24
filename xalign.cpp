#include "myutils.h"
#include "mx.h"
#include "xprof.h"
#include "xprofdata.h"
#include "xtrainer.h"
#include "outputfiles.h"
#include "omplock.h"

float SW(const Mx<float> &SMx, 
  Mx<float> &a_FwdM,
  Mx<float> &a_FwdD,
  Mx<float> &a_FwdI,
  Mx<char> &a_TBM,
  Mx<char> &a_TBD,
  Mx<char> &a_TBI,
  uint &Starti, uint &Startj, string &Path);

static uint g_QueryCount;
static uint g_HitCount;

static double GetScorePosPair(
  const XProfData &XD1, uint Pos1, 
  const XProfData &XD2, uint Pos2,
  const vector<vector<vector<double> > > &FeatureIndexToLogOddsMx)
	{
	const uint FeatureCount = XProf::GetFeatureCount();

	double Sum = 0;
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		uint iValue1 = XD1.GetIntFeature(FeatureIndex, Pos1);
		uint iValue2 = XD2.GetIntFeature(FeatureIndex, Pos2);
		if (iValue1 != UINT_MAX && iValue2 != UINT_MAX)
			{
			double Score =
			  FeatureIndexToLogOddsMx[FeatureIndex][iValue1][iValue2];
			Sum += Score;
			}
		}
	return Sum/FeatureCount;
	}

double XAlign(Mx<float> &SMx, 
  Mx<float> &a_FwdM,
  Mx<float> &a_FwdD,
  Mx<float> &a_FwdI,
  Mx<char> &a_TBM,
  Mx<char> &a_TBD,
  Mx<char> &a_TBI,
  const XProfData &Prof1, const XProfData &Prof2,
  const vector<vector<vector<double> > > &FeatureIndexToLogOddsMx,
  double MinScore)
	{
	const uint L1 = Prof1.GetSeqLength();
	const uint L2 = Prof2.GetSeqLength();
	SMx.Alloc("SMx", L1, L2);

	for (uint Pos1 = 0; Pos1 < L1; ++Pos1)
		{
		for (uint Pos2 = 0; Pos2 < L2; ++Pos2)
			{
			float Score = (float)
			  GetScorePosPair(Prof1, Pos1, Prof2, Pos2,
				FeatureIndexToLogOddsMx);
			SMx.Put(Pos1, Pos2, Score);
			}
		}

	uint Start0, Start1;
	string Path;
	float Score = SW(SMx, a_FwdM, a_FwdD, a_FwdI,
	  a_TBM, a_TBD, a_TBI, Start0, Start1, Path);
	if (Score < MinScore)
		return Score;
	++g_HitCount;
	uint Cols = SIZE(Path);

	const string &A = Prof1.m_Seq;
	const string &B = Prof2.m_Seq;
	const uint ColCount = SIZE(Path);
	uint i = Start0;
	uint j = Start1;
	string ARow;
	string BRow;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == 'D')
			{
			ARow += A[i];
			++i;
			}
		else
			ARow += '-';

		if (c == 'M' || c == 'I')
			{
			BRow += B[j];
			++j;
			}
		else
			BRow += '-';
		}

	if (g_ftsv)
		{
		Lock();
		FILE *f = g_ftsv;
		fprintf(f, "%.1f", Score);
		fprintf(f, "\t%s", Prof1.m_Label.c_str());
		fprintf(f, "\t%s", Prof2.m_Label.c_str());
		fprintf(f, "\t%u", ColCount);
		fprintf(f, "\t%.2f", Score/ColCount);
		if (!opt_norows)
			{
			fprintf(f, "\t%s", ARow.c_str());
			fprintf(f, "\t%s", BRow.c_str());
			}
		fprintf(f, "\n");
		Unlock();
		}
	return Score;
	}

static void Thread(FILE *fQ, double MinScore,
  const vector<XProfData *> &DBXDs,
  const vector<vector<vector<double> > > &FeatureIndexToLogOddsMx)
	{
	const uint DBN = SIZE(DBXDs);
	Mx<float> SMx;
	Mx<float> a_FwdM;
	Mx<float> a_FwdD;
	Mx<float> a_FwdI;
	Mx<char> a_TBM;
	Mx<char> a_TBD;
	Mx<char> a_TBI;

	XProfData QP;
	for (;;)
		{
		bool Ok;
#pragma omp critical
		{
		Ok = QP.FromCfv(fQ);
		}

		if (!Ok)
			return;

#pragma omp critical
		{
		++g_QueryCount;
		if (g_QueryCount%1 == 0)
			Progress("%u queries, %u hits\r", g_QueryCount, g_HitCount);
		}

		for (uint i = 0; i < DBN; ++i)
			{
			const XProfData &DP = *DBXDs[i];
			XAlign(SMx, a_FwdM, a_FwdD, a_FwdI,
			  a_TBM, a_TBD, a_TBI,
				QP, DP, FeatureIndexToLogOddsMx, MinScore);
			}
		}
	}

void cmd_xalign()
	{
	const string &QueryFileName = opt_xalign;
	const string &DBFileName = opt_db;

	double MinScore = 100;
	if (optset_minscore)
		MinScore = opt_minscore;

	vector<XProfData *> DBXDs;
	Progress("Reading db... ");
	FILE *fDB = OpenStdioFile(DBFileName);
	for (;;)
		{
		XProfData *XD = new XProfData;
		bool Ok = XD->FromCfv(fDB);
		if (!Ok)
			break;
		DBXDs.push_back(XD);
		}
	CloseStdioFile(fDB);
	fDB = 0;
	const uint DBN = SIZE(DBXDs);
	Progress(" %u profiles\n", DBN);

	vector<vector<vector<double> > > FeatureIndexToLogOddsMx;
	XTrainer::LogOddsFromTsv(opt_logodds, FeatureIndexToLogOddsMx);

	XProfData QP;
	FILE *fQ = OpenStdioFile(QueryFileName);

	uint ThreadCount = GetRequestedThreadCount();
#pragma omp parallel num_threads(ThreadCount)
	Thread(fQ, MinScore, DBXDs, FeatureIndexToLogOddsMx);

	ProgressLog("%u queries, %u hits\n", g_QueryCount, g_HitCount);
	}

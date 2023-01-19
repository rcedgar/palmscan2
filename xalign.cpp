#include "myutils.h"
#include "mx.h"
#include "xprof.h"
#include "xprofdata.h"
#include "xtrainer.h"

float SW(const Mx<float> &SMx, uint &Starti, uint &Startj, string &Path);

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

static void XAlign(const XProfData &Prof1, const XProfData &Prof2,
  const vector<vector<vector<double> > > &FeatureIndexToLogOddsMx)
	{
	Mx<float> SMx;
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
	float Score = SW(SMx, Start0, Start1, Path);
	uint Cols = SIZE(Path);
	Log("Score = %.1f, Cols=%u, Path = %s\n", Score, Cols, Path.c_str());

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
	Log("A  >%s\n", Prof1.m_Label.c_str());
	Log("B  >%s\n", Prof2.m_Label.c_str());
	Log("A  %s\n", ARow.c_str());
	Log("B  %s\n", BRow.c_str());
	}

void cmd_xalign()
	{
	const string &QueryFileName = opt_xalign;
	const string &DBFileName = opt_db;

	vector<XProfData *> DBXDs;
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

	vector<vector<vector<double> > > FeatureIndexToLogOddsMx;
	XTrainer::LogOddsFromTsv(opt_logodds, FeatureIndexToLogOddsMx);

	XProfData QP;
	FILE *fQ = OpenStdioFile(QueryFileName);
	for (;;)
		{
		bool Ok = QP.FromCfv(fQ);
		if (!Ok)
			break;

		Progress("%s\r", QP.m_Label.c_str());
		for (uint i = 0; i < DBN; ++i)
			{
			const XProfData &DP = *DBXDs[i];
			XAlign(QP, DP, FeatureIndexToLogOddsMx);
			}
		}
	Progress("\n");
	}

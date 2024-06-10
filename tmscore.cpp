#include "myutils.h"
#include "tma.h"
#include "seqdb.h"
#include "outputfiles.h"

double TMA::CalcTMScore(const PDBChain &Q, const PDBChain &R,
  const string &RowQ, const string &RowR)
	{
	m_Q = &Q;
	m_R = &R;

	const uint QL = Q.GetSeqLength();
	const uint RL = R.GetSeqLength();

	char *seqx = new char[QL+1];
	char *seqy = new char[RL+1];

	memcpy(seqx, Q.m_Seq.c_str(), QL+1);
	memcpy(seqy, R.m_Seq.c_str(), RL+1);

	double **xa = 0;
	double **ya = 0;
	NewArray(&xa, QL, 3);
	NewArray(&ya, RL, 3);

	int xlen = (int) QL;
	int ylen = (int) RL;

	for (int i = 0; i < xlen; ++i)
		{
		xa[i][0] = Q.m_Xs[i];
		xa[i][1] = Q.m_Ys[i];
		xa[i][2] = Q.m_Zs[i];
		}

	for (int i = 0; i < ylen; ++i)
		{
		ya[i][0] = R.m_Xs[i];
		ya[i][1] = R.m_Ys[i];
		ya[i][2] = R.m_Zs[i];
		}

	double d0A = DBL_MAX;
	double d0B = DBL_MAX;
	return -1;
	}

void cmd_tmscore()
	{
	const string &QueryFileName = opt_tmscore;
	const string &RefFileName = opt_ref;
	const string &AlnFileName = opt_input;

	SeqDB Aln;
	Aln.FromFasta(AlnFileName, true);
	asserta(Aln.IsAligned());
	Aln.SetLabelToIndex();

	TMA T;

	PDBChain Q;
	PDBChain R;
	Q.FromCal(QueryFileName);
	R.FromCal(RefFileName);

	const string &QLabel = Q.m_Label;
	const string &RLabel = R.m_Label;

	string RowQ;
	string RowR;
	Aln.GetSeqByLabel(QLabel, RowQ);
	Aln.GetSeqByLabel(RLabel, RowR);
	asserta(SIZE(RowQ) == SIZE(RowR));

	string QAcc;
	string RAcc;
	Q.GetAcc(QAcc);
	R.GetAcc(RAcc);

	double TM = T.AlignChains(Q, R);
	ProgressLog("%.4f  %s  %s\n", TM, QAcc.c_str(), RAcc.c_str());

	double TM2 = T.CalcTMScore(Q, R, RowQ, RowR);
	ProgressLog("TM2=%.4f\n", TM2);
	}

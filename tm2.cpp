#include "myutils.h"
#include "tma.h"
#include "outputfiles.h"

void cmd_tm2()
	{
	const string &QueryFileName = opt_tm2;
	const string &RefFileName = opt_ref;

	TMA T;

	PDBChain Q;
	PDBChain R;
	Q.FromCal(QueryFileName);
	R.FromCal(RefFileName);

	string QAcc;
	string RAcc;
	Q.GetAcc(QAcc);
	R.GetAcc(RAcc);

	double TM = T.AlignChains(Q, R);
	ProgressLog("%.4f  %s  %s\n", TM, QAcc.c_str(), RAcc.c_str());
	}

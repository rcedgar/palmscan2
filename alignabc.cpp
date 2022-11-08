#include "myutils.h"
#include "pdb.h"
#include "trisearcher.h"
#include "abcxyz.h"

void cmd_alignabc()
	{
	const string &QueryFN = opt_alignabc;
	const string &RefFN = opt_ref;

	PDB Q;
	PDB R;

	Q.FromFile(QueryFN);
	Q.LogMe();

	R.FromFile(RefFN);
	R.LogMe();

	TriSearcher TS;
	TS.Radius = 1.5;
	TS.MaxTriRMSD2 = 1.5;
	TS.NABmin = 10;
	TS.NABmax = 80;
	TS.NBCmin = 10;
	TS.NBCmax = 80;
	TS.NACmin = 80;
	TS.NACmax = 200;

	TS.Search(Q, R);
	TS.SetTriForm();
	TS.LogMe();

	//uint QPosA, QPosB, QPosC;
	//double MotifRMSD2;
	//bool Found = TS.GetTopHit(QPosA, QPosB, QPosC, MotifRMSD2);
	//asserta(Found);

	//string Path;
	//TS.AlignPalm(QPosA, QPosB, QPosC, Path);
	}

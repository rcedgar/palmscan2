#include "myutils.h"
#include "tma.h"
#include "seqdb.h"
#include "outputfiles.h"

int TMA::TMalign_main_score(
	const string &rowa, const string &rowb,
	double** xa, double** ya,
	const char* seqx, const char* secy,
	double& TM1, double& TM2,
	double& d0A, double& d0B,
	string& seqM, string& seqxA, string& seqyA,
	const int xlen, const int ylen)
	{
	double D0_MIN;        //for d0
	double Lnorm;         //normalization length
	double score_d8, d0, d0_search, dcu0;//for TMscore search
	double t[3], u[3][3]; //Kabsch translation vector and rotation matrix
	double** xtm, ** ytm;  // for TMscore search engine
	double** xt;          //for saving the superposed version of r_1 or xtm
	double** r1, ** r2;    // for Kabsch rotation

	/***********************/
	/* allocate memory     */
	/***********************/
	int minlen = min(xlen, ylen);
	//NewArray(&val, xlen+1, ylen+1);
	NewArray(&xtm, minlen, 3);
	NewArray(&ytm, minlen, 3);
	NewArray(&xt, xlen, 3);
	NewArray(&r1, minlen, 3);
	NewArray(&r2, minlen, 3);

	/***********************/
	/*    parameter set    */
	/***********************/
	parameter_set4search(xlen, ylen, D0_MIN, Lnorm,
		score_d8, d0, d0_search, dcu0);
	if (0)
		{
		Log("parameter_set4search:\n");
		Log("   D0_MIN  %.3g\n", D0_MIN);
		Log("   Lnorm   %.3g\n", Lnorm);
		Log("score_d8   %.3g\n", score_d8);
		Log("      d0   %.3g\n", d0);
		Log("d0_search  %.3g\n", d0_search);
		Log("    dcu0   %.3g\n", dcu0);
		}

	int i;
	int* invmap0 = new int[ylen + 1];
	int* invmap = new int[ylen + 1];
	double TM, TMmax = -1;
	for (i = 0; i < ylen; i++) invmap0[i] = -1;

	double ddcc = 0.4;
	if (Lnorm <= 40) ddcc = 0.1;   //Lnorm was setted in parameter_set4search
	double local_d0_search = d0_search;

	for (int j = 0; j < ylen; j++)// Set aligned position to be "-1"
		invmap[j] = -1;

	int i1 = -1;// in C version, index starts from zero, not from one
	int i2 = -1;
	int L1 = (int) rowa.size();
	int L2 = (int) rowb.size();
	int L = min(L1, L2);// Get positions for aligned residues
	for (int kk1 = 0; kk1 < L; kk1++)
		{
		if (rowa[kk1] != '-') i1++;
		if (rowb[kk1] != '-')
			{
			i2++;
			if (i2 >= ylen || i1 >= xlen) kk1 = L;
			else if (rowa[kk1] != '-') invmap[i2] = i1;
			}
		}

	//--------------- 2. Align proteins from original alignment
	double prevD0_MIN = D0_MIN;// stored for later use
	double prevLnorm = Lnorm;
	double prevd0 = d0;
	double rmsd_ali = DBL_MAX;
	int L_ali = INT_MAX;
	double TM_ali = standard_TMscore(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
		invmap, L_ali, rmsd_ali, D0_MIN, Lnorm, d0, d0_search, score_d8,
		t, u, -2);
	D0_MIN = prevD0_MIN;
	Lnorm = prevLnorm;
	d0 = prevd0;
	int simplify_step = 40;
	int score_sum_method = 8;
	TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
		invmap, t, u, simplify_step, score_sum_method,
		local_d0_search, true, Lnorm, score_d8, d0);
	if (TM > TMmax)
		{
		TMmax = TM;
		for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
		}

	//check if the initial alignment is generated approriately
	bool flag = false;
	for (i = 0; i < ylen; i++)
		{
		if (invmap0[i] >= 0)
			{
			flag = true;
			break;
			}
		}
	if (!flag)
		Die("No alignment");

	simplify_step = 40;
	score_sum_method = 8;
	TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
		invmap0, t, u, simplify_step, score_sum_method, local_d0_search,
		false, Lnorm, score_d8, d0);

	//select pairs with dis<d8 for final TMscore computation and output alignment
	int k = 0;
	int* m1, * m2;
	double d;
	m1 = new int[xlen]; //alignd index in x
	m2 = new int[ylen]; //alignd index in y
	do_rotation(xa, xt, xlen, t, u);
	k = 0;
	int n_ali = 0;
	for (int j = 0; j < ylen; j++)
		{
		i = invmap0[j];
		if (i >= 0)//aligned
			{
			n_ali++;
			d = sqrt(dist(&xt[i][0], &ya[j][0]));
			m1[k] = i;
			m2[k] = j;

			xtm[k][0] = xa[i][0];
			xtm[k][1] = xa[i][1];
			xtm[k][2] = xa[i][2];

			ytm[k][0] = ya[j][0];
			ytm[k][1] = ya[j][1];
			ytm[k][2] = ya[j][2];

			r1[k][0] = xt[i][0];
			r1[k][1] = xt[i][1];
			r1[k][2] = xt[i][2];
			r2[k][0] = ya[j][0];
			r2[k][1] = ya[j][1];
			r2[k][2] = ya[j][2];

			k++;
			}
		}
	int n_ali8 = k;

	double rmsd;
	double Lnorm_0 = ylen;

	//normalized by length of structure A
	parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, -2);
	d0A = d0;
	double d0_0 = d0A;
	local_d0_search = d0_search;
	double t0[3], u0[3][3];
	{
	const int simplify_step = 1;
	const int score_sum_method = 0;
	TM1 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
		simplify_step, score_sum_method, &rmsd, local_d0_search,
		Lnorm, score_d8, d0);
	}
	double TM_0 = TM1;

	//normalized by length of structure B
	parameter_set4final(xlen + 0.0, D0_MIN, Lnorm, d0, d0_search, -2);
	d0B = d0;
	local_d0_search = d0_search;
	{
	const int simplify_step = 1;
	const int score_sum_method = 0;
	TM2 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t, u,
		simplify_step, score_sum_method, &rmsd, local_d0_search,
		Lnorm, score_d8, d0);
	}

	return 0; // zero for no exception
	}

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
	string seqM;
	int iResult = TMalign_main_score(
		RowQ, RowR,
	  xa, ya, seqx, seqy, m_TM1, m_TM2, d0A, d0B,
	  seqM, m_QRow, m_RRow, xlen, ylen);

	return m_TM1;
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

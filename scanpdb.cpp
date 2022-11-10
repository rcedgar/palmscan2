#include "myutils.h"
#include "pdbchain.h"
#include "rdrpmodel.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"

void cmd_scanpdb()
	{
	const string &QueryFN = opt_scanpdb;

	Die("TODO");
	//PDBChain Q;
	//Q.FromFile(QueryFN);
	//Q.LogMe();

	//const string &QSeq = Q.m_Seq;
	//const string &QLabel = Q.m_Label;

	//asserta(optset_model);
	//const string &ModelFileName = opt_model;
	//RdRpModel Model;
	//Model.FromModelFile(ModelFileName);

	//RdRpSearcher RS;
	//RS.Init(Model);

	//RdRpSearcher::InitOutput();
	//RS.Search(QLabel, QSeq);
	//RS.WriteOutput();
	//RdRpSearcher::CloseOutput();

	//uint APos = RS.GetMotifPos(0);
	//uint BPos = RS.GetMotifPos(1);
	//uint CPos = RS.GetMotifPos(2);

	//if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
	//	Die("Missing motif(s)");

	//Q.m_MotifPosVec.clear();
	//Q.m_MotifPosVec.push_back(APos);
	//Q.m_MotifPosVec.push_back(BPos);
	//Q.m_MotifPosVec.push_back(CPos);

	//vector<vector<double> > MotifCoords;
	//Q.GetMotifCoords(MotifCoords);
	//LogMx("MotifCoords", MotifCoords);

	//vector<double> t;
	//vector<vector<double> > R;
	//GetTriForm(MotifCoords, t, R);

	//LogVec("t", t);
	//LogMx("R", R);

	//Q.ToCal(opt_calout);

	//const uint QL = Q.GetSeqLength();
	//vector<double> Pt;
	//vector<double> XPt;
	//for (uint Pos = 0; Pos < QL; ++Pos)
	//	{
	//	Q.GetPt(Pos, Pt);
	//	XFormPt(Pt, t, R, XPt);
	//	Q.SetPt(Pos, XPt);
	//	}
	//Q.ToCal(opt_caloutx);
	//Q.ToPDB(opt_pdboutx);
	}

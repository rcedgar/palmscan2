#include "myutils.h"
#include "pdbchain.h"
#include "rdrpmodel.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "calreader.h"
#include "omplock.h"

void ReadPDBs(const string &FileName, vector<PDBChain *> &Structures);

void cmd_search3d_pssms()
	{
	const string &QueryFN = opt_search3d_pssms;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;
	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

	RdRpSearcher::InitOutput();
	uint ThreadCount = GetRequestedThreadCount();

	vector<vector<PDBChain *> > QVecs(ThreadCount);
	
	uint ThreadFinishedCount = 0;
	uint DoneCount = 0;
	uint HitCount = 0;
	vector<PDBChain> Qs(ThreadCount);
	vector<bool> ThreadDone(ThreadCount);
	vector<RdRpSearcher> RSs(ThreadCount);
	for (uint i = 0; i < ThreadCount; ++i)
		RSs[i].Init(Model);

	CalReader CR;
	CR.Open(QueryFN);

#pragma omp parallel num_threads(ThreadCount)
	for (;;)
		{
		uint ThreadIndex = GetThreadIndex();
		PDBChain& Q = Qs[ThreadIndex];
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			{
			Lock("Done");
			if (!ThreadDone[ThreadIndex])
				{
				ThreadDone[ThreadIndex] = true;
				++ThreadFinishedCount;
				}
			Unlock("Done");
			if (ThreadFinishedCount == ThreadCount)
				break;
			}
		Lock("Done");
		++DoneCount;
		Unlock("Done");
		if (DoneCount%100 == 0)
			{
			Lock("Progress");
			string sPct;
			CR.GetPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r", sPct, HitCount, DoneCount);
			Unlock("Progress");
			}

		const string &QSeq = Q.m_Seq;
		string QLabel;
		Q.GetLabel(QLabel);

		RdRpSearcher &RS = RSs[ThreadIndex];
		RS.Search(QLabel, QSeq);
		LockOutput();
		RS.WriteOutput();
		UnlockOutput();

		uint APos = RS.GetMotifPos(0);
		uint BPos = RS.GetMotifPos(1);
		uint CPos = RS.GetMotifPos(2);

		if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
			{
			Log("Motif(s) not found in >%s\n", QLabel.c_str());
			continue;
			}
		if (CPos < APos)
			{
			Log("Permuted domain >%s\n", QLabel.c_str());
			continue;
			}

		Lock("Done");
		++HitCount;
		Unlock("Done");

		Q.m_MotifPosVec.clear();
		Q.m_MotifPosVec.push_back(APos);
		Q.m_MotifPosVec.push_back(BPos);
		Q.m_MotifPosVec.push_back(CPos);

		vector<vector<double> > MotifCoords;
		Q.GetMotifCoords(MotifCoords);
//		LogMx("MotifCoords", MotifCoords);

		vector<double> t;
		vector<vector<double> > R;
		GetTriForm(MotifCoords, t, R);

		//LogVec("t", t);
		//LogMx("R", R);

		const uint QL = Q.GetSeqLength();
		vector<double> Pt;
		vector<double> XPt;
		for (uint Pos = 0; Pos < QL; ++Pos)
			{
			Q.GetPt(Pos, Pt);
			XFormPt(Pt, t, R, XPt);
			Q.SetPt(Pos, XPt);
			}
		uint PPL = CPos + CL - APos;
		LockOutput();
		Q.ToCalSeg(g_fppc, APos, PPL);
		UnlockOutput();
		}
	Progress("100.0%% done, %u / %u hits\r", HitCount, DoneCount);
	}

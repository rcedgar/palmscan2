#include "myutils.h"
#include "pdbchain.h"
#include "trisearcher.h"
#include "abcxyz.h"
#include "sort.h"
#include "tshitmgr.h"
#include "calreader.h"
#include "omplock.h"

void GetFileNames(const string &SpecFileName, vector<string> &FileNames);
void ReadPDBs(const vector<string> &FileNames, vector<PDBChain *> &Structures);
void InitTS(TriSearcher &TS);
void Search1(TriSearcher &TS, TSHitMgr &HM,
  PDBChain &Q, vector<PDBChain *> &RefPDBs);

void cmd_search3dcal()
	{
	const string &QueryFN = opt_search3dcal;
	const string &RefFN = opt_ref;

	CalReader CR;
	CR.Open(QueryFN);

	vector<string> RefFileNames;
	GetFileNames(RefFN, RefFileNames);
	Progress("Read reference structures...");
	vector<PDBChain *> RefPDBs;
	ReadPDBs(RefFileNames, RefPDBs);
	for (uint i = 0; i < SIZE(RefPDBs); ++i)
		RefPDBs[i]->SetSS();
	Progress(" done.\n");

	const uint RefN = SIZE(RefPDBs);

	uint ThreadCount = GetRequestedThreadCount();
	vector<vector<PDBChain *> > QVecs(ThreadCount);
	vector<TriSearcher> TSs(ThreadCount);
	vector<TSHitMgr> HMs(ThreadCount);
	for (uint i = 0; i < ThreadCount; ++i)
		InitTS(TSs[i]);

	TriSearcher::OpenOutputFiles();
	TSHitMgr::OpenOutputFiles();

	uint ThreadFinishedCount = 0;
	uint DoneCount = 0;
	uint HitCount = 0;
	vector<PDBChain> Qs(ThreadCount);
	vector<bool> ThreadDone(ThreadCount);

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
		++DoneCount;
		if (DoneCount%100 == 0)
			{
			Lock("Progress");
			Progress("%u done, %u hits\r", DoneCount, HitCount);
			Unlock("Progress");
			}
		TriSearcher &TS = TSs[ThreadIndex];
		TSHitMgr &HM = HMs[ThreadIndex];

		Q.SetSS();
		Search1(TS, HM, Q, RefPDBs);
		if (!HM.m_Hits.empty())
			++HitCount;
		}
	Progress("%u done, %u hits\n", DoneCount, HitCount);

	TriSearcher::CloseOutputFiles();
	TSHitMgr::CloseOutputFiles();
	}

#include "myutils.h"
#include "pdbchain.h"
#include "trisearcher.h"
#include "abcxyz.h"
#include "sort.h"
#include "tshitmgr.h"
#include "omplock.h"

void GetFileNames(const string &SpecFileName, vector<string> &FileNames);
void ReadPDBs(const vector<string> &FileNames, vector<PDBChain *> &Structures);
void InitTS(TriSearcher& TS);

void Search1(TriSearcher &TS, TSHitMgr &HM,
  PDBChain &Q, vector<PDBChain *> &RefPDBs)
	{
	HM.SetQuery(Q);

	vector<uint> PosAs;
	vector<uint> PosBs;
	vector<uint> PosCs;
	vector<double> MotifRMSD2s;
	vector<string> RefLabels;
	const uint RefN = SIZE(RefPDBs);
	for (uint iR = 0; iR < RefN; ++iR)
		{
		const PDBChain &R = *RefPDBs[iR];

		if (opt_self && Q.m_ChainLabel == R.m_ChainLabel)
			continue;

		TS.Search(Q, R);
		TS.WriteOutput();
		TSHit TH;
		bool Ok = TS.GetTopHit(TH);
		if (Ok)
			HM.Add(TH);
		}

	HM.WriteOutput();
	}

void cmd_search3d()
	{
	const string &QueryFN = opt_search3d;
	const string &RefFN = opt_ref;

	vector<string> QueryFileNames;
	vector<string> RefFileNames;
	GetFileNames(QueryFN, QueryFileNames);
	GetFileNames(RefFN, RefFileNames);

	Progress("Read reference structures...");
	vector<PDBChain *> RefPDBs;
	ReadPDBs(RefFileNames, RefPDBs);
	for (uint i = 0; i < SIZE(RefPDBs); ++i)
		RefPDBs[i]->SetSS();
	Progress(" done.\n");

	const uint QueryN = SIZE(QueryFileNames);
	const uint RefN = SIZE(RefPDBs);

	uint ThreadCount = GetRequestedThreadCount();
	vector<vector<PDBChain *> > QVecs(ThreadCount);
	vector<TriSearcher> TSs(ThreadCount);
	vector<TSHitMgr> HMs(ThreadCount);
	for (uint i = 0; i < ThreadCount; ++i)
		InitTS(TSs[i]);

	TriSearcher::OpenOutputFiles();
	TSHitMgr::OpenOutputFiles();

	uint DoneCount = 0;

#pragma omp parallel for num_threads(ThreadCount)
	for (int iQ = 0; iQ < (int) QueryN; ++iQ)
		{
		Lock("iQloop");
		ProgressStep(DoneCount++, QueryN, "Searching");
		Unlock("iQloop");

		uint ThreadIndex = GetThreadIndex();
		vector<PDBChain *> &QVec = QVecs[ThreadIndex];
		TriSearcher &TS = TSs[ThreadIndex];
		TSHitMgr &HM = HMs[ThreadIndex];

		const string &QueryFileName = QueryFileNames[iQ];
		PDBChain::ReadChainsFromFile(QueryFileName, QVec);
		uint ChainCount = SIZE(QVec);
		for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
			{
			PDBChain &Q = *QVec[ChainIndex];
			Q.SetSS();
			Search1(TS, HM, Q, RefPDBs);
			}
		}

	TriSearcher::CloseOutputFiles();
	TSHitMgr::CloseOutputFiles();
	}

#include "myutils.h"
#include "pdb.h"
#include "trisearcher.h"
#include "abcxyz.h"
#include "sort.h"
#include "omplock.h"

void GetFileNames(const string &SpecFileName, vector<string> &FileNames);
void ReadPDBs(const vector<string> &FileNames, vector<PDB *> &Structures);

void InitTS(TriSearcher &TS)
	{
	TS.Radius = 1.5;
	TS.MaxTriRMSD2 = 1.5;
	TS.NABmin = 10;
	TS.NABmax = 80;
	TS.NBCmin = 10;
	TS.NBCmax = 80;
	TS.NACmin = 80;
	TS.NACmax = 200;
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
	vector<PDB *> RefPDBs;
	ReadPDBs(RefFileNames, RefPDBs);
	Progress(" done.\n");

	const uint QueryN = SIZE(QueryFileNames);
	const uint RefN = SIZE(RefPDBs);

	uint ThreadCount = GetRequestedThreadCount();
	vector<PDB> Qs(ThreadCount);
	vector<TriSearcher> TSs(ThreadCount);
	for (uint i = 0; i < ThreadCount; ++i)
		InitTS(TSs[i]);

	uint DoneCount = 0;

#pragma omp parallel for num_threads(ThreadCount)
	for (int iQ = 0; iQ < (int) QueryN; ++iQ)
		{
		Lock();
		ProgressStep(DoneCount++, QueryN, "Searching");
		Unlock();

		uint ThreadIndex = GetThreadIndex();
		PDB &Q = Qs[ThreadIndex];
		TriSearcher &TS = TSs[ThreadIndex];

		const string &QueryFileName = QueryFileNames[iQ];
		Q.FromFile(QueryFileName);

		vector<uint> PosAs;
		vector<uint> PosBs;
		vector<uint> PosCs;
		vector<double> MotifRMSD2s;
		vector<string> RefLabels;
		for (uint iR = 0; iR < RefN; ++iR)
			{
			const PDB &R = *RefPDBs[iR];

			TS.Search(Q, R);
#if 0
			Log("\n");
			Log("__________________________________\n");
			Log("Q>%s\n", Q.m_Label.c_str());
			Log("R>%s\n", R.m_Label.c_str());
			TS.LogMe();
#endif
			uint PosA = UINT_MAX;
			uint PosB = UINT_MAX;
			uint PosC = UINT_MAX;
			double MotifRMSD2 = DBL_MAX;
			bool Found = TS.GetTopHit(PosA, PosB, PosC, MotifRMSD2);
			if (Found)
				{
				PosAs.push_back(PosA);
				PosBs.push_back(PosB);
				PosCs.push_back(PosC);
				MotifRMSD2s.push_back(MotifRMSD2);
				RefLabels.push_back(R.m_Label);
				}
			}

		const uint N = SIZE(MotifRMSD2s);
		if (N == 0)
			{
			Lock();
			Log("\n");
			Log("Q>%s  %u hits\n", Q.m_Label.c_str(), N);
			Unlock();
			continue;
			}

		vector<uint> HitOrder;
		HitOrder.resize(N);
		QuickSortOrder(MotifRMSD2s.data(), N, HitOrder.data());

#if 1
		Lock();
		Log("\n");
		Log("Q>%s  %u hits\n", Q.m_Label.c_str(), N);
		Log(" PosA   PosB   PosC    MotifD2\n");
		for (uint k = 0; k < min(N, 4u); ++k)
			{
			uint i = HitOrder[k];
			uint PosA = PosAs[i];
			uint PosB = PosBs[i];
			uint PosC = PosCs[i];
			double MotifRMSD2 = MotifRMSD2s[i];
			const string &RefLabel = RefLabels[i];

			string A, B, C;
			Q.GetMotifSeq(PosA, AL, false, A);
			Q.GetMotifSeq(PosB, BL, false, B);
			Q.GetMotifSeq(PosC, CL, false, C);

			Log("%5u", PosA);
			Log("  %5u", PosB);
			Log("  %5u", PosC);
			Log("  %9.1f", MotifRMSD2);
			Log("  %s", A.c_str());
			Log("  %s", B.c_str());
			Log("  %s", C.c_str());
			Log("  %s", RefLabel.c_str());
			Log("\n");
			}
		Unlock();
#endif
		}
	}

#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "sort.h"
#include "calreader.h"
#include "ppcaligner.h"
#include "searchparams.h"
#include "sort.h"
#include "omplock.h"
#include "outputfiles.h"

void GetFileNames(const string &SpecFileName, vector<string> &FileNames);
void ReadPDBs(const vector<string> &FileNames, vector<PDBChain *> &Structures);

uint PpcSearch1(PpcAligner &PA, PDBChain &Q, vector<PDBChain*> &RefPDBs)
	{
	const uint RefN = SIZE(RefPDBs);
	PA.SetQuery(Q);
	vector<double> RMSDs;
	vector<string> RefLabels;
	for (uint iR = 0; iR < RefN; ++iR)
		{
		const PDBChain &R = *RefPDBs[iR];

		if (opt_self && Q.m_ChainLabel == R.m_ChainLabel)
			continue;

		PA.SetRef(R);
		double MotifRMSD = PA.GetMotifRMSD();
		if (MotifRMSD <= MaxMotifRMSD)
			{
			RMSDs.push_back(MotifRMSD);
			string RefLabel;
			R.GetLabel(RefLabel);
			RefLabels.push_back(RefLabel);
			}
		}
	uint HitCount = SIZE(RMSDs);
	if (HitCount == 0)
		return 0;
	if (g_ftsv == 0)
		return HitCount;

	vector<uint> Order(HitCount);
	QuickSortOrder(RMSDs.data(), HitCount, Order.data());

	string QLabel;
	Q.GetLabel(QLabel);

	LockOutput();
	for (uint k = 0; k < HitCount; ++k)
		{
		uint i = Order[k];
		fprintf(g_ftsv, "%s", QLabel.c_str());
		fprintf(g_ftsv, "\t%s", RefLabels[i].c_str());
		fprintf(g_ftsv, "\t%.2f\n", RMSDs[i]);
		}
	UnlockOutput();
	return HitCount;
	}

void ReadPpc(const string &FN, vector<PDBChain *> &Chains)
	{
	Chains.clear();
	CalReader CR;
	CR.Open(FN);
	for (;;)
		{
		PDBChain *Chain = new PDBChain;
		bool Ok = CR.GetNext(*Chain);
		if (!Ok)
			return;
		Chains.push_back(Chain);
		}
	}

void cmd_search3d_ppc()
	{
	const string &QueryFN = opt_search3d_ppc;
	const string &RefFN = opt_ref;

	CalReader CR;
	CR.Open(QueryFN);

	vector<PDBChain *> RefPDBs;
	ReadPpc(RefFN, RefPDBs);

	const uint RefN = SIZE(RefPDBs);

	uint ThreadCount = GetRequestedThreadCount();
	vector<vector<PDBChain *> > QVecs(ThreadCount);
	vector<PpcAligner> PAs(ThreadCount);

	uint ThreadFinishedCount = 0;
	uint DoneCount = 0;
	uint HitCount = 0;
	vector<PDBChain> Qs(ThreadCount);
	vector<bool> ThreadDone(ThreadCount);

#pragma omp parallel num_threads(ThreadCount)
	for (;;)
		{
		uint ThreadIndex = GetThreadIndex();
		PDBChain &Q = Qs[ThreadIndex];
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
		PpcAligner &PA = PAs[ThreadIndex];

		uint n = PpcSearch1(PA, Q, RefPDBs);
		if (n > 0)
			++HitCount;
		}
	Progress("100.0% done, %u / %u hits\n", DoneCount, HitCount);
	}

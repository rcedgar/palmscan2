#include "myutils.h"
#include "pssm.h"
#include "rdrpsearcher.h"
#include "seqinfo.h"
#include "fastaseqsource.h"
#include "omplock.h"
#include <time.h>

static uint g_QueryCount;
static uint g_FoundCount;

vector<string> g_ExcludeNames;

void SetExcludes()
	{
	if (!optset_exclude)
		return;

	string NamesStr = string(opt_exclude);
	Split(NamesStr, g_ExcludeNames, '+');
	const uint N = SIZE(g_ExcludeNames);
	for (uint i = 0; i < N; ++i)
		ProgressLog("  Exclude %s\n", g_ExcludeNames[i].c_str());
	}

void SeqToUpper(string &Seq)
	{
	const uint QL = SIZE(Seq);
	for (uint i = 0; i < QL; ++i)
		Seq[i] = toupper(Seq[i]);
	}

void cmd_search_pssms()
	{
	string QueryFileName = opt_search_pssms;

	asserta(optset_model);
	const string &ModelFileName = opt_model;
	bool Trace = opt_trace;

	if (!opt_notrunclabels)
		opt_trunclabels = true;

	SetExcludes();

	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

	const uint ThreadCount = GetRequestedThreadCount();
	vector<ObjMgr *> OMs;
	vector<RdRpSearcher *> Mods;
	for (int i = 0; i < int(ThreadCount); ++i)
		{
		ObjMgr *OM = new ObjMgr;
		OMs.push_back(OM);
		RdRpSearcher *RS = new RdRpSearcher;
		RS->Init(Model);
		Mods.push_back(RS);
		}

	FASTASeqSource *SS = new FASTASeqSource;
	SS->Open(QueryFileName);

	RdRpSearcher::InitOutput();
	uint LastElapsedSecs = 0;
	uint CurrElapsedSecs = 0;
	ProgressStep(0, 1000, "Searching");
#pragma omp parallel num_threads(ThreadCount)
	{
	int ThreadIndex = GetThreadIndex();
	RdRpSearcher &RS = *Mods[ThreadIndex];
	RS.m_Trace = Trace;
	ObjMgr *OM = OMs[ThreadIndex];
	for (;;)
		{
		SeqInfo *QSI = OM->GetSeqInfo();

		if (ThreadIndex == 0)
			{
			if (g_QueryCount%100 == 0)
				{
				CurrElapsedSecs = GetElapsedSecs();
				if (CurrElapsedSecs > LastElapsedSecs)
					{
					uint Pct10 = SS->GetPctDoneX10();
					double HitPct = GetPct(g_FoundCount, g_QueryCount);
					ProgressStep(Pct10, 1000, "Searching %u / %u hits (%.1f%%)",
					  g_FoundCount, g_QueryCount, HitPct);
					LastElapsedSecs = CurrElapsedSecs;
					}
				}
			}

		bool Ok = SS->GetNext(QSI);
		if (!Ok)
			break;

		const string Label = string(QSI->m_Label);
		string Seq;
		for (uint i = 0; i < QSI->m_L; ++i)
			Seq += char(toupper(QSI->m_Seq[i]));

		RS.Search(Label, Seq);
		RS.WriteOutput();
		OM->Down(QSI);

		Lock();
		++g_QueryCount;
		if (RS.m_TopPalmHit.m_Score > 0)
			++g_FoundCount;
		Unlock();
		}
	}

	double HitPct = GetPct(g_FoundCount, g_QueryCount);
	ProgressStep(999, 1000, "Searching %u/%u hits (%.1f%%)",
	  g_FoundCount, g_QueryCount, HitPct);
	}

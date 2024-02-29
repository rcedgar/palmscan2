#include "myutils.h"
#include "tma.h"
#include "outputfiles.h"
#include <map>

static FILE *g_fout;

void Shuffle(vector<unsigned> &v);
void TMAlignPair(TMA &T, const PDBChain &Q, const PDBChain &R);

static const uint MAXFAM = 32;
static const double MINTM = 0.5;

static void DoFam(TMA &T, const string &Fam,
  vector<PDBChain *> Chains, vector<uint> &ChainIndexes)
	{
	uint N = SIZE(ChainIndexes);
#pragma omp critical
	{
	Log("[%4u] %s\n", N, Fam.c_str());
	}
	if (N == 0)
		{
		Warning("Empty family %s", Fam.c_str());
		return;
		}
	if (N == 1)
		return;
	if (N > MAXFAM)
		{
		Shuffle(ChainIndexes);
		ChainIndexes.resize(MAXFAM);
		N = MAXFAM;
		}
	for (uint i = 0; i < N; ++i)
		{
		uint ChainIndexi = ChainIndexes[i];
		const PDBChain &Q = *Chains[ChainIndexi];
		for (uint j = i + 1; j < N; ++j)
			{
			uint ChainIndexj = ChainIndexes[j];
			const PDBChain &R = *Chains[ChainIndexj];
			double TM = T.AlignChains(Q, R);
			if (TM < MINTM)
				continue;
#pragma omp critical
			{
			uint Ids = 0;
			uint Cols = 0;
			const string &QRow = T.m_QRow;
			const string &RRow = T.m_RRow;
			const uint ColCount = SIZE(QRow);
			asserta(SIZE(RRow) == ColCount);
			for (uint Col = 0; Col < ColCount; ++Col)
				{
				char q = QRow[Col];
				char r = RRow[Col];
				if (q == '-' || r == '-')
					continue;
				++Cols;
				if (q == r)
					++Ids;
				}
			double PctId = GetPct(Ids, Cols);
			if (g_ftsv)
				fprintf(g_ftsv, "%s\t%s\t%.4f\t%.1f\n",
				  Q.m_Label.c_str(), R.m_Label.c_str(), TM, PctId);
			if (g_fout)
				{
				fprintf(g_fout, ">%s/%.4f/%.1f\n", Q.m_Label.c_str(), TM, PctId);
				fprintf(g_fout, "%s\n", QRow.c_str());
				fprintf(g_fout, ">%s\n", R.m_Label.c_str());
				fprintf(g_fout, "%s\n", T.m_RRow.c_str());
				}
			}
			}
		}
	}

void cmd_tm_scop()
	{
	const string &InputFileName = opt_tm_scop;
	g_fout = CreateStdioFile(opt_output);
	vector<PDBChain *> Chains;
	ReadChains(InputFileName, Chains);
	const uint ChainCount = SIZE(Chains);
	map<string, vector<uint> > FamToChainIndexes;
	vector<string> Fams;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const string &Label = Chains[ChainIndex]->m_Label;
		vector<string> Fields;
		Split(Label, Fields, '/');
		asserta(SIZE(Fields) == 2);
		const string &Fam = Fields[1];
		map<string, vector<uint> >::const_iterator p = FamToChainIndexes.find(Fam);
		if (p == FamToChainIndexes.end())
			{
			Fams.push_back(Fam);
			vector<uint> Empty;
			FamToChainIndexes[Fam] = Empty;
			}
		FamToChainIndexes[Fam].push_back(ChainIndex);
		}
	const uint FamCount = SIZE(Fams);
	uint ThreadCount = GetRequestedThreadCount();

	vector<TMA *> Ts;
	for (uint i = 0; i < ThreadCount; ++i)
		{
		TMA *T = new TMA;
		Ts.push_back(T);
		}

	uint FamIndex = 0;
#pragma omp parallel num_threads(ThreadCount)
	for (;;)
		{
		uint MyFamIndex = UINT_MAX;
#pragma omp critical
		{
		if (FamIndex < FamCount)
			{
			MyFamIndex = FamIndex;
			ProgressStep(FamIndex, (uint) FamCount, "Aligning");
			++FamIndex;
			}
		}
		if (MyFamIndex == UINT_MAX)
			break;

		const string &Fam = Fams[MyFamIndex];
		vector<uint> &ChainIndexes = FamToChainIndexes[Fam];
		uint ThreadIndex = GetThreadIndex();
		TMA &T = *Ts[ThreadIndex];
		DoFam(T, Fam, Chains, ChainIndexes);
		}
	CloseStdioFile(g_fout);
	}

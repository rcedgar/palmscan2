#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "quarts.h"
#include "binner.h"

void cmd_distmx()
	{
	const string &QFN = opt_distmx;
	asserta(g_ftsv != 0);

	vector<PDBChain *> Qs;
	ReadChains(QFN, Qs);

	const uint N = SIZE(Qs);
	
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Q = *Qs[i];
		uint QL = Q.GetSeqLength();
		const char *LabelQ = Q.m_Label.c_str();
		const string &SeqQ = Q.m_Seq;
		for (uint i = 0; i < QL; ++i)
			{
			char ci = SeqQ[i];
			for (uint j = i+1; j < QL; ++j)
				{
				char cj = SeqQ[j];
				double d = Q.GetDist(i, j);
				fprintf(g_ftsv, "%s", LabelQ);
				fprintf(g_ftsv, "\t%u", i);
				fprintf(g_ftsv, "\t%u", j);
				fprintf(g_ftsv, "\t%c", ci);
				fprintf(g_ftsv, "\t%c", cj);
				fprintf(g_ftsv, "\t%.1f", d);
				fprintf(g_ftsv, "\n");
				}
			}
		}
	}

void cmd_backbone_dist()
	{
	const string &QFN = opt_backbone_dist;

	vector<PDBChain *> Qs;
	ReadChains(QFN, Qs);

	const uint N = SIZE(Qs);

	vector<double> v;
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Q = *Qs[i];
		uint QL = Q.GetSeqLength();
		for (uint i = 0; i + 1 < QL; ++i)
			{
			double d = Q.GetDist(i, i+1);
			if (d < 2 || d > 4)
				continue;
			v.push_back(d);
			}
		}
	QuartsDouble QD;
	GetQuartsDouble(v, QD);
	QD.LogMe();
	}

double GetNENDistance(const PDBChain &Q, uint Pos)
	{
	double NEN = 999;
	uint QL = Q.GetSeqLength();
	for (uint i = 0; i < QL; ++i)
		{
		if (i + 5 >= Pos || Pos + 5 <= i)
			continue;
		double d = Q.GetDist(i, Pos);
		NEN = min(d, NEN);
		}
	return NEN;
	}

//NEN: Min=2.61, LoQ=5.49, Med=7.72, HiQ=9.84, Max=20, Avg=8.11, StdDev=3.1
// NN: Min=2.37, LoQ=3.79, Med=3.81, HiQ=3.82, Max=4, Avg=3.8, StdDev=0.0492
void cmd_nen_dist()
	{
	const string &QFN = opt_nen_dist;

	vector<PDBChain *> Qs;
	ReadChains(QFN, Qs);

	const uint N = SIZE(Qs);

	vector<double> vNEN;
	vector<double> vNN;
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Q = *Qs[i];
		uint QL = Q.GetSeqLength();
		for (uint i = 0; i < QL; ++i)
			{
			double d = GetNENDistance(Q, i);
			if (d < 20)
				vNEN.push_back(d);

			if (i > 0)
				{
				d = Q.GetDist(i, i-1);
				if (d > 2 && d < 4)
					vNN.push_back(d);
				}
			}
		}
	QuartsDouble QD;
	GetQuartsDouble(vNEN, QD);
	Log("NEN: ");
	QD.LogMe();

	GetQuartsDouble(vNN, QD);
	Log(" NN: ");
	QD.LogMe();

	Binner<double> B(vNEN, 32);
	B.ToTsv("nen_dist.tsv");
	}

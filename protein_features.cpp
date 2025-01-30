#include "myutils.h"
#include "pdbchain.h"
#include "quarts.h"
#include "binner.h"

double GetNENDistance(const PDBChain &Q, uint Pos);

static void GetCenter(const PDBChain &Q,
					  double &x, double &y, double &z)
	{
	x = 0;
	y = 0;
	z = 0;
	const uint L = Q.GetSeqLength();
	if (L == 0)
		return;
	for (uint i = 0; i < L; ++i)
		{
		x += Q.m_Xs[i];
		y += Q.m_Ys[i];
		z += Q.m_Zs[i];
		}
	x /= L;
	y /= L;
	z /= L;
	}

static double GetMDL(const PDBChain &Q)
	{
	double cx, cy, cz;
	GetCenter(Q, cx, cy, cz);

	double MaxDist = 0;
	const uint L = Q.GetSeqLength();
	if (L == 0)
		return 0;
	for (uint i = 0; i < L; ++i)
		{
		for (uint j = i+1; j < L; ++j)
			{
			double d = Q.GetDist(i, j);
			MaxDist = max(MaxDist, d);
			}
		}
	return MaxDist/L;
	}

static vector<double> s_NENMeds;
static vector<double> s_NENStdDevs;
static vector<double> s_MDLs;

static void ProteinFeatures(FILE *f, PDBChain &Q)
	{
	if (f == 0)
		return;

	double MDL = GetMDL(Q);

	vector<double> NENs;
	uint QL = Q.GetSeqLength();
	for (uint i = 0; i < QL; ++i)
		{
		double d = GetNENDistance(Q, i);
		if (d < 20)
			NENs.push_back(d);
		}
	QuartsDouble QD;
	GetQuartsDouble(NENs, QD);

	double NENMed = QD.Med;
	double NENStdDev = QD.StdDev;

	s_NENMeds.push_back(NENMed);
	s_NENStdDevs.push_back(NENStdDev);
	s_MDLs.push_back(MDL);

	fprintf(f, "%s", Q.m_Label.c_str());
	fprintf(f, "\t%.3g", NENMed);
	fprintf(f, "\t%.3g", NENStdDev);
	fprintf(f, "\t%.3g", MDL);
	fprintf(f, "\n");
	}

void cmd_protein_features()
	{
	const string &InputFN = g_Arg1;
	FILE *f = CreateStdioFile(opt_output);
	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);
	const uint N = SIZE(Chains);
	for (uint i = 0; i < N; ++i)
		ProteinFeatures(f, *Chains[i]);
	CloseStdioFile(f);

	Binner<double> B_NENMed(s_NENMeds, 16);
	Binner<double> B_NENStdDev(s_NENStdDevs, 16);
	Binner<double> B_MDL(s_MDLs, 16);

	B_NENMed.ToTsv("nenmed_bins.tsv");
	B_NENStdDev.ToTsv("nenstddev_bins.tsv");
	B_MDL.ToTsv("mdl_bins.tsv");
	}

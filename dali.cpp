#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include <map>

// Holm2019 eq(A3)
static double Getm(double L)
	{
	if (L > 400)
		return Getm(400) + L - 400;

	double y = 7.95;
	y += 0.71*L;
	y += 0.000259*L*L;
	y += 0.00000192*L*L*L;
	return y;
	}

/***
DaliLite v5
===========================
comparemodules.f, line 1436
===========================
        enveloperadius=20.0
        x=1/(enveloperadius*enveloperadius)
        do i=0,100
                wght(i)=exp(-x*i*i)
        end do
***/
static double Weight(double y)
	{
	const double D = 20.0;
	const double x = 1.0/(D*D);
	double w = exp(-x*y*y);
	return w;
	}

/***
DaliLite v5
===========================
comparemodules.f, line 1397
  a, b are integer distances in units of 1/10 Angstrom,
  so multiply by 10 to get Angstroms.
===========================
        function dpscorefun(a,b) result(s)
        implicit none
        include 'parsizes.for'
        real s
        integer*2 a,b
c
        real x,y,d0
        logical lela
        parameter(lela=.true.)
        parameter(d0=0.20)
c !!!   elastic uses weights !!!
        x=float(abs(a-b))/10
        if(lela) then
                y=float(a+b)/20
                if(y.gt.100) then
                        s=0.0
                else
                        if(y.gt.0) then
                          s=wght(nint(y))*(d0-x/y)
                        else
                          s=wght(nint(y))*d0
                        end if
                end if
        end if
***/
static double dpscorefun(double a, double b)
	{
	double Score = 0;
	const double d0 = 0.20;
	double x = fabs(a - b);
	double y = (a + b)/2;
	if (y > 100)
		Score = 0;
	else
		{
		if (y > 0)
			Score = Weight(y)*(d0 - x/y);
		else
			Score = Weight(y)*d0;
		}
	return Score;
	}

// Holm2019 eq(A1)
double GetDALI(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs)
	{
	const uint Lali = SIZE(PosQs);
	asserta(SIZE(PosTs) == Lali);

	const double Theta = 0.2;
	const double D = 20;
	double Sum = 0;
	for (uint i = 0; i < Lali; ++i)
		{
		uint PosQi = PosQs[i];
		uint PosTi = PosTs[i];
		for (uint j = 0; j < Lali; ++j)
			{
			if (i == j)
				Sum += Theta;
			else
				{
				uint PosQj = PosQs[j];
				uint PosTj = PosTs[j];

				double dij_Q = Q.GetDist(PosQi, PosQj);
				double dij_T = T.GetDist(PosTi, PosTj);
				double x = dpscorefun(dij_Q, dij_T);
				Sum += x;
				}
			}
		}
	return Sum;
	}

/***
"./src/comparemodules.f" line 1473
        function zscore(l1,l2,score) result(z)
        implicit none
        real z,score
        integer l1,l2
c
        real n12,mean,sigma,x
c
        n12=sqrt(float(l1*l2))
        x=min(n12,400.0)
        mean=7.9494+0.70852*x+2.5895e-4*x*x-1.9156e-6*x*x*x
        if(n12.gt.400.0) mean=mean+(n12-400.0)*1.0              ! hack !
        sigma=0.50*mean
        z=(score-mean)/max(1.0,sigma)

        return
        end function zscore

***/

static double GetDALIZ(double DALI, uint QL, uint TL)
	{
	double n12 = sqrt(QL*TL);
	double x = min(n12, 400.0);
	double mean = 7.9494 + 0.70852*x + 2.5895e-4*x*x - 1.9156e-6*x*x*x;
	if (n12 > 400)
		mean += n12 - 400.0;
	double sigma = 0.5*mean;
	double z = (DALI - mean)/max(1.0, sigma);
	return z;
	}

#if 0
// Per Holm2019, source code is different
double GetDALIZ(double DALI, uint QL, uint TL)
	{
	double L = sqrt(QL*TL);
	double m = Getm(L);
	double sigma = 0.5*m;
	double Z = (DALI - m)/sigma;
	return Z;
	}
#endif

static void GetPosVecs(const string &QRow, const string &TRow,
  vector<uint> &PosQs, vector<uint> &PosTs)
	{
	PosQs.clear();
	PosTs.clear();
	uint ColCount = SIZE(QRow);
	asserta(SIZE(TRow) == ColCount);
	uint PosQ = 0;
	uint PosT = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char q = QRow[Col];
		char t = TRow[Col];

		if (q != '-' && t != '-' &&
		  isupper(q) && isupper(t))
			{
			PosQs.push_back(PosQ);
			PosTs.push_back(PosT);
			}

		if (q != '-')
			++PosQ;
		if (t != '-')
			++PosT;
		}
	}

void cmd_dali()
	{
	const string &QFN = opt_dali;
	const string &TsvFN = opt_ref;

	vector<PDBChain *> Chains;
	ReadChains(QFN, Chains);

	const uint NQ = SIZE(Chains);
	map<string, uint> LabelToIndex;
	for (uint i = 0; i < NQ; ++i)
		{
		const string &Label = Chains[i]->m_Label;
		LabelToIndex[Label] = i;
		}

	vector<string> Fields;
	Split(opt_fieldnrs, Fields, ',');
	assert(SIZE(Fields) == 4);

	uint Qfn = StrToUint(Fields[0]);
	uint Tfn = StrToUint(Fields[1]);
	uint QRowfn = StrToUint(Fields[2]);
	uint TRowfn = StrToUint(Fields[3]);

	FILE *f = OpenStdioFile(TsvFN);
	string Line;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');

		const string &Query = Fields[Qfn];
		const string &Target = Fields[Tfn];
		map<string, uint>::const_iterator pQ = LabelToIndex.find(Query);
		map<string, uint>::const_iterator pT = LabelToIndex.find(Target);
		if (pQ == LabelToIndex.end() || pT == LabelToIndex.end())
			continue;

		uint iQ = pQ->second;
		uint iT = pT->second;

		const PDBChain &Q = *Chains[iQ];
		const PDBChain &T = *Chains[iT];
		asserta(Q.m_Label == Query && T.m_Label == Target);

		const uint QL = Q.GetSeqLength();
		const uint TL = T.GetSeqLength();

		const string &QRow = Fields[QRowfn];
		const string &TRow = Fields[TRowfn];

		vector<uint> PosQs;
		vector<uint> PosTs;
		GetPosVecs(QRow, TRow, PosQs, PosTs);
		double DALI = GetDALI(Q, T, PosQs, PosTs);
		double Z = GetDALIZ(DALI, QL, TL);
		if (g_ftsv)
			{
			fprintf(g_ftsv, "%s", Q.m_Label.c_str());
			fprintf(g_ftsv, "\t%s", T.m_Label.c_str());
			fprintf(g_ftsv, "\t%.3g", Z);
			fprintf(g_ftsv, "\n");
			}
		}
	CloseStdioFile(f);
	}

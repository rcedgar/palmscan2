#include "myutils.h"
#include "cmp.h"
#include <map>

static double g_Gate_Rdrp_Probs[256];
static map<string, double> g_GDD_To_Rdrp_Prob;

static void Init_Gate_RdRp(int c, double P)
	{
	asserta(c >  0);
	g_Gate_Rdrp_Probs[toupper(c)] = P;
	g_Gate_Rdrp_Probs[tolower(c)] = P;
	}

static void Init_GDD_RdRp(const string &GDD, double P)
	{
	asserta(SIZE(GDD) == 3);
	g_GDD_To_Rdrp_Prob[GDD] = P;
	}

static bool Init()
	{
	for (uint i = 0; i < 256; ++i)
		g_Gate_Rdrp_Probs[i] = 0.5;

#define GATE_RDRP(c, P)	Init_Gate_RdRp(c, P);
GATE_RDRP('D', 0.9790) // gate  D       3068    24      0.6743  0.0144  P(RdRp) = 0.9790
GATE_RDRP('N', 0.6124) // gate  N       594     162     0.1314  0.0832  P(RdRp) = 0.6124
GATE_RDRP('S', 0.8850) // gate  S       257     10      0.0575  0.0075  P(RdRp) = 0.8850
GATE_RDRP('G', 0.9020) // gate  G       204     5       0.0459  0.0050  P(RdRp) = 0.9020
GATE_RDRP('T', 0.9193) // gate  T       202     3       0.0454  0.0040  P(RdRp) = 0.9193
GATE_RDRP('C', 0.8583) // gate  C       105     3       0.0241  0.0040  P(RdRp) = 0.8583
GATE_RDRP('E', 0.7296) // gate  E       44      3       0.0108  0.0040  P(RdRp) = 0.7296
GATE_RDRP('F', 0.0043) // gate  F       1       609     0.0013  0.3059  P( RT ) = 0.9957
GATE_RDRP('Y', 0.0059) // gate  Y       0       364     0.0011  0.1839  P( RT ) = 0.9941
#undef GATE_RDRP

#define GDD(c, P)	Init_GDD_RdRp(c, P);
GDD("GDD", 0.9826) // gdd       GDD     3587    23      0.7882  0.0140  P(RdRp) = 0.9826
GDD("SDD", 0.9440) // gdd       SDD     531     9       0.1176  0.0070  P(RdRp) = 0.9440
GDD("GDN", 0.9626) // gdd       GDN     346     1       0.0770  0.0030  P(RdRp) = 0.9626
GDD("ADD", 0.0047) // gdd       ADD     4       834     0.0020  0.4180  P( RT ) = 0.9953
GDD("VDD", 0.0107) // gdd       VDD     0       199     0.0011  0.1016  P( RT ) = 0.9893
GDD("IDD", 0.0422) // gdd       IDD     0       45      0.0011  0.0249  P( RT ) = 0.9578
GDD("MDD", 0.1039) // gdd       MDD     0       14      0.0011  0.0095  P( RT ) = 0.8961
GDD("TDD", 0.1221) // gdd       TDD     1       14      0.0013  0.0095  P( RT ) = 0.8779
#undef GDD_RDRP

	return true;
	}

static bool g_InitDone = Init();

double CMP::GetRdRpProb_Gate(int Gate)
	{
	double P = g_Gate_Rdrp_Probs[Gate];
	return P;
	}

double CMP::GetRdRpProb_GDD(const string &GDD)
	{
	map<string, double>::const_iterator p = g_GDD_To_Rdrp_Prob.find(GDD);
	if (p == g_GDD_To_Rdrp_Prob.end())
		return 0.1;
	double P = p->second;
	return P;
	}

double CMP::GetRdRpProb(int Gate, const string &GDD)
	{
	double P_gate = GetRdRpProb_Gate(Gate);
	double P_gdd = GetRdRpProb_GDD(GDD);
	double P = min(P_gate, P_gdd);
	return P;
	}
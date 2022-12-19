#include "myutils.h"
#include "ppsp.h"

static double g_Gate_Rdrp_Probs[256];
static double g_CmfX_Rdrp_Probs[256];

static void Init_Gate_RdRp(int c, double P)
	{
	asserta(c >  0);
	g_Gate_Rdrp_Probs[toupper(c)] = P;
	g_Gate_Rdrp_Probs[tolower(c)] = P;
	}

static void Init_CmfX_RdRp(int c, double P)
	{
	asserta(c > 0);
	g_CmfX_Rdrp_Probs[toupper(c)] = P;
	g_CmfX_Rdrp_Probs[tolower(c)] = P;
	}

static bool Init()
	{
	for (uint i = 0; i < 256; ++i)
		{
		g_Gate_Rdrp_Probs[i] = 0.5;
		g_CmfX_Rdrp_Probs[i] = 0.5;
		}

#define GATE_RDRP(c, P)	Init_Gate_RdRp(c, P);
GATE_RDRP('D', 0.9790) // gate  D       3068    24      0.6743  0.0144  P(RdRp) = 0.9790
GATE_RDRP('N', 0.6124) // gate  N       594     162     0.1314  0.0832  P(RdRp) = 0.6124
GATE_RDRP('F', 0.0043) // gate  F       1       609     0.0013  0.3059  P( RT ) = 0.9957
GATE_RDRP('Y', 0.0059) // gate  Y       0       364     0.0011  0.1839  P( RT ) = 0.9941
GATE_RDRP('S', 0.8850) // gate  S       257     10      0.0575  0.0075  P(RdRp) = 0.8850
GATE_RDRP('G', 0.9020) // gate  G       204     5       0.0459  0.0050  P(RdRp) = 0.9020
GATE_RDRP('T', 0.9193) // gate  T       202     3       0.0454  0.0040  P(RdRp) = 0.9193
GATE_RDRP('C', 0.8583) // gate  C       105     3       0.0241  0.0040  P(RdRp) = 0.8583
GATE_RDRP('E', 0.7296) // gate  E       44      3       0.0108  0.0040  P(RdRp) = 0.7296
#undef GATE_RDRP

#define CMFX_RDRP(c, P)	Init_CmfX_RdRp(c, P);
CMFX_RDRP('G', 0.9704) // C-X   G       3939    48      0.8655  0.0264  P(RdRp) = 0.9704
CMFX_RDRP('A', 0.0057) // C-X   A       6       840     0.0024  0.4210  P( RT ) = 0.9943
CMFX_RDRP('S', 0.9292) // C-X   S       531     13      0.1176  0.0090  P(RdRp) = 0.9292
CMFX_RDRP('V', 0.0107) // C-X   V       0       199     0.0011  0.1016  P( RT ) = 0.9893
CMFX_RDRP('I', 0.0422) // C-X   I       0       45      0.0011  0.0249  P( RT ) = 0.9578
CMFX_RDRP('M', 0.0992) // C-X   M       0       15      0.0011  0.0100  P( RT ) = 0.9008
CMFX_RDRP('T', 0.1221) // C-X   T       1       14      0.0013  0.0095  P( RT ) = 0.8779
#undef CMFX_RDRP

	return true;
	}

static bool g_InitDone = Init();

double PPSP::GetRdRpProb_Gate(int Gate)
	{
	double P = g_Gate_Rdrp_Probs[Gate];
	return P;
	}

double PPSP::GetRdRpProb_CmfX(int X)
	{
	double P = g_CmfX_Rdrp_Probs[X];
	return P;
	}

double PPSP::GetRdRpProb(int Gate, int CmfX)
	{
	double P_gate = GetRdRpProb_Gate(Gate);
	double P_cmfx = GetRdRpProb_CmfX(CmfX);
	double P = min(P_gate, P_cmfx);
	return P;
	}

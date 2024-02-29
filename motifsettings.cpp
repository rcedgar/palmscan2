#include "myutils.h"
#include "motifsettings.h"

#if MOTIF_SETTINGS

uint g_NA;
uint g_NB;
uint g_NC;

uint g_GoodNA;
uint g_GoodNB;
uint g_GoodNC;

uint g_PSSMNA;
uint g_PSSMNB;
uint g_PSSMNC;

uint g_PSSMGoodNA;
uint g_PSSMGoodNB;
uint g_PSSMGoodNC;

void ClearMotifGoodCounts()
	{
	g_NA = 0;
	g_NB = 0;
	g_NC = 0;
	g_GoodNA = 0;
	g_GoodNB = 0;
	g_GoodNC = 0;
	g_PSSMNA = 0;
	g_PSSMNB = 0;
	g_PSSMNC = 0;
	g_PSSMGoodNA = 0;
	g_PSSMGoodNB = 0;
	g_PSSMGoodNC = 0;
	}

void MotifSettingsToFile(FILE *f)
	{
	if (f == 0)
		return;
	fprintf(f, "MotifSettings,");
	fprintf(f, "A,%u,%u,", g_LA, g_OffAd);
	fprintf(f, "B,%u,%u,", g_LB, g_OffBg);
	fprintf(f, "C,%u,%u\n", g_LC, g_OffCd);
	}

void MotifSettingsFromLine(const string &Line)
	{
	vector<string> Fields;
	Split(Line, Fields, ',');
	asserta(SIZE(Fields) == 10);
	asserta(Fields[0] == "MotifSettings");
	asserta(Fields[1] == "A");
	asserta(Fields[4] == "B");
	asserta(Fields[7] == "C");
	asserta(StrToUint(Fields[2]) == g_LA);
	asserta(StrToUint(Fields[3]) == g_OffAd);
	asserta(StrToUint(Fields[5]) == g_LB);
	asserta(StrToUint(Fields[6]) == g_OffBg);
	asserta(StrToUint(Fields[8]) == g_LC);
	asserta(StrToUint(Fields[9]) == g_OffCd);
	}

void WriteMotifSettings(FILE *f)
	{
	if (f == 0)
		return;

	fprintf(f, "Motif settings\n");
	fprintf(f, "A length %u  offset D %u  -left %u  +right %u\n",
	  g_LA, g_OffAd, g_LeftSubA, g_RightAddA);

	fprintf(f, "B length %u  offset G %u  -left %u  +right %u\n",
	  g_LB, g_OffBg, g_LeftSubB, g_RightAddB);

	fprintf(f, "C length %u  offset D %u  -left %u  +right %u\n",
	  g_LC, g_OffCd, g_LeftSubC, g_RightAddC);

	if (g_NA > 0)
		fprintf(f, "  GoodA %u/%u(%.1f%%)\n",
		  g_GoodNA, g_NA, GetPct(g_GoodNA, g_NA));

	if (g_NB > 0)
		fprintf(f, "  GoodB %u/%u(%.1f%%)\n",
		  g_GoodNB, g_NB, GetPct(g_GoodNB, g_NB));

	if (g_NC > 0)
		fprintf(f, "  GoodC %u/%u(%.1f%%)\n",
		  g_GoodNC, g_NC, GetPct(g_GoodNC, g_NC));

	if (g_PSSMNA > 0)
		fprintf(f, "  PSSMGoodA %u/%u(%.1f%%)\n",
		  g_PSSMGoodNA, g_PSSMNA, GetPct(g_PSSMGoodNA, g_NA));

	if (g_PSSMGoodNB > 0)
		fprintf(f, "  PSSMGoodB %u/%u(%.1f%%)\n",
		  g_PSSMGoodNB, g_PSSMNB, GetPct(g_PSSMGoodNB, g_NB));

	if (g_PSSMGoodNC > 0)
		fprintf(f, "  PSSMGoodC %u/%u(%.1f%%)\n",
		  g_PSSMGoodNC, g_PSSMNC, GetPct(g_PSSMGoodNC, g_NC));
	}

void CheckA(const string &A)
	{
	if (A == "" || A == ".")
		return;
	++g_NA;
	asserta(SIZE(A) == g_LA);
	if (toupper(A[g_OffAd]) == 'D')
		++g_GoodNA;
	}

void CheckB(const string &B)
	{
	if (B == "" || B == ".")
		return;
	++g_NB;
	asserta(SIZE(B) == g_LB);
	if (toupper(B[g_OffBg]) == 'G')
		++g_GoodNB;
	}

void CheckC(const string &C)
	{
	if (C == "" || C == ".")
		return;
	++g_NC;
	asserta(SIZE(C) == g_LC);
	if (toupper(C[g_OffCd]) == 'D')
		++g_GoodNC;
	}

void CheckABC3(const string &A, const string &B, const string &C)
	{
	CheckA(A);
	CheckB(B);
	CheckC(C);
	}

void CheckABC(const vector<string> &ABC)
	{
	asserta(SIZE(ABC) == 3);
	CheckA(ABC[0]);
	CheckB(ABC[1]);
	CheckC(ABC[2]);
	}

void PSSMCheckA(const string &A)
	{
	if (A == "" || A == ".")
		return;
	asserta(SIZE(A) == 12);
	++g_PSSMNA;
	}

void PSSMCheckB(const string &B)
	{
	if (B == "" || B == ".")
		return;
	asserta(SIZE(B) == 14);
	++g_PSSMNB;
	}

void PSSMCheckC(const string &C)
	{
	if (C == "" || C == ".")
		return;
	asserta(SIZE(C) == 8);
	++g_PSSMNC;
	}

void PSSMCheckABC(const string &A, const string &B, const string &C)
	{
	PSSMCheckA(A);
	PSSMCheckB(A);
	PSSMCheckC(A);
	}

#endif // MOTIF_SETTINGS
#pragma once

#define	MOTIF_SETTINGS	0

const uint g_LeftSubA = 2;
const uint g_LeftSubB = 0;
const uint g_LeftSubC = 2;
const uint g_RightAddA = 4;
const uint g_RightAddB = 8;
const uint g_RightAddC = 4;

const uint g_LA = g_LeftSubA + 8 + g_RightAddA;
const uint g_LB = g_LeftSubB + 12 + g_RightAddB;
const uint g_LC = g_LeftSubC + 8 + g_RightAddC;

const uint g_OffAd = 3 + g_LeftSubA;
const uint g_OffBg = 1 + g_LeftSubB;
const uint g_OffCd = 3 + g_LeftSubC;

#if MOTIF_SETTINGS

extern uint g_NA;
extern uint g_NB;
extern uint g_NC;

extern uint g_GoodNA;
extern uint g_GoodNB;
extern uint g_GoodNC;

extern uint g_PSSMNA;
extern uint g_PSSMNB;
extern uint g_PSSMNC;

extern uint g_PSSMGoodNA;
extern uint g_PSSMGoodNB;
extern uint g_PSSMGoodNC;

void WriteMotifSettings(FILE *f);
void MotifSettingsToFile(FILE *f);

void CheckA(const string &A);
void CheckB(const string &B);
void CheckC(const string &C);
void CheckABC3(const string &A, const string &B, const string &C);
void CheckABC(const vector<string> &ABC);

void PSSMCheckA(const string &A);
void PSSMCheckB(const string &B);
void PSSMCheckC(const string &C);
void PSSMCheckABC(const string &A, const string &B, const string &C);

void MotifSettingsFromLine(const string &Line);
void ClearMotifGoodCounts();
#else

#pragma warning(disable: 4390) // empty controlled statement

#define WriteMotifSettings(f) /* empty */

#define MotifSettingsToFile(f) fprintf(f, "MotifSettings,-\n")

#define CheckA(A) /* empty */
#define CheckB(B) /* empty */
#define CheckC(C) /* empty */
#define CheckABC3(A, B, C) /* empty */
#define CheckABC(ABC) /* empty */

#define PSSMCheckA(A) /* empty */
#define PSSMCheckB(B) /* empty */
#define PSSMCheckC(C) /* empty */
#define PSSMCheckABC(A, B, C) /* empty */

#define MotifSettingsFromLine(Line) /* empty */
#define ClearMotifGoodCounts() /* empty */

#endif // 0

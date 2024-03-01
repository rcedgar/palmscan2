#include "myutils.h"
#include "alpha.h"
#include "xprof.h"

static double g_3di_Medians[20][6] = {
//         Ang_m2_p  Ang_m3_p     ED_p4     ED_m4        NU        ND
/* A */ {      43.4,    53.13,    10.46,     11.3,       13,       18}, // A
/* C */ {     107.3,    35.28,    5.975,    6.371,       13,       14}, // C
/* D */ {     38.43,    51.56,    9.148,    10.75,        6,       13}, // D
/* E */ {     20.73,    33.32,    12.07,    12.36,       19,       23}, // E
/* F */ {     43.56,    52.46,    11.51,    11.01,       15,       17}, // F
/* G */ {     43.19,    55.48,     10.9,    10.95,       12,       17}, // G
/* H */ {     45.73,    59.19,    10.57,    10.12,       13,       17}, // H
/* I */ {     19.82,    34.02,    12.93,    11.73,       16,       18}, // I
/* K */ {     18.88,    29.96,    13.09,    11.56,       17,       18}, // K
/* L */ {     109.8,    33.27,    6.312,    6.345,       16,       17}, // L
/* M */ {     15.49,    26.11,     13.1,    12.34,       17,       19}, // M
/* N */ {     107.5,    51.89,     6.92,    6.607,       13,       16}, // N
/* P */ {     90.53,    93.07,    9.345,    9.761,        3,       14}, // P
/* Q */ {     95.97,    62.86,     8.91,    8.455,       13,       17}, // Q
/* R */ {     88.59,    68.53,    9.457,    7.738,       13,       16}, // R
/* S */ {     110.1,    32.72,    6.305,     6.36,       15,       17}, // S
/* T */ {     29.37,    39.92,    11.98,    11.36,       15,       17}, // T
/* V */ {     109.9,    35.31,     6.39,    6.323,        4,       16}, // V
/* W */ {     25.99,    35.09,    12.79,    11.55,       16,       18}, // W
/* Y */ {     19.96,    29.55,    12.15,    11.58,       20,       20}, // Y
};

static double g_AminoSubstMx[20][20] = {
//              A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
/* A */ {  0.7892,  0.1351, -0.2798, -0.1031, -0.2414, -0.0349, -0.2369, -0.1676, -0.1370, -0.1271, -0.0086, -0.2398, -0.1653, -0.1084, -0.1764,  0.0487, -0.0833,  0.0265, -0.3102, -0.2366}, // A
/* C */ {  0.1351,  2.8679, -0.7950, -0.8209, -0.0963, -0.4257, -0.3716, -0.1128, -0.6827, -0.1382, -0.0714, -0.4635, -0.5833, -0.5846, -0.5769, -0.1758, -0.1935,  0.0764, -0.4376, -0.2110}, // C
/* D */ { -0.2798, -0.7950,  1.3116,  0.4138, -0.7560, -0.0447,  0.0030, -0.9037,  0.0569, -0.7993, -0.5213,  0.4160,  0.0158,  0.1543, -0.0422,  0.1364, -0.0389, -0.8380, -0.5691, -0.5238}, // D
/* E */ { -0.1031, -0.8209,  0.4138,  0.9592, -0.6587, -0.2708,  0.0125, -0.7331,  0.3118, -0.6282, -0.4008,  0.0835, -0.0172,  0.4309,  0.2035,  0.0633, -0.0136, -0.6137, -0.5379, -0.4107}, // E
/* F */ { -0.2414, -0.0963, -0.7560, -0.6587,  1.2745, -0.6269, -0.0786,  0.2820, -0.6191,  0.3802,  0.2613, -0.5062, -0.4567, -0.4582, -0.4711, -0.4339, -0.2906,  0.1430,  0.5755,  0.7346}, // F
/* G */ { -0.0349, -0.4257, -0.0447, -0.2708, -0.6269,  1.5076, -0.2520, -0.8630, -0.2445, -0.7884, -0.4925,  0.1138, -0.1065, -0.2652, -0.3117,  0.0689, -0.2759, -0.6958, -0.6092, -0.5309}, // G
/* H */ { -0.2369, -0.3716,  0.0030,  0.0125, -0.0786, -0.2520,  1.6489, -0.4908,  0.0295, -0.3840, -0.2077,  0.2345, -0.2399,  0.1520,  0.1141,  0.0069, -0.0576, -0.4222, -0.0164,  0.3342}, // H
/* I */ { -0.1676, -0.1128, -0.9037, -0.7331,  0.2820, -0.8630, -0.4908,  0.9921, -0.6172,  0.5415,  0.3471, -0.7603, -0.5145, -0.5376, -0.5366, -0.5888, -0.2626,  0.6891, -0.1811, -0.0746}, // I
/* K */ { -0.1370, -0.6827,  0.0569,  0.3118, -0.6191, -0.2445,  0.0295, -0.6172,  0.9358, -0.5124, -0.2459,  0.1218, -0.0239,  0.3776,  0.5009,  0.0435, -0.0147, -0.5126, -0.5373, -0.3487}, // K
/* L */ { -0.1271, -0.1382, -0.7993, -0.6282,  0.3802, -0.7884, -0.3840,  0.5415, -0.5124,  0.8825,  0.4602, -0.6394, -0.4740, -0.4272, -0.4367, -0.5368, -0.3327,  0.3400,  0.0073,  0.0090}, // L
/* M */ { -0.0086, -0.0714, -0.5213, -0.4008,  0.2613, -0.4925, -0.2077,  0.3471, -0.2459,  0.4602,  1.0332, -0.3459, -0.3377, -0.1225, -0.2351, -0.2537, -0.1153,  0.1976, -0.0226,  0.0489}, // M
/* N */ { -0.2398, -0.4635,  0.4160,  0.0835, -0.5062,  0.1138,  0.2345, -0.7603,  0.1218, -0.6394, -0.3459,  1.1595, -0.0039,  0.1817,  0.0616,  0.2404,  0.0884, -0.6361, -0.4154, -0.2360}, // N
/* P */ { -0.1653, -0.5833,  0.0158, -0.0172, -0.4567, -0.1065, -0.2399, -0.5145, -0.0239, -0.4740, -0.3377, -0.0039,  1.6216, -0.1384, -0.1318,  0.0821, -0.0559, -0.4024, -0.5042, -0.3784}, // P
/* Q */ { -0.1084, -0.5846,  0.1543,  0.4309, -0.4582, -0.2652,  0.1520, -0.5376,  0.3776, -0.4272, -0.1225,  0.1817, -0.1384,  0.8004,  0.3332,  0.1053,  0.0435, -0.4875, -0.3505, -0.2277}, // Q
/* R */ { -0.1764, -0.5769, -0.0422,  0.2035, -0.4711, -0.3117,  0.1141, -0.5366,  0.5009, -0.4367, -0.2351,  0.0616, -0.1318,  0.3332,  1.0291,  0.0233,  0.0042, -0.4699, -0.3203, -0.1899}, // R
/* S */ {  0.0487, -0.1758,  0.1364,  0.0633, -0.4339,  0.0689,  0.0069, -0.5888,  0.0435, -0.5368, -0.2537,  0.2404,  0.0821,  0.1053,  0.0233,  0.7257,  0.3676, -0.4472, -0.3899, -0.2676}, // S
/* T */ { -0.0833, -0.1935, -0.0389, -0.0136, -0.2906, -0.2759, -0.0576, -0.2626, -0.0147, -0.3327, -0.1153,  0.0884, -0.0559,  0.0435,  0.0042,  0.3676,  0.8505, -0.0762, -0.2945, -0.2150}, // T
/* V */ {  0.0265,  0.0764, -0.8380, -0.6137,  0.1430, -0.6958, -0.4222,  0.6891, -0.5126,  0.3400,  0.1976, -0.6361, -0.4024, -0.4875, -0.4699, -0.4472, -0.0762,  0.8727, -0.1976, -0.1044}, // V
/* W */ { -0.3102, -0.4376, -0.5691, -0.5379,  0.5755, -0.6092, -0.0164, -0.1811, -0.5373,  0.0073, -0.0226, -0.4154, -0.5042, -0.3505, -0.3203, -0.3899, -0.2945, -0.1976,  2.5990,  0.6562}, // W
/* Y */ { -0.2366, -0.2110, -0.5238, -0.4107,  0.7346, -0.5309,  0.3342, -0.0746, -0.3487,  0.0090,  0.0489, -0.2360, -0.3784, -0.2277, -0.1899, -0.2676, -0.2150, -0.1044,  0.6562,  1.4538}, // Y
};

vector<vector<double> > XProf::g_BinLos(XFEATS);
vector<vector<double> > XProf::g_Scores(XFEATS);

uint XProf::Get3di(const vector<double> &FeatureValues)
	{
	asserta(SIZE(FeatureValues) == XFEATS);
	vector<uint> Bins(XFEATS);
	double BestScore = -999;
	uint BestLetter = UINT_MAX;
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			double Value = FeatureValues[FeatureIndex];
			double MedValue = g_3di_Medians[Letter][FeatureIndex];
			double Diff = GetDiff(FeatureIndex, Value, MedValue);
			uint Bin = GetFeatureBin(FeatureIndex, Diff);
			Bins[FeatureIndex] = Bin;
			}
		double Score = GetScore('X', 'X', Bins);
		if (Score > BestScore)
			{
			BestScore = Score;
			BestLetter = Letter;
			}
		}
	return BestLetter;
	}

double XProf::GetDiff(uint FeatureIndex, double Value1, double Value2)
	{
// Angles Ang_m2_p2=0, Ang_m3_p3=1
	if (FeatureIndex <= 1)
		{
		double d1 = fabs(Value1 - Value2);
		double d2 = fabs(fabs(Value1 - Value2) - 180);
		double Diff = min(d1, d2);
		asserta(Diff <= 90);
		return Diff;
		}
	return fabs(Value1 - Value2);
	}

double XProf::GetScore(char Amino1, char Amino2, const vector<uint> &Bins)
	{
	asserta(SIZE(Bins) == XFEATS);
	double Score = 0;

	uint Letter1 = g_CharToLetterAmino[(byte) Amino1];
	uint Letter2 = g_CharToLetterAmino[(byte) Amino2];
	if (Letter1 < 20 && Letter2 < 20)
		Score += g_AminoSubstMx[Letter1][Letter2];

	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		{
		uint Bin = Bins[FeatureIndex];
		Score += g_Scores[FeatureIndex][Bin];
		}
	return Score;
	}

uint XProf::GetFeatureBin(uint FeatureIndex, double Diff)
	{
	asserta(FeatureIndex < XFEATS);
	for (uint Bin = 0; Bin < XBINS-1; ++Bin)
		{
		if (Diff < g_BinLos[FeatureIndex][Bin+1])
			return Bin;
		}
	return XBINS-1;
	}

void XProf::FeatureScoreBin(const string &FeatureName, uint Bin,
  double BinLo, double Score)
	{
	uint FeatureIndex = XProf::GetFeatureIndex(FeatureName);
	asserta(FeatureIndex < XFEATS);
	asserta(Bin < XBINS);
	asserta(g_BinLos[FeatureIndex][Bin] == DBL_MAX);
	asserta(g_Scores[FeatureIndex][Bin] == DBL_MAX);
	g_BinLos[FeatureIndex][Bin] = BinLo;
	g_Scores[FeatureIndex][Bin] = Score;
	}

void XProf::InitScoreTable()
	{
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		{
		g_BinLos[FeatureIndex].resize(XBINS, DBL_MAX);
		g_Scores[FeatureIndex].resize(XBINS, DBL_MAX);
		}

	FeatureScoreBin("Ang_m2_p2", 0, 0, 1.38);
	FeatureScoreBin("Ang_m2_p2", 1, 1, 1.23);
	FeatureScoreBin("Ang_m2_p2", 2, 2, 1.16);
	FeatureScoreBin("Ang_m2_p2", 3, 3.68, 1.12);
	FeatureScoreBin("Ang_m2_p2", 4, 5, 0.931);
	FeatureScoreBin("Ang_m2_p2", 5, 7.1, 0.825);
	FeatureScoreBin("Ang_m2_p2", 6, 10, 0.528);
	FeatureScoreBin("Ang_m2_p2", 7, 14.2, 0.235);
	FeatureScoreBin("Ang_m2_p2", 8, 20.5, -0.115);
	FeatureScoreBin("Ang_m2_p2", 9, 30.6, -1.17);

	FeatureScoreBin("Ang_m3_p3", 0, 0, 1.22);
	FeatureScoreBin("Ang_m3_p3", 1, 1.2, 1.01);
	FeatureScoreBin("Ang_m3_p3", 2, 2.7, 1.01);
	FeatureScoreBin("Ang_m3_p3", 3, 4.3, 0.875);
	FeatureScoreBin("Ang_m3_p3", 4, 6.2, 0.724);
	FeatureScoreBin("Ang_m3_p3", 5, 8.6, 0.523);
	FeatureScoreBin("Ang_m3_p3", 6, 11.7, 0.293);
	FeatureScoreBin("Ang_m3_p3", 7, 16.1, -0.0062);
	FeatureScoreBin("Ang_m3_p3", 8, 22.8, -0.317);
	FeatureScoreBin("Ang_m3_p3", 9, 33, -0.985);

	FeatureScoreBin("ED_p4", 0, 0, 1.37);
	FeatureScoreBin("ED_p4", 1, 0.07, 1.33);
	FeatureScoreBin("ED_p4", 2, 0.13, 1.2);
	FeatureScoreBin("ED_p4", 3, 0.2, 1.16);
	FeatureScoreBin("ED_p4", 4, 0.3, 1.08);
	FeatureScoreBin("ED_p4", 5, 0.44, 0.88);
	FeatureScoreBin("ED_p4", 6, 0.6, 0.623);
	FeatureScoreBin("ED_p4", 7, 0.9, 0.267);
	FeatureScoreBin("ED_p4", 8, 1.37, -0.143);
	FeatureScoreBin("ED_p4", 9, 2.19, -1.2);

	FeatureScoreBin("ED_m4", 0, 0, 1.37);
	FeatureScoreBin("ED_m4", 1, 0.07, 1.35);
	FeatureScoreBin("ED_m4", 2, 0.13, 1.21);
	FeatureScoreBin("ED_m4", 3, 0.2, 1.16);
	FeatureScoreBin("ED_m4", 4, 0.3, 1.09);
	FeatureScoreBin("ED_m4", 5, 0.44, 0.885);
	FeatureScoreBin("ED_m4", 6, 0.6, 0.626);
	FeatureScoreBin("ED_m4", 7, 0.9, 0.27);
	FeatureScoreBin("ED_m4", 8, 1.38, -0.153);
	FeatureScoreBin("ED_m4", 9, 2.2, -1.2);

	FeatureScoreBin("NU", 0, 0, 1.08);
	FeatureScoreBin("NU", 1, 1, 0.941);
	FeatureScoreBin("NU", 2, 2, 0.788);
	FeatureScoreBin("NU", 3, 3, 0.6);
	FeatureScoreBin("NU", 4, 4, 0.394);
	FeatureScoreBin("NU", 5, 5, 0.181);
	FeatureScoreBin("NU", 6, 6, -0.0334);
	FeatureScoreBin("NU", 7, 7, -0.242);
	FeatureScoreBin("NU", 8, 8, -0.436);
	FeatureScoreBin("NU", 9, 9, 0);

	FeatureScoreBin("ND", 0, 0, 0.776);
	FeatureScoreBin("ND", 1, 1, 0.692);
	FeatureScoreBin("ND", 2, 2, 0.538);
	FeatureScoreBin("ND", 3, 3, 0.35);
	FeatureScoreBin("ND", 4, 4, 0.119);
	FeatureScoreBin("ND", 5, 5, -0.103);
	FeatureScoreBin("ND", 6, 6, -0.33);
	FeatureScoreBin("ND", 7, 7, -0.553);
	FeatureScoreBin("ND", 8, 8, -0.767);
	FeatureScoreBin("ND", 9, 9, 0);

	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		{
		for (uint Bin = 0; Bin < XBINS; ++Bin)
			{
			asserta(g_BinLos[FeatureIndex][Bin] != DBL_MAX);
			asserta(g_Scores[FeatureIndex][Bin] != DBL_MAX);
			}
		}
	}

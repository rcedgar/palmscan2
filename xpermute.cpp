#include "myutils.h"
#include "xprof.h"
#include "xbinner.h"
#include "alpha.h"

static double g_aaMx[20][20] = {
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

static double g_XMx[20][20] = {
//              A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
/* A */ {  1.1344, -4.4433, -2.2949, -2.3964, -2.3509, -0.7956, -0.6346, -0.8931, -0.4830, -2.8896, -1.1612,  0.1859, -1.8941, -0.8284, -3.2736, -1.4376, -2.0107, -1.4270, -3.9159, -2.5608}, // A
/* C */ { -4.4433,  2.1461,  0.0494, -1.7059, -0.0803, -2.5452, -2.6322, -2.7716, -4.1066, -0.5011, -2.0272, -2.9995, -0.4640, -1.7583,  0.9157, -2.7929, -0.2626, -2.2279,  1.1122, -1.7550}, // C
/* D */ { -2.2949,  0.0494,  1.9120,  0.1387,  0.4489, -0.5293, -0.7278, -0.4099, -1.8489,  0.1487,  0.0609, -0.6819,  1.1144,  0.1313, -0.6101, -0.9407, -0.9231, -1.7436,  0.7455, -0.1176}, // D
/* E */ { -2.3964, -1.7059,  0.1387,  2.0774, -0.3370,  0.1351,  0.1370,  0.4459, -0.4653,  0.4144,  0.0436, -0.7546,  0.4964, -0.3379, -1.1935,  0.5588, -0.0732,  0.3065, -1.0129,  0.6006}, // E
/* F */ { -2.3509, -0.0803,  0.4489, -0.3370,  1.8147, -1.3457, -0.3421, -0.7367, -1.3744,  0.4322,  0.5261, -0.7432, -0.6184, -0.1423,  0.7209, -0.4603,  0.2190, -0.6620,  0.6610, -0.9585}, // F
/* G */ { -0.7956, -2.5452, -0.5293,  0.1351, -1.3457,  2.3498, -0.8686,  1.0271, -1.1409, -1.0370, -0.5839, -0.2796,  0.2347,  0.7324, -2.1704,  0.2997, -1.6562, -0.9414, -2.0041,  0.6590}, // G
/* H */ { -0.6346, -2.6322, -0.7278,  0.1370, -0.3421, -0.8686,  1.7082,  0.2585,  0.9374, -0.3916,  0.5381,  0.5650, -0.6399, -0.5289, -1.4165,  0.4990,  0.2327,  0.9114, -1.7905, -0.5376}, // H
/* I */ { -0.8931, -2.7716, -0.4099,  0.4459, -0.7367,  1.0271,  0.2585,  1.9752,  0.3799, -0.5961,  0.3954,  0.3371, -0.2316,  0.6598, -1.6915,  0.6802, -0.9519, -0.5336, -2.1324,  0.3313}, // I
/* K */ { -0.4830, -4.1066, -1.8489, -0.4653, -1.3744, -1.1409,  0.9374,  0.3799,  2.5859, -1.8304, -0.1995,  0.1436, -1.5782, -1.1148, -2.5351,  0.7528, -0.9082,  0.7440, -3.1089, -0.2956}, // K
/* L */ { -2.8896, -0.5011,  0.1487,  0.4144,  0.4322, -1.0370, -0.3916, -0.5961, -1.8304,  1.9401,  0.2770, -0.9030, -0.4912, -0.6888,  0.1154, -0.8148,  0.7315, -0.1956, -0.5161, -1.2009}, // L
/* M */ { -1.1612, -2.0272,  0.0609,  0.0436,  0.5261, -0.5839,  0.5381,  0.3954, -0.1995,  0.2770,  1.8958,  0.5848, -0.3648,  0.5996, -0.6673,  0.4342,  0.2330,  0.0873, -1.1280, -0.4286}, // M
/* N */ {  0.1859, -2.9995, -0.6819, -0.7546, -0.7432, -0.2796,  0.5650,  0.3371,  0.1436, -0.9030,  0.5848,  1.6014, -0.5690,  0.1605, -1.5773,  0.0633, -0.6175, -0.3140, -2.1178, -1.0709}, // N
/* P */ { -1.8941, -0.4640,  1.1144,  0.4964, -0.6184,  0.2347, -0.6399, -0.2316, -1.5782, -0.4912, -0.3648, -0.5690,  2.4983,  0.7221, -1.6433, -0.6694, -1.4917, -1.5373,  0.2856,  0.0542}, // P
/* Q */ { -0.8284, -1.7583,  0.1313, -0.3379, -0.1423,  0.7324, -0.5289,  0.6598, -1.1148, -0.6888,  0.5996,  0.1605,  0.7221,  2.0854, -1.0041, -0.1261, -0.6989, -0.9746, -1.4601, -0.5689}, // Q
/* R */ { -3.2736,  0.9157, -0.6101, -1.1935,  0.7209, -2.1704, -1.4165, -1.6915, -2.5351,  0.1154, -0.6673, -1.5773, -1.6433, -1.0041,  1.9117, -1.4394,  0.6053, -0.8082,  0.1444, -1.7439}, // R
/* S */ { -1.4376, -2.7929, -0.9407,  0.5588, -0.4603,  0.2997,  0.4990,  0.6802,  0.7528, -0.8148,  0.4342,  0.0633, -0.6694, -0.1261, -1.4394,  2.1669, -0.3666,  0.9136, -2.2159,  1.2628}, // S
/* T */ { -2.0107, -0.2626, -0.9231, -0.0732,  0.2190, -1.6562,  0.2327, -0.9519, -0.9082,  0.7315,  0.2330, -0.6175, -1.4917, -0.6989,  0.6053, -0.3666,  2.1098,  0.8026, -0.8720, -0.9996}, // T
/* V */ { -1.4270, -2.2279, -1.7436,  0.3065, -0.6620, -0.9414,  0.9114, -0.5336,  0.7440, -0.1956,  0.0873, -0.3140, -1.5373, -0.9746, -0.8082,  0.9136,  0.8026,  2.6166, -2.5629,  0.1339}, // V
/* W */ { -3.9159,  1.1122,  0.7455, -1.0129,  0.6610, -2.0041, -1.7905, -2.1324, -3.1089, -0.5161, -1.1280, -2.1178,  0.2856, -1.4601,  0.1444, -2.2159, -0.8720, -2.5629,  2.1357, -1.1580}, // W
/* Y */ { -2.5608, -1.7550, -0.1176,  0.6006, -0.9585,  0.6590, -0.5376,  0.3313, -0.2956, -1.2009, -0.4286, -1.0709,  0.0542, -0.5689, -1.7439,  1.2628, -0.9996,  0.1339, -1.1580,  3.1071}, // Y
};

static int GetSign(double Score)
	{
	if (Score < 0.2)
		return -1;
	else if (Score > 0.2)
		return 1;
	asserta(Score >= 0.2 && Score <= 0.2);
	return 0;
	}

static double GetSumAbsDiffFract(const vector<uint> &Order)
	{
	double SumDiffs = 0;
	double SumScores = 0;
	for (uint ia = 0; ia < 20; ++ia)
		{
		uint ix = Order[ia];
		for (uint ja = 0; ja <= ia; ++ja)
			{
			uint jx = Order[ja];
			double Scorea = g_aaMx[ia][ja];
			double Scorex = g_XMx[ix][jx];
			SumScores += fabs(Scorea + Scorex)/2;
			SumDiffs += fabs(Scorea - Scorex);
			}
		}
	return SumDiffs/SumScores;
	}

static double GetCorrectSignFract(const vector<uint> &Order)
	{
	uint n = 0;
	uint N = 0;
	for (uint ia = 0; ia < 20; ++ia)
		{
		uint ix = Order[ia];
		for (uint ja = 0; ja < ia; ++ja)
			{
			uint jx = Order[ja];
			double Scorea = g_aaMx[ia][ja];
			double Scorex = g_XMx[ix][jx];
			int Signa = GetSign(Scorea);
			int Signx = GetSign(Scorex);
			++N;
			if (Signa == Signx)
				++n;
			}
		}
	return double(n)/N;
	}

void cmd_xpermute()
	{
	XBinner::InitCentroids();

	vector<uint> Order;
	for (uint i = 0; i < 20; ++i)
		Order.push_back(i);

	vector<uint> BestOrder = Order;
	double BestCorrectSignFract = GetCorrectSignFract(Order);
	double BestSumAbsDiffFract = GetSumAbsDiffFract(Order);
	ProgressLog("Init signs %.3g, diffs %.3g\n",
	  BestCorrectSignFract, BestSumAbsDiffFract);

	uint Iter = 0;
	uint LastImprovedIter = 0;
	uint Improvements = 0;
	const uint MAXITERS = 1024*1024;
	for (uint Iter = 0; Iter < MAXITERS; ++Iter)
		{
		bool DoShuffle = (Iter < 1024);
		ProgressStep(Iter, MAXITERS, "Last %u, signs %.3g, diff %.3g, improves %u shuffle %c",
		  LastImprovedIter, BestCorrectSignFract, BestSumAbsDiffFract, Improvements, yon(DoShuffle));
		if (DoShuffle)
			Shuffle(Order);
		else
			{
			uint i1 = randu32()%20;
			uint i2 = randu32()%19;
			if (i2 == i1)
				i2 = (i2+1)%20;
			swap(Order[i1], Order[i2]);
			}
		double CorrectSignFract = GetCorrectSignFract(Order);
		double SumAbsDiffFract = GetSumAbsDiffFract(Order);
		if (CorrectSignFract > BestCorrectSignFract ||
		  (CorrectSignFract == BestCorrectSignFract && SumAbsDiffFract > BestSumAbsDiffFract))
			{
			++Improvements;
			BestOrder = Order;
			BestCorrectSignFract = CorrectSignFract;
			BestSumAbsDiffFract = SumAbsDiffFract;
			LastImprovedIter = Iter;
			}
		}

	Log("// CorrectSignFract %.3g, SumAbsDiffFract %.3g\n",
	  BestCorrectSignFract, BestSumAbsDiffFract);
	Log("static double g_PermutedMx[20][20] = {\n");
	Log("//      ");
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		Log("        %c", c);
		}
	Log("\n");
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		Log("/* %c */ {", c);
		for (uint j = 0; j < 20; ++j)
			{
			double Score = g_XMx[Order[i]][Order[j]];
			Log(" %7.4f", Score);
			if (j != 19)
				Log(",");
			}
		Log("}, // %c\n", c);
		}
	Log("};\n");

	const vector<vector<double> > &ValuesVec = XBinner::m_ValuesVec;
	for (uint i = 0; i < 20; ++i)
		{
		const vector<double> &v = ValuesVec[Order[i]];
		Log("DefineCentroid(%2u", i);
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			Log(", %8.2f", v[FeatureIndex]);
		Log("); // xpermute[%u]\n", Order[i]);
		}
	}

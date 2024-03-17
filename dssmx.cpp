#include "myutils.h"
#include "dss.h"

// NUDX -- Sphere density
static double Mx0[8][8] = {
//                A        C        D        E        F        G        H        I
/*   A */ {  1.7305,  0.7576, -0.3659, -1.3799, -2.2241, -2.9314, -3.6245, -1.4450}, // A
/*   C */ {  0.7576,  1.1985,  0.6384, -0.2670, -1.1905, -1.9621, -2.6912, -2.6354}, // C
/*   D */ { -0.3659,  0.6384,  1.0001,  0.5656, -0.2358, -1.0797, -1.8985, -2.6200}, // D
/*   E */ { -1.3799, -0.2670,  0.5656,  0.8852,  0.5430, -0.1948, -1.0835, -2.1497}, // E
/*   F */ { -2.2241, -1.1905, -0.2358,  0.5430,  0.8602,  0.5192, -0.2804, -1.4212}, // F
/*   G */ { -2.9314, -1.9621, -1.0797, -0.1948,  0.5192,  0.8960,  0.5143, -0.6096}, // G
/*   H */ { -3.6245, -2.6912, -1.8985, -1.0835, -0.2804,  0.5143,  1.0356,  0.4442}, // H
/*   I */ { -1.4450, -2.6354, -2.6200, -2.1497, -1.4212, -0.6096,  0.4442,  1.4337}, // I
};
static double Freqs0[8] = {
	0.1023, // A
	0.1177, // C
	0.1181, // D
	0.1199, // E
	0.1248, // F
	0.1306, // G
	0.1362, // H
	0.1505, // I
};
static uint N0 = 8;

// SSX -- secondary structure at Pos
static double Mx1[4][4] = {
//              A        C        D        E
/* A */ {  0.9236, -3.0191, -0.6890, -1.1623}, // A
/* C */ { -3.0191,  1.1476, -1.6956, -0.5019}, // C
/* D */ { -0.6890, -1.6956,  1.6978, -0.0853}, // D
/* E */ { -1.1623, -0.5019, -0.0853,  0.6855}, // E
};
static double Freqs1[4] = {
    0.3331, // A
    0.2424, // C
    0.08769, // D
    0.3369, // E
};
static uint N1 = 4;

// SSX2 -- secondary structure at closest non-local (Pos2)
static double Mx2[4][4] = {
//              A        C        D        E
/* A */ {  1.0094, -1.4946, -0.1778, -0.6197}, // A
/* C */ { -1.4946,  0.7252, -0.9605, -0.3759}, // C
/* D */ { -0.1778, -0.9605,  1.2658, -0.0244}, // D
/* E */ { -0.6197, -0.3759, -0.0244,  0.5106}, // E
};
static double Freqs2[4] = {
    0.243, // A
    0.3306, // C
    0.09882, // D
    0.3275, // E
};
static uint N2 = 4;

// SSD2 -- distance to Pos2
static double Mx3[4][4] = {
//              A        C        D        E
/* A */ {  1.0366, -0.2202, -1.5137, -2.1158}, // A
/* C */ { -0.2202,  0.5117, -0.1705, -1.0896}, // C
/* D */ { -1.5137, -0.1705,  0.6648,  0.0602}, // D
/* E */ { -2.1158, -1.0896,  0.0602,  1.1004}, // E
};
static double Freqs3[4] = {
	0.2346, // A
	0.3269, // C
	0.2347, // D
	0.2038, // E
};
static uint N3 = 4;

// CMAA -- compressed alphabet AST,C,DN,EHKQR,FWY,G,ILMV,P
static double Mx4[8][8] = {
//                A        C        D        E        F        G        H        I
/*   A */ {  0.3931, -0.0364, -0.0344, -0.0424, -0.2849, -0.0617, -0.2274, -0.0553}, // A
/*   C */ { -0.0364,  2.8679, -0.6430, -0.6420, -0.1889, -0.4257, -0.0576, -0.5833}, // C
/*   D */ { -0.0344, -0.6430,  0.9369,  0.1440, -0.5174,  0.0246, -0.7373,  0.0076}, // D
/*   E */ { -0.0424, -0.6420,  0.1440,  0.4986, -0.3658, -0.2703, -0.5085, -0.0827}, // E
/*   F */ { -0.2849, -0.1889, -0.5174, -0.3658,  1.0811, -0.5858,  0.1095, -0.4328}, // F
/*   G */ { -0.0617, -0.4257,  0.0246, -0.2703, -0.5858,  1.5076, -0.7524, -0.1065}, // G
/*   H */ { -0.2274, -0.0576, -0.7373, -0.5085,  0.1095, -0.7524,  0.6280, -0.4509}, // H
/*   I */ { -0.0553, -0.5833,  0.0076, -0.0827, -0.4328, -0.1065, -0.4509,  1.6216}, // I
};
static double Freqs4[8] = {
    0.1969, // A
    0.01418, // C
    0.09424, // D
    0.2336, // E
    0.09094, // F
    0.07123, // G
    0.2577, // H
    0.04122, // I
};
static uint N4 = 8;

uint DSS::GetIdx_WithCMAA(uint i0, uint i1, uint i2, uint i3, uint i4)
	{
	return i0 + i1*N0 + i2*N0*N1 + i3*N0*N1*N2 + i4*N0*N1*N2*N3;
	}

uint DSS::GetIdx_NoCMAA(uint i0, uint i1, uint i2, uint i3)
	{
	return i0 + i1*N0 + i2*N0*N1 + i3*N0*N1*N2;
	}

double DSS::GetScore_WithCMAA(
  uint i0, uint i1, uint i2, uint i3, uint i4,
  uint j0, uint j1, uint j2, uint j3, uint j4)
	{
	if (i0 >= N0 || j0 > N0) return 0;
	if (i1 >= N1 || j1 > N1) return 0;
	if (i2 >= N2 || j2 > N2) return 0;
	if (i3 >= N3 || j3 > N3) return 0;
	if (i4 >= N4 || j4 > N4) return 0;

	double Score0 = Mx0[i0][j0];
	double Score1 = Mx1[i1][j1];
	double Score2 = Mx2[i2][j2];
	double Score3 = Mx3[i3][j3];
	double Score4 = Mx4[i4][j4];
	return Score0 + Score1 + Score2 + Score3 + Score4;
	}

double DSS::GetScore_NoCMAA(
  uint i0, uint i1, uint i2, uint i3,
  uint j0, uint j1, uint j2, uint j3)
	{
	if (i0 >= N0 || j0 > N0) return 0;
	if (i1 >= N1 || j1 > N1) return 0;
	if (i2 >= N2 || j2 > N2) return 0;
	if (i3 >= N3 || j3 > N3) return 0;

	double Score0 = Mx0[i0][j0];
	double Score1 = Mx1[i1][j1];
	double Score2 = Mx2[i2][j2];
	double Score3 = Mx3[i3][j3];
	return Score0 + Score1 + Score2 + Score3;
	}

#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"

double get_norm(const coords_t a)
	{
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
	}

coords_t normalize(const coords_t &a)
	{
	coords_t c;
	double norm = get_norm(a);
	assert(norm > 1e-3);
	c.x = a.x / norm;
	c.y = a.y / norm;
	c.z = a.z / norm;
	assert(feq(get_norm(c), 1));
	return c;
	}

coords_t subtract(const coords_t &a, const coords_t &b)
	{
	coords_t c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return c;
	}

static coords_t add(const coords_t a, const coords_t b)
	{
	coords_t c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return c;
	}

static double dot_product(const coords_t a, const coords_t b)
	{
	return a.x*b.x + a.y*b.y + a.z*b.z;
	}

coords_t cross_product(const coords_t &a, const coords_t &b)
	{
	coords_t c;
	c.x = a.y*b.z - a.z*b.y;
	c.y = a.z*b.x - a.x*b.z;
	c.z = a.x*b.y - a.y*b.x;
	return c;
	}

double get_angle(const coords_t &a, const coords_t &b)
	{
	double moda = get_norm(a);
	double modb = get_norm(b);
	double dotab = dot_product(a, b);
	assert(moda > 0);
	assert(modb > 0);
	double theta = acos(dotab/(moda*modb));
	return theta;
	}

//static void GetAngles(coords_t a, coords_t b, coords_t c,
//					  double &theta_bc, double &theta_vc)
//	{
//	theta_bc = get_angle(b, c);
//	coords_t v = cross_product(a, b);
//	theta_vc = get_angle(v, c);
//	}

// https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
// axis defines the axis for rotation
coords_t rotate_around_vector(const coords_t &v,
							  const coords_t &axis, double theta_rad)
	{
	coords_t k = normalize(axis);

	assert(feq(get_norm(k), 1));

	coords_t k_cross_v = cross_product(v, k);

	double sin_theta = sin(theta_rad);
	double cos_theta = cos(theta_rad);

	coords_t v_cos_theta;
	v_cos_theta.x = v.x*cos_theta;
	v_cos_theta.y = v.y*cos_theta;
	v_cos_theta.z = v.z*cos_theta;

	coords_t k_cross_v_sin_theta;
	k_cross_v_sin_theta.x = k_cross_v.x*sin_theta;
	k_cross_v_sin_theta.y = k_cross_v.y*sin_theta;
	k_cross_v_sin_theta.z = k_cross_v.z*sin_theta;

	double k_dot_v = dot_product(k, v);
	coords_t kkvc;
	kkvc.x = k.x*k_dot_v*(1 - cos_theta);
	kkvc.y = k.y*k_dot_v*(1 - cos_theta);
	kkvc.z = k.z*k_dot_v*(1 - cos_theta);

	coords_t v_rot;
	v_rot.x = v_cos_theta.x + k_cross_v_sin_theta.x + kkvc.x;
	v_rot.y = v_cos_theta.y + k_cross_v_sin_theta.y + kkvc.y;
	v_rot.z = v_cos_theta.z + k_cross_v_sin_theta.z + kkvc.z;

	return v_rot;
	}

//static coords_t GetVec(const PDBChain &Chain, uint Pos)
//	{
//	asserta(Pos > 0);
//	coords_t c0 = Chain.GetCoords(Pos-1);
//	coords_t c = Chain.GetCoords(Pos);
//	coords_t v;
//	v.x = c.x - c0.x;
//	v.y = c.y - c0.y;
//	v.z = c.z - c0.z;
//	return v;
//	}

static void LogCoord(double x)
	{
	if (x > -1e-6 && x < 1e-6)
		x = 0;
	double fint;
	double ffract = modf(x, &fint);
	double f = fabs(ffract);
	if (f < 1e-4 || f > 0.9999)
		Log("%.0f", x);
	else
		Log("%.2f", x);
	}

void LogCoords(const char *Name, coords_t c)
	{
	Log("%s=(", Name);
	LogCoord(c.x);
	Log(", ");
	LogCoord(c.y);
	Log(", ");
	LogCoord(c.z);
	Log(")");
	Log("\n");
	}

static void calculate_theta_phi(coords_t A, coords_t B, coords_t C, coords_t D,
								  double &theta_rad, double &phi_rad)
	{
	coords_t ab = subtract(B, A);
	coords_t bc = subtract(C, B);
	coords_t cd2 = subtract(D, C);
	coords_t vert_abc = cross_product(ab, bc);
	theta_rad = get_angle(bc, cd2);

	coords_t vert_bcd2 = cross_product(bc, cd2);

	phi_rad = PI - get_angle(vert_abc, vert_bcd2);
	}

static double randangle()
	{
	double x = 0.1*PI*(randu32()%7829)/7829.0;
	return x;
	}

static double randcoord()
	{
	double x = 10.0*(randu32()%7829)/7829.0;
	return x;
	}

coords_t get_unit_cd(coords_t A, coords_t B, coords_t C,
			   double theta_rad, double phi_rad)
	{
	coords_t ab = subtract(B, A);
	coords_t bc = subtract(C, B);
	coords_t vert_abc = cross_product(ab, bc);

	coords_t cd1 = rotate_around_vector(bc, vert_abc, theta_rad);
	coords_t cd2 = rotate_around_vector(cd1, bc, phi_rad);
	coords_t unit_cd = normalize(cd2);
	return unit_cd;
	}

static bool TestAnglesRand1()
	{
	coords_t A(randcoord(), randcoord(), randcoord());
	coords_t B(randcoord(), randcoord(), randcoord());
	coords_t C(randcoord(), randcoord(), randcoord());

	coords_t ab = subtract(B, A);
	coords_t bc = subtract(C, B);
	coords_t vert_abc = cross_product(ab, bc);

	double theta_rad = randangle();
	double phi_rad = randangle();
	coords_t cd1 = rotate_around_vector(bc, vert_abc, theta_rad);
	coords_t cd2 = rotate_around_vector(cd1, bc, phi_rad);
	coords_t D = add(C, cd2);

	double theta2_rad, phi2_rad;
	calculate_theta_phi(A, B, C, D, theta2_rad, phi2_rad);

	double theta_deg = degrees(theta_rad);
	double phi_deg = degrees(phi_rad);

	double theta2_deg = degrees(theta2_rad);
	double phi2_deg = degrees(phi2_rad);

	ProgressLog("theta %.1f, %.1f   phi %.1f, %.1f, %.1f\n",
		theta_deg, theta2_deg, phi_deg, phi2_deg, 180 - phi2_deg);

	bool Ok = feq(theta_rad, theta2_rad) && feq(phi_rad, phi2_rad);
	return Ok;
	}

static void TestAngles()
	{
	coords_t A(0, 0, 0);
	coords_t B(1, 0, 1);
	coords_t C(2, 0, 0);
	LogCoords("A", A);
	LogCoords("B", B);
	LogCoords("C", C);

	coords_t ab = subtract(B, A);
	coords_t bc = subtract(C, B);
	LogCoords("ab", ab);
	LogCoords("bc", bc);

	coords_t unit_ab = normalize(ab);
	coords_t unit_bc = normalize(bc);
	LogCoords("unit_ab", unit_ab);
	LogCoords("unit_bc", unit_bc);

	coords_t vert_abc = cross_product(ab, bc);
	LogCoords("vert_abc", vert_abc);

	coords_t unit_vert_abc = normalize(vert_abc);
	LogCoords("unit_vert_abc", unit_vert_abc);

	double theta_deg = 32;
	double phi_deg = 10;
	double theta_rad = radians(theta_deg);
	double phi_rad = radians(phi_deg);
	coords_t cd1 = rotate_around_vector(unit_bc, unit_vert_abc, theta_rad);
	Log("theta_deg = %.1f\n", theta_deg);
	LogCoords("cd1", cd1);

	coords_t unit_cd1 = normalize(cd1);
	LogCoords("unit_cd1", unit_cd1);

	coords_t cd2 = rotate_around_vector(unit_cd1, unit_bc, phi_rad);
	Log("phi_deg = %.1f\n", phi_deg);
	coords_t unit_cd2 = normalize(cd2);
	LogCoords("unit_cd2", unit_cd2);

	coords_t D = add(C, unit_cd2);
	LogCoords("D", D);

	double theta2_rad, phi2_rad;
	calculate_theta_phi(A, B, C, D, theta2_rad, phi2_rad);
	Log("theta2_deg = %.1f\n", degrees(theta_rad));
	Log("phi2_deg = %.1f\n", degrees(phi_rad));
	}

static void TestAnglesRand()
	{
	uint N = 1000;
	uint n = 0;
	for (uint i = 0; i < N; ++i)
		{
		bool Ok = TestAnglesRand1();
		if (Ok) ++n;
		}
	ProgressLog("%u / %u ok\n", n, N);
	return;
	}


static void Angles(FILE *f, const PDBChain &Chain)
	{
	static uint Counter;
	const uint L = Chain.GetSeqLength();
	for (uint Pos = 4; Pos < L; ++Pos)
		{
		coords_t A = Chain.GetCoords(Pos-3);
		coords_t B = Chain.GetCoords(Pos-2);
		coords_t C = Chain.GetCoords(Pos-1);
		coords_t D = Chain.GetCoords(Pos);
		double theta_rad, phi_rad;
		calculate_theta_phi(A, B, C, D, theta_rad, phi_rad);
		double theta_deg = degrees(theta_rad);
		double phi_deg = degrees(phi_rad);
		if (Counter++%5 == 0)
			Log("\n");
		Log(" {%5.1f, %5.1f},", theta_deg, phi_deg);
		if (f != 0)
			fprintf(f, "%.1f\t%.1f\n", theta_deg, phi_deg);
		}
	}

void cmd_angles()
	{
	const string &InputFN = g_Arg1;
	FILE *f = CreateStdioFile(opt_output);
	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);
	const uint N = SIZE(Chains);
	for (uint i = 0; i < N; ++i)
		Angles(f, *Chains[i]);
	CloseStdioFile(f);
	}

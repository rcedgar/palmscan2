#include "myutils.h"
#include "abcxyz.h"
#include "pdbchain.h"

static void mul_mx_v(double mx[3][3], double v[3], double result[3])
	{
	for (int i = 0; i < 3; ++i)
		{
		double sum = 0;
		for (int j = 0; j < 3; ++j)
			sum += mx[i][j]*v[i];
		result[i] = sum;
		}
	}

static void get_mx_rotx(double theta, double mx[3][3])
	{
	double s = sin(theta);
	double c = cos(theta);

	mx[0][0] = 1;
	mx[0][1] = 0;
	mx[0][2] = 0;

	mx[1][0] = 0;
	mx[1][1] = c;
	mx[1][2] = s;
	
	mx[2][0] = 0;
	mx[2][1] = -s;
	mx[2][2] = c;
	}

static void get_mx_roty(double theta, double mx[3][3])
	{
	double s = sin(theta);
	double c = cos(theta);

	mx[0][0] = c;
	mx[0][1] = 0;
	mx[0][2] = s;

	mx[1][0] = 0;
	mx[1][1] = 1;
	mx[1][2] = 0;
	
	mx[2][0] = -s;
	mx[2][1] = 0;
	mx[2][2] = c;
	}

static void get_mx_rotz(double theta, double mx[3][3])
	{
	double s = sin(theta);
	double c = cos(theta);

	mx[0][0] = c;
	mx[0][1] = -s;
	mx[0][2] = 0;

	mx[1][0] = s;
	mx[1][1] = c;
	mx[1][2] = 0;
	
	mx[2][0] = 0;
	mx[2][1] = 0;
	mx[2][2] = 1;
	}

void cmd_rotate3d()
	{
	Die("TODO");
	}

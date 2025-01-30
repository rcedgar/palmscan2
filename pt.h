#pragma once

template <typename T> class tpl_pt
	{
public:
	T x;
	T y;
	T z;

public:
	tpl_pt()
		{
		x = 0;
		y = 0;
		z = 0;
		}

	tpl_pt(int ax, int ay, int az)
		{
		x = ax;
		y = ay;
		z = az;
		}

	tpl_pt(const tpl_pt &rhs)
		{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		}

	bool operator<(const tpl_pt &rhs) const
		{
		if (x < rhs.x)
			return true;
		if (x > rhs.x)
			return false;
		if (y < rhs.y)
			return true;
		if (y > rhs.y)
			return false;
		return z < rhs.z;
		}
	};

typedef tpl_pt<int> intpt_t;
typedef tpl_pt<double> coords_t;

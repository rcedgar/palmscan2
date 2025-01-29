#pragma once

class PDBChain;
#include <set>

class pt_t
	{
public:
	int x;
	int y;
	int z;

public:
	pt_t()
		{
		x = 0;
		y = 0;
		z = 0;
		}

	pt_t(int ax, int ay, int az)
		{
		x = ax;
		y = ay;
		z = az;
		}

	pt_t(const pt_t &rhs)
		{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		}

	bool operator<(const pt_t &rhs) const
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

class FakeChain
	{
public:
	PDBChain m_Chain;
	double m_Size = 5;
	vector<const PDBChain *> m_Frags;
	vector<pt_t> m_AppendPts;

	vector<pt_t> m_CubePts;
	set<int64> m_HashSet;

	void GetCube(double x0, double y0, double z0,
				 double x1, double y1, double z1,
				 pt_t &p) const;
	void GetCubeSet(const PDBChain &Chain,
					set<pt_t> &CubeSet) const;
	bool HasOverlap(pt_t p0, const set<pt_t> &Cubes) const;
	bool TryAppendFrag(const PDBChain &Frag);
	uint64 GetHash(pt_t p) const;
	uint64 GetHash(int x, int y, int z) const;
	void GetAppendPt(pt_t &p) const;

private:
	void AppendFrag(const PDBChain &Frag, pt_t p,
					const set<pt_t> &Cubes);
	};

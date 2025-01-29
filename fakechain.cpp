#include "myutils.h"
#include "pdbchain.h"
#include "fakechain.h"

void FakeChain::GetCube(double x0, double y0, double z0,
				 double x1, double y1, double z1,
				 pt_t &p) const
	{
	p.x = int(round((x1 - x0)/m_Size));
	p.y = int(round((y1 - y0)/m_Size));
	p.z = int(round((y1 - y0)/m_Size));
	}

void FakeChain::GetCubeSet(const PDBChain &Chain,
						   set<pt_t> &CubeSet) const
	{
	CubeSet.clear();
	const uint L = Chain.GetSeqLength();
	const double x0 = Chain.m_Xs[0];
	const double y0 = Chain.m_Ys[0];
	const double z0 = Chain.m_Zs[0];
	for (uint i = 0; i < L; ++i)
		{
		double x = Chain.m_Xs[i];
		double y = Chain.m_Ys[i];
		double z = Chain.m_Zs[i];
		pt_t Cube;
		GetCube(x0, y0, y0,
				x, y, x, Cube);

		CubeSet.insert(Cube);
		}
	}

void FakeChain::GetAppendPt(pt_t &p) const
	{
	uint n = SIZE(m_CubePts);
	if (n == 0)
		{
		p.x = 0;
		p.y = 0;
		p.z = 0;
		return;
		}
	asserta(n >= 2);
	pt_t pt1 = m_CubePts[n-2];
	pt_t pt2 = m_CubePts[n-1];

	p.x = pt2.x + (pt2.x - pt1.x);
	p.y = pt2.y + (pt2.y - pt1.y);
	p.z = pt2.z + (pt2.z - pt1.z);
	}

uint64 FakeChain::GetHash(int x, int y, int z) const
	{
	pt_t p(x, y, z);
	return GetHash(p);
	}

uint64 FakeChain::GetHash(pt_t p) const
	{
	static int64 a = 6364136223846793005;
	static int64 c = 1442695040888963407;
	int64 h = int64(p.x + 100*(p.y + 100) + 10000*(p.z + 100));
	h = h*a + c;
	return h;
	}

bool FakeChain::HasOverlap(pt_t p0, const set<pt_t> &Cubes) const
	{
	pt_t p = p0;
	for (set<pt_t>::const_iterator iter = Cubes.begin();
		 iter != Cubes.end(); ++iter)
		{
		pt_t Cube = *iter;
		Cube.x += p.x;
		Cube.y += p.y;
		Cube.z += p.z;
		int64 h = GetHash(Cube);
		if (m_HashSet.find(h) != m_HashSet.end())
			return true;
		}
	return false;
	}

bool FakeChain::TryAppendFrag(const PDBChain &Frag)
	{
	set<pt_t> CubeSet;
	GetCubeSet(Frag, CubeSet);
	pt_t appendpt;
	GetAppendPt(appendpt);
	if (HasOverlap(appendpt, CubeSet))
		return false;
	AppendFrag(Frag, appendpt, CubeSet);
	Log("Append ok %s\n", Frag.m_Seq.c_str());
	return true;
	}

void FakeChain::AppendFrag(const PDBChain &Frag, pt_t p,
					const set<pt_t> &CubeSet)
	{
	m_Frags.push_back(&Frag);
	m_AppendPts.push_back(p);
	for (set<pt_t>::const_iterator iter = CubeSet.begin();
		 iter != CubeSet.end(); ++iter)
		{
		pt_t Cube = *iter;
		Cube.x += p.x;
		Cube.y += p.y;
		Cube.z += p.z;
		m_CubePts.push_back(Cube);
		int64 h = GetHash(Cube);
		if (m_HashSet.find(h) != m_HashSet.end())
			Die("Overlap");
		m_HashSet.insert(h);
		}

	m_Chain.m_Seq += Frag.m_Seq;
	m_Chain.m_Xs.insert(m_Chain.m_Xs.end(), Frag.m_Xs.begin(), Frag.m_Xs.end());
	m_Chain.m_Ys.insert(m_Chain.m_Ys.end(), Frag.m_Ys.begin(), Frag.m_Ys.end());
	m_Chain.m_Zs.insert(m_Chain.m_Zs.end(), Frag.m_Zs.begin(), Frag.m_Zs.end());
	}

#include "myutils.h"
#include "pdbchain.h"
#include "fakechain.h"

void FakeChain::LogMe() const
	{
	const uint L = m_Chain.GetSeqLength();
	Log("FakeChain::LogMe() >%s(%u)\n",
		m_Chain.m_Label.c_str(), L);
	const string &Seq = m_Chain.m_Seq;
	int64 lasth = 0;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		coords_t c;
		m_Chain.GetCoords(Pos, c);

		Log("%5u |%c|  X=%5.1f  Y=%5.1f  Z=%5.1f",
			Pos, Seq[Pos], c.x, c.y, c.z);

		intpt_t Cube;
		GetCube(c, Cube);
		Log("  Cube=(%4d, %4d, %4d)",
			Cube.x, Cube.y, Cube.z);

		int64 h = GetHash(Cube);
		if (h != lasth)
			{
			Log("  %16llx", h);
		
			map<int64, vector<uint> >::const_iterator iter = m_HashToPosVec.find(h);
			asserta(iter != m_HashToPosVec.end());

			const vector<uint> &PosVec = iter->second;
			const uint n = SIZE(PosVec);
			Log(" posvec(%u)", n);
			for (uint i = 0; i < n; ++i)
				Log(" %u", PosVec[i]);
			}
		lasth = h;
		Log("\n");
		}
	}

bool FakeChain::IsOccupied(coords_t c) const
	{
	intpt_t Cube;
	GetCube(c, Cube);
	int64 h = GetHash(Cube);
	return m_HashSet.find(h) != m_HashSet.end();
	}

void FakeChain::GetCube(coords_t c, intpt_t &Cube) const
	{
	Cube.x = int(round(c.x)/m_Size);
	Cube.y = int(round(c.y)/m_Size);
	Cube.z = int(round(c.z)/m_Size);
	}

bool FakeChain::GetAppendCoords(coords_t &Coords) const
	{
	const uint L = m_Chain.GetSeqLength();
	if (L == 0)
		{
		Coords.x = 0;
		Coords.y = 0;
		Coords.z = 0;
		return true;
		}
	asserta(L >= 2);
	coords_t c1;
	coords_t c2;
	m_Chain.GetCoords(L-2, c1);
	m_Chain.GetCoords(L-1, c2);

	double dx = c2.x - c1.x;
	double dy = c2.y - c1.y;
	double dz = c2.z - c1.z;

	for (int Try = 0; Try < 100; ++Try)
		{
		double fx = randu32()%200/100.0;
		double fy = randu32()%200/100.0;
		double fz = randu32()%200/100.0;
		Coords.x = c2.x + dx*fx;
		Coords.y = c2.y + dy*fx;
		Coords.z = c2.z + dz*fx;
		if (!IsOccupied(Coords))
			return true;
		}
	return false;
	}

uint64 FakeChain::GetHash(int x, int y, int z) const
	{
	intpt_t p(x, y, z);
	return GetHash(p);
	}

uint64 FakeChain::GetHash(intpt_t p) const
	{
	static int64 a = 6364136223846793005;
	static int64 c = 1442695040888963407;
	int64 h = int64(p.x + 100*(p.y + 100) + 10000*(p.z + 100));
	h = h*a + c;
	return h;
	}

uint FakeChain::FindOverlap(const PDBChain &Chain) const
	{
	const uint L = Chain.GetSeqLength();
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		coords_t Coords;
		Chain.GetCoords(Pos, Coords);

		intpt_t Cube;
		GetCube(Coords, Cube);
		int64 h = GetHash(Cube);
		if (m_HashSet.find(h) != m_HashSet.end())
			return Pos;
		}
	return UINT_MAX;
	}

bool FakeChain::TryAppendFrag(const PDBChain &Frag)
	{
	const uint L = Frag.GetSeqLength();
	asserta(L >= 4);
	asserta(Frag.m_Xs[0] == 0);
	asserta(Frag.m_Ys[0] == 0);
	asserta(Frag.m_Zs[0] == 0);

	coords_t a;
	bool Ok = GetAppendCoords(a);
	if (!Ok)
		return false;

	PDBChain *NewFrag = new PDBChain;
	*NewFrag = Frag;
	NewFrag->SetOrigin(a.x, a.y, a.z);

	uint OvPos = FindOverlap(*NewFrag);
	if (OvPos != UINT_MAX)
		{
		coords_t c;
		NewFrag->GetCoords(OvPos, c);

		intpt_t Cube;
		GetCube(c, Cube);

		Log("\n");
		Log("AppendCoords %.1f, %.1f, %.1f\n", a.x, a.y, a.z);
		Log("OvPos = %u", OvPos);
		Log(" x=%.1f, y=%.1f, z=%.1f", c.x, c.y, c.z);
		Log(" cube=%d,%d,%d", Cube.x, Cube.y, Cube.z);
		Log("\n");
		delete NewFrag;
		return false;
		}
	AppendFrag(*NewFrag, a);
	Log("Append ok %s\n", Frag.m_Seq.c_str());
	return true;
	}

void FakeChain::AppendFrag(const PDBChain &Frag, coords_t a)
	{
	const uint L = Frag.GetSeqLength();
	uint FragIdx = SIZE(m_Frags);
	asserta(SIZE(m_AppendCoordsVec) == FragIdx);

	m_Frags.push_back(&Frag);
	m_AppendCoordsVec.push_back(a);
	set<int64> NewSet;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		coords_t Coords;
		Frag.GetCoords(Pos, Coords);

		intpt_t Cube;
		GetCube(Coords, Cube);
		int64 h = GetHash(Cube);
		NewSet.insert(h);
		}
	for (set<int64>::const_iterator iter = NewSet.begin();
		 iter != NewSet.end(); ++iter)
		{
		int64 h = *iter;
		if (m_HashSet.find(h) != m_HashSet.end())
			Die("Overlap");
		m_HashSet.insert(h);
		}

	for (uint FragPos = 0; FragPos < L; ++FragPos)
		{
		coords_t Coords;
		Frag.GetCoords(FragPos, Coords);

		intpt_t Cube;
		GetCube(Coords, Cube);
		int64 h = GetHash(Cube);

		map<int64, vector<uint> >::iterator iter = m_HashToPosVec.find(h);
		if (iter == m_HashToPosVec.end())
			{
			vector<uint> EmptyPosVec;
			m_HashToPosVec[h] = EmptyPosVec;
			}
		uint Pos = SIZE(m_Chain.m_Xs);
		m_HashToPosVec[h].push_back(Pos);

		m_Chain.m_Xs.push_back(Coords.x);
		m_Chain.m_Ys.push_back(Coords.y);
		m_Chain.m_Zs.push_back(Coords.z);

		m_PosToFragIdx.push_back(FragIdx);
		m_PosToFragPos.push_back(FragPos);
		}
	m_Chain.m_Seq += Frag.m_Seq;
	}

void FakeChain::Validate() const
	{
	const uint FragCount = SIZE(m_Frags);
	asserta(SIZE(m_AppendCoordsVec) == FragCount);

	const uint L = m_Chain.GetSeqLength();
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		coords_t c;
		m_Chain.GetCoords(Pos, c);

		uint FragIdx = m_PosToFragIdx[Pos];
		uint FragPos = m_PosToFragPos[Pos];
		asserta(FragIdx < FragCount);
		const PDBChain &Frag = *m_Frags[FragIdx];
		const uint FragL = Frag.GetSeqLength();
		asserta(FragPos < FragL);

		intpt_t Cube;
		GetCube(c, Cube);
		int64 h = GetHash(Cube);
		asserta(m_HashSet.find(h) != m_HashSet.end());
		map<int64, vector<uint> >::const_iterator iter = m_HashToPosVec.find(h);
		asserta(iter != m_HashToPosVec.end());
		const vector<uint> &PosVec = iter->second;
		const uint n = SIZE(PosVec);
		asserta(n > 0);
		bool Found = 0;
		for (uint i = 0; i < n; ++i)
			{
			if (PosVec[i] == Pos)
				{
				Found = true;
				break;
				}
			}
		asserta(Found);
		}
	}

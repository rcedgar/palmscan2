#include "myutils.h"
#include "pdbchain.h"
#include "fakechain.h"

static const double MEDIAN_CALPHA_DIST = 3.81;

void RotateChain(PDBChain &Chain, double alpha, double beta, double gamma);

static double get_rand_radians()
	{
	const double TWOPI = 2*3.1415926535;
	double zero_to_one = (randu32()%1500450271)/1500450271.0;
	assert(zero_to_one >= 0 && zero_to_one < 1);
	double radians = zero_to_one*TWOPI;
	return radians;
	}

static double get_dist(const coords_t &c1, const coords_t &c2)
	{
	double dx = c1.x - c2.x;
	double dy = c1.y - c2.y;
	double dz = c1.z - c2.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
	}

// x = r * sin(theta) * cos(phi)
// y = r * sin(theta) * sin(phi)
// z = r * cos(theta)
static void GetRandomCoordsByDistance(const coords_t &c, double d, coords_t &rc)
	{
	double theta = get_rand_radians();
	double phi = get_rand_radians();
	double sin_theta = sin(theta);
	double cos_theta = cos(theta);
	double sin_phi = sin(phi);
	double cos_phi = cos(phi);

	rc.x = c.x + d*sin_theta*cos_phi;
	rc.y = c.y + d*sin_theta*sin_phi;
	rc.z = c.z + d*cos_theta;

#if DEBUG
	{
	double dist = get_dist(c, rc);
	double diff = fabs(dist - d);
	assert(diff < 0.1);
	}
#endif
	}

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
	coords_t cterm;
	m_Chain.GetCoords(L-1, cterm);
	for (int Try = 0; Try < 100; ++Try)
		{
		GetRandomCoordsByDistance(cterm, MEDIAN_CALPHA_DIST, Coords);
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
		{
		Log("AppendCoords failed\n");
		return false;
		}

	PDBChain *NewFrag = new PDBChain;
	*NewFrag = Frag;
	NewFrag->SetOrigin(a.x, a.y, a.z);

	uint OvPos = FindOverlap(*NewFrag);
	if (OvPos != UINT_MAX)
		{
		if (1)
			{
			coords_t c;
			NewFrag->GetCoords(OvPos, c);

			intpt_t Cube;
			GetCube(c, Cube);

			Log("\n");
			LogMe();
			Log("AppendCoords %.1f, %.1f, %.1f\n", a.x, a.y, a.z);
			Log("OvPos = %u", OvPos);
			Log(" x=%.1f, y=%.1f, z=%.1f", c.x, c.y, c.z);
			Log(" cube=%d,%d,%d", Cube.x, Cube.y, Cube.z);
			Log("\n");
			}
		delete NewFrag;
		Log("Found overlap\n");
		return false;
		}
	AppendFrag(*NewFrag, a);
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

bool FakeChain::AppendBest(const vector<PDBChain *> &Frags, uint Iters)
	{
	const uint FragCount = SIZE(Frags);
	if (m_Chain.GetSeqLength() == 0)
		{
		uint FragIdx = randu32()%FragCount;
		const PDBChain &Frag = *Frags[FragIdx];
		coords_t c0;
		c0.x = 0;
		c0.y = 0;
		c0.z = 0;
		AppendFrag(Frag, c0);
		return true;
		}

	uint BestFragIdx = UINT_MAX;
	double BestDiameter = DBL_MAX;
	double best_alpha = DBL_MAX;
	double best_beta = DBL_MAX;
	double best_gamma = DBL_MAX;
	for (uint Iter = 0; Iter < Iters; ++Iter)
		{
		FakeChain FC2 = *this;
		Ps(FC2.m_Chain.m_Label, "FC2.%u", Iter);
		uint FragIdx = randu32()%FragCount;
		const PDBChain &Frag = *Frags[FragIdx];
		PDBChain *NewFrag = new PDBChain;
		*NewFrag = Frag;
		double alpha = get_rand_radians();
		double beta = get_rand_radians();
		double gamma = get_rand_radians();
		RotateChain(*NewFrag, alpha, beta, gamma);
		bool Ok = FC2.TryAppendFrag(*NewFrag);
		delete NewFrag;
		if (!Ok)
			continue;

		double Diameter = m_Chain.GetDiameter();
		if (Diameter < BestDiameter)
			{
			BestDiameter = Diameter;
			BestFragIdx = FragIdx;
			best_alpha = alpha;
			best_beta = beta;
			best_gamma = gamma;
			}
		}

	if (BestFragIdx == UINT_MAX)
		return false;

	const PDBChain &Frag = *Frags[BestFragIdx];
	PDBChain *NewFrag = new PDBChain;
	*NewFrag = Frag;
	RotateChain(*NewFrag, best_alpha, best_beta, best_gamma);
	bool Ok = TryAppendFrag(*NewFrag);
	asserta(Ok);
	return true;
	}

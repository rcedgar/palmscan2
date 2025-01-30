#include "myutils.h"
#include "pdbchain.h"
#include "fakechain.h"

void GetRandomAnglePair_Radians(double &rad_bc, double &rad_vc);
double get_norm(const coords_t a);

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

void FakeChain::GetNEN_Plus(coords_t c, uint &Pos, double &Dist) const
	{
	Pos = UINT_MAX;
	Dist = DBL_MAX;
	const uint L = m_Chain.GetSeqLength();
	for (uint Pos2 = Pos + 4; Pos2 < L; ++Pos2)
		{
		coords_t cPos = GetCoords(Pos2);
		double d = get_dist(c, cPos);
		if (d < Dist)
			{
			Dist = d;
			Pos = Pos2;
			}
		}
	}

void FakeChain::GetNEN_Minus(coords_t c, uint &Pos, double &Dist) const
	{
	Pos = UINT_MAX;
	Dist = DBL_MAX;
	const uint L = m_Chain.GetSeqLength();
	for (int Pos2 = int(Pos) - 4; Pos2 >= 0; --Pos2)
		{
		coords_t cPos = GetCoords(uint(Pos2));
		double d = get_dist(c, cPos);
		if (d < Dist)
			{
			Dist = d;
			Pos = Pos2;
			}
		}
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
		coords_t c = m_Chain.GetCoords(Pos);

		Log("%5u |%c|  X=%5.1f  Y=%5.1f  Z=%5.1f\n",
			Pos, Seq[Pos], c.x, c.y, c.z);
		}
	}

bool FakeChain::IsOccupied(coords_t c, uint &Pos) const
	{
	Pos = UINT_MAX;
	const uint L = m_Chain.GetSeqLength();
	for (uint Pos2 = 0; Pos2 < L; ++Pos2)
		{
		coords_t cPos2 = GetCoords(Pos2);
		double d = get_dist(c, cPos2);
		if (d < m_MinNENDist)
			{
			Pos = Pos2;
			return true;
			}
		}
	return false;
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
	coords_t cterm1 = m_Chain.GetCoords(L-2);
	coords_t cterm = m_Chain.GetCoords(L-1);
	for (int Try = 0; Try < 100; ++Try)
		{
		double theta_bc, theta_vc;
		GetRandomAnglePair_Radians(theta_bc, theta_vc);
		Die("TODO");
		//Coords = NextCoords(cterm1, cterm, theta_bc, theta_vc);
		uint OvPos;
		if (!IsOccupied(Coords, OvPos))
			return true;
		}
	return false;
	}

uint FakeChain::FindOverlap(const PDBChain &Chain) const
	{
	const uint L = Chain.GetSeqLength();
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		coords_t Coords = Chain.GetCoords(Pos);
		uint OvPos;
		if (IsOccupied(Coords, OvPos))
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
			coords_t c = NewFrag->GetCoords(OvPos);

			Log("\n");
			LogMe();
			Log("AppendCoords %.1f, %.1f, %.1f\n", a.x, a.y, a.z);
			Log("OvPos = %u", OvPos);
			Log(" x=%.1f, y=%.1f, z=%.1f", c.x, c.y, c.z);
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

	for (uint FragPos = 0; FragPos < L; ++FragPos)
		{
		coords_t Coords = Frag.GetCoords(FragPos);
		uint OvPos;
		asserta(!IsOccupied(Coords, OvPos));

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
		coords_t c = m_Chain.GetCoords(Pos);

		uint FragIdx = m_PosToFragIdx[Pos];
		uint FragPos = m_PosToFragPos[Pos];
		asserta(FragIdx < FragCount);
		const PDBChain &Frag = *m_Frags[FragIdx];
		const uint FragL = Frag.GetSeqLength();
		asserta(FragPos < FragL);
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

coords_t FakeChain::NextCoords(
	const coords_t cterm2, const coords_t cterm1, const coords_t cterm,
	double theta_bc, double theta_vc) const
	{
	coords_t a, b;
	a.x = cterm1.x - cterm2.x;
	a.y = cterm1.y - cterm2.y;
	a.z = cterm1.z - cterm2.z;

	b.x = cterm.x - cterm1.x;
	b.y = cterm.y - cterm1.y;
	b.z = cterm.z - cterm1.z;

#if DEBUG
	{
	double moda = get_norm(a);
	double modb = get_norm(b);
	asserta(moda > 3.5 && moda < 4);
	asserta(modb > 3.5 && modb < 4);
	}
#endif

	coords_t nextc;
	nextc.x = cterm.x + (b.x - a.x)*m_CADist*sin(theta_bc);
	nextc.y = cterm.y + (b.y - a.y)*m_CADist*cos(theta_bc);
	Die("TODO");
	return a;
	}

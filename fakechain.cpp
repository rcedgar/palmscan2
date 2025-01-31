#include "myutils.h"
#include "pdbchain.h"
#include "fakechain.h"

void GetRandomAnglePair_Radians(double &rad_bc, double &rad_vc);
double get_norm(const coords_t a);
void LogCoords(const char *Name, coords_t c);

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

void FakeChain::GetNEN_Plus(uint Pos, uint &NENPos, double &Dist) const
	{
	NENPos = UINT_MAX;
	Dist = DBL_MAX;
	coords_t cPos = GetCoords(Pos);
	const uint L = m_Chain.GetSeqLength();
	for (uint Pos2 = Pos + 4; Pos2 < L; ++Pos2)
		{
		coords_t cPos2 = GetCoords(Pos2);
		double d = get_dist(cPos, cPos2);
		if (d < Dist)
			{
			Dist = d;
			NENPos = Pos2;
			}
		}
	}

void FakeChain::GetNEN_Minus(uint Pos, uint &NENPos, double &Dist) const
	{
	NENPos = UINT_MAX;
	Dist = DBL_MAX;
	coords_t cPos = GetCoords(Pos);
	const uint L = m_Chain.GetSeqLength();
	for (int Pos2 = int(Pos) - 4; Pos2 >= 0; --Pos2)
		{
		coords_t cPos2 = GetCoords(Pos2);
		double d = get_dist(cPos, cPos2);
		if (d < Dist)
			{
			Dist = d;
			NENPos = Pos2;
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

		Log("%5u |%c|  X=%5.1f  Y=%5.1f  Z=%5.1f",
			Pos, Seq[Pos], c.x, c.y, c.z);

		if (Pos > 0)
			{
			double dp, dm;
			uint NEN_Plus_Pos, NEN_Minus_Pos;
			GetNEN_Plus(Pos, NEN_Plus_Pos, dp);
			if (NEN_Plus_Pos != UINT_MAX)
				Log(" NEN+%4u=%4.1f", NEN_Plus_Pos, dp);

			GetNEN_Minus(Pos, NEN_Minus_Pos, dm);
			if (NEN_Minus_Pos != UINT_MAX)
				Log("   NEN-%4u=%4.1f", NEN_Minus_Pos, dm);

			double d = m_Chain.GetDist(Pos, Pos-1);
			if (d < 3.7 || d > 3.9)
				Log(" **CA<%.1f>", d);
			}

		Log("\n");
		}
	}

uint FakeChain::FindCollision(coords_t c, uint Lo, uint Hi) const
	{
	const uint L = m_Chain.GetSeqLength();
	for (uint Pos = Lo; Pos < Hi; ++Pos)
		{
		coords_t c2 = GetCoords(Pos);
		double d = get_dist(c, c2);
		if (d < m_MinNENDist)
			return Pos;
		}
	return UINT_MAX;
	}

coords_t get_unit_cd(coords_t A, coords_t B, coords_t C,
			   double theta_rad, double phi_rad);

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
	asserta(L >= 3);
	coords_t A = m_Chain.GetCoords(L-3);
	coords_t B = m_Chain.GetCoords(L-2);
	coords_t C = m_Chain.GetCoords(L-1);
	for (int Try = 0; Try < 100; ++Try)
		{
		double theta_rad, phi_rad;
		GetRandomAnglePair_Radians(theta_rad, phi_rad);
		
		coords_t unit_cd = get_unit_cd(A, B, C, theta_rad, phi_rad);
		
		coords_t cd;
		cd.x = unit_cd.x*m_CADist;
		cd.y = unit_cd.y*m_CADist;
		cd.z = unit_cd.z*m_CADist;
		
		Coords.x = C.x + cd.x;
		Coords.y = C.y + cd.y;
		Coords.z = C.z + cd.z;

#if DEBUG
		{
		double dx = C.x - Coords.x;
		double dy = C.y - Coords.y;
		double dz = C.z - Coords.z;
		double d = sqrt(dx*dx + dy*dy + dz*dz);
		assert(d > 3.7 && d < 3.9);
		}
#endif

		uint CollPos = FindCollision(Coords, 0, L-2);
		if (CollPos == UINT_MAX)
			return true;
		if (0)
			{
			LogMe();
			LogCoords("A", A);
			LogCoords("B", B);
			LogCoords("C", C);
			LogCoords("D", Coords);

			coords_t ov = GetCoords(CollPos);
			Log("OvPos %u\n", CollPos);
			LogCoords("ov", ov);
			}
		}
	return false;
	}

bool FakeChain::FitOk(const PDBChain &Frag,
					  uint &CollisionFakePos,
					  uint &CollisionFragPos) const
	{
	CollisionFakePos = UINT_MAX;
	CollisionFragPos = UINT_MAX;

	const uint FakeL = m_Chain.GetSeqLength();
	const uint FragL = Frag.GetSeqLength();
	for (uint FragPos = 0; FragPos < FragL; ++FragPos)
		{
		coords_t FragCoords = Frag.GetCoords(FragPos);
		for (uint FakePos = 0; FakePos < FakeL; ++FakePos)
			{
			coords_t FakeCoords = GetCoords(FakePos);
			double d = get_dist(FragCoords, FakeCoords);
			if (FragPos == 0 && FakePos + 1 == FakeL)
				{
				if (d < 3.5 || d > 4)
					{
					Warning("FitOk append");
					return false;
					}
				continue;
				}
			if (d < m_MinNENDist)
				{
				CollisionFakePos = FakePos;
				CollisionFragPos = FragPos;
				if (0)
					{
					Log("FitOk failed\n");
					Log("Fake pos %u coords ", FakePos);
					LogCoords("", GetCoords(FakePos));
					Log(" Frag pos %u coords", FragPos);
					LogCoords("", Frag.GetCoords(FragPos));
					Log("\n");
					}
				return false;
				}
			}
		}
	return true;
	}

const PDBChain *FakeChain::CreateFrag(uint LibIdx,
							const coords_t &AppendCoords,
							double alpha,
							double beta,
							double gamma) const
	{
	asserta(LibIdx < SIZE(m_Library));
	PDBChain *NewFrag = new PDBChain;
	*NewFrag = *m_Library[LibIdx];
	RotateChain(*NewFrag, alpha, beta, gamma);
	const uint FragL = NewFrag->GetSeqLength();
	asserta(FragL > 4);
	asserta(NewFrag->m_Xs[0] == 0);
	asserta(NewFrag->m_Ys[0] == 0);
	asserta(NewFrag->m_Zs[0] == 0);
	for (uint i = 0; i < FragL; ++i)
		{
		NewFrag->m_Xs[i] += AppendCoords.x;
		NewFrag->m_Ys[i] += AppendCoords.y;
		NewFrag->m_Zs[i] += AppendCoords.z;
		}
	return NewFrag;
	}

void FakeChain::AppendFrag(uint LibIdx,
						   const coords_t &AppendCoords,
						   double alpha, double beta, double gamma)
	{
	const PDBChain *NewFrag = CreateFrag(LibIdx, AppendCoords,
										  alpha, beta, gamma);
	const uint FragL = NewFrag->GetSeqLength();
	m_Frags.push_back(NewFrag);
	m_LibIdxs.push_back(LibIdx);
	m_AppendCoordsVec.push_back(AppendCoords);
	m_Alphas.push_back(alpha);
	m_Betas.push_back(beta);
	m_Gammas.push_back(gamma);

	for (uint FragPos = 0; FragPos < FragL; ++FragPos)
		{
		coords_t Coords = NewFrag->GetCoords(FragPos);
		m_Chain.m_Xs.push_back(Coords.x);
		m_Chain.m_Ys.push_back(Coords.y);
		m_Chain.m_Zs.push_back(Coords.z);
		}
	m_Chain.m_Seq += NewFrag->m_Seq;
	}

void FakeChain::Validate() const
	{
	Die("TODO");
	}

double FakeChain::GetQualityScore(const PDBChain &Chain) const
	{
	return Chain.GetDiameter();
	}

double FakeChain::GetQualityScoreFrag(const PDBChain &Frag) const
	{
	PDBChain Cat = m_Chain;
	Cat.m_Xs.insert(Cat.m_Xs.end(), Frag.m_Xs.begin(), Frag.m_Xs.end());
	Cat.m_Ys.insert(Cat.m_Ys.end(), Frag.m_Ys.begin(), Frag.m_Ys.end());
	Cat.m_Zs.insert(Cat.m_Zs.end(), Frag.m_Zs.begin(), Frag.m_Zs.end());
	return GetQualityScore(Cat);
	}

bool FakeChain::BestFit(uint Iters,
				 uint &BestLibIdx,
				 double &BestAlpha,
				 double &BestBeta,
				 double &BestGamma,
				 coords_t &BestAppendCoords) const
	{
	BestLibIdx = UINT_MAX;
	BestAlpha = DBL_MAX;
	BestGamma = DBL_MAX;
	BestBeta = DBL_MAX;

	const uint LibSize = SIZE(m_Library);
	asserta(m_Chain.GetSeqLength() > 0);

	double BestDiameter = DBL_MAX;
	for (uint Iter = 0; Iter < Iters; ++Iter)
		{
		uint LibIdx = randu32()%LibSize;
		coords_t AppendCoords;
		bool Ok = GetAppendCoords(AppendCoords);
		if (!Ok)
			return false;
		double alpha = get_rand_radians();
		double beta = get_rand_radians();
		double gamma = get_rand_radians();

		const PDBChain *NewFrag =
			CreateFrag(LibIdx, AppendCoords, alpha, beta, gamma);

		uint CollFakePos, CollFragPos;
		Ok = FitOk(*NewFrag, CollFakePos, CollFragPos);
		if (Ok)
			{
			double Diameter = GetQualityScoreFrag(*NewFrag);
			if (Diameter < BestDiameter)
				{
				BestDiameter = Diameter;
				BestLibIdx = LibIdx;
				BestAlpha = alpha;
				BestBeta = beta;
				BestGamma = gamma;
				BestAppendCoords = AppendCoords;
				}
			}
		delete NewFrag;
		}

	if (BestLibIdx == UINT_MAX)
		return false;
	return true;
	}

void FakeChain::Init()
	{
	m_Chain.Clear();
	m_Frags.clear();
	m_LibIdxs.clear();
	m_AppendCoordsVec.clear();
	m_Alphas.clear();
	m_Betas.clear();
	m_Gammas.clear();

	const uint LibSize = SIZE(m_Library);
	uint LibIdx = randu32()%LibSize;
	coords_t ZeroCoords;
	AppendFrag(LibIdx, ZeroCoords, 0, 0, 0);
	}

bool FakeChain::MakeFake(uint L)
	{
	Init();

	uint LibIdx;
	double Alpha;
	double Beta;
	double Gamma;
	coords_t AppendCoords;

	for (uint Try = 0; Try < 100; ++Try)
		{

		bool Ok = BestFit(100, LibIdx, Alpha, Beta, Gamma, AppendCoords);
		if (!Ok)
			break;
		if (m_Chain.GetSeqLength() >= L)
			return true;

		AppendFrag(LibIdx, AppendCoords, Alpha, Beta, Gamma);
		}

	return false;
	}

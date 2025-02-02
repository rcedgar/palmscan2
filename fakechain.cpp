#include "myutils.h"
#include "pdbchain.h"
#include "fakechain.h"

void GetRandomAnglePair_Radians(double &rad_bc, double &rad_vc);
double get_norm(const coords_t a);
void LogCoords(const char *Name, coords_t c);
double GetMDL(const PDBChain &Q);
double GetNENMed(const PDBChain &Q);

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

FakeChain::~FakeChain()
	{
	DeleteFrags();
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
	uint FragIdx = 0;
	uint FragPos = 0;
	const PDBChain *CurrFrag = m_Frags[0];
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		coords_t c = m_Chain.GetCoords(Pos);

		Log("%5u |%c|  X=%5.1f  Y=%5.1f  Z=%5.1f",
			Pos, Seq[Pos], c.x, c.y, c.z);

		Log(" {%3d,%2d}", FragIdx, FragPos);
		Log(" %.1f,%.1f,%.1f",
			CurrFrag->m_Xs[FragPos],
			CurrFrag->m_Ys[FragPos],
			CurrFrag->m_Zs[FragPos]);
		++FragPos;
		if (FragPos >= CurrFrag->GetSeqLength())
			{
			++FragIdx;
			FragPos = 0;
			if (Pos+1 < L)
				{
				asserta(FragIdx < SIZE(m_Frags));
				CurrFrag = m_Frags[FragIdx];
				}
			}

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

double FakeChain::FitOk(const PDBChain &Frag,
					  uint &CollisionFakePos,
					  uint &CollisionFragPos) const
	{
	CollisionFakePos = UINT_MAX;
	CollisionFragPos = UINT_MAX;

	const uint FakeL = m_Chain.GetSeqLength();
	const uint FragL = Frag.GetSeqLength();
	double MaxDist = 0;
	for (uint FragPos = 0; FragPos < FragL; ++FragPos)
		{
		coords_t FragCoords = Frag.GetCoords(FragPos);
		for (uint FakePos = 0; FakePos < FakeL; ++FakePos)
			{
			coords_t FakeCoords = GetCoords(FakePos);
			double d = get_dist(FragCoords, FakeCoords);
			MaxDist = max(d, MaxDist);
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
				return DBL_MAX;
				}
			}
		}
	return MaxDist;
	}

const PDBChain *FakeChain::CreateFrag(uint LibIdx,
							const coords_t &AppendCoords,
							double alpha,
							double beta,
							double gamma) const
	{
	asserta(LibIdx < SIZE(*m_Library));
	PDBChain *NewFrag = new PDBChain;
	*NewFrag = *(*m_Library)[LibIdx];
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
	const uint N = SIZE(m_LibIdxs);
	if (N == 0)
		return;

	FakeChain FC;
	FC.m_Library = m_Library;
	uint LibIdx0 = m_LibIdxs[0];
	FC.Init(LibIdx0);
	for (uint i = 1; i < N; ++i)
		{
		FC.AppendFrag(m_LibIdxs[i],
					  m_AppendCoordsVec[i],
					  m_Alphas[i],
					  m_Betas[i],
					  m_Gammas[i]);
		}
	const uint L = m_Chain.GetSeqLength();
	asserta(FC.m_Chain.GetSeqLength() == L);
	for (uint i = 0; i < L; ++i)
		{
		asserta(feq(m_Chain.m_Xs[i], FC.m_Chain.m_Xs[i]));
		asserta(feq(m_Chain.m_Ys[i], FC.m_Chain.m_Ys[i]));
		asserta(feq(m_Chain.m_Zs[i], FC.m_Chain.m_Zs[i]));
		}
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

	const uint LibSize = SIZE(*m_Library);
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
		double MaxDist = FitOk(*NewFrag, CollFakePos, CollFragPos);
		if (MaxDist != DBL_MAX)
			{
			//double Diameter = GetQualityScoreFrag(*NewFrag);
			//if (Diameter < BestDiameter)
			if (MaxDist < BestDiameter)
				{
				//BestDiameter = Diameter;
				BestDiameter = MaxDist;
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

void FakeChain::DeleteFrags()
	{
	const uint n = SIZE(m_Frags);
	for (uint i = 0; i < n; ++i)
		delete m_Frags[i];
	}

void FakeChain::Init(uint LibIdx)
	{
	DeleteFrags();
	m_Chain.Clear();
	m_Frags.clear();
	m_LibIdxs.clear();
	m_AppendCoordsVec.clear();
	m_Alphas.clear();
	m_Betas.clear();
	m_Gammas.clear();
	m_MDL = DBL_MAX;

	const uint LibSize = SIZE(*m_Library);
	coords_t ZeroCoords;
	AppendFrag(LibIdx, ZeroCoords, 0, 0, 0);
	}

bool FakeChain::MakeFake(uint L)
	{
	uint LibSize = SIZE(*m_Library);
	asserta(LibSize > 0);
	Init(randu32()%LibSize);

	uint LibIdx;
	double Alpha;
	double Beta;
	double Gamma;
	coords_t AppendCoords;

	for (uint Try = 0; Try < 100; ++Try)
		{
		bool Ok = BestFit(20, LibIdx, Alpha, Beta, Gamma, AppendCoords);
		if (!Ok)
			break;
		if (m_Chain.GetSeqLength() >= L)
			{
			m_MDL = GetMDL(m_Chain);
			m_NENMed = GetNENMed(m_Chain);
			return true;
			}

		AppendFrag(LibIdx, AppendCoords, Alpha, Beta, Gamma);
		}

	return false;
	}

void FakeChain::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, "%s", m_Chain.m_Label.c_str());
	const uint N = SIZE(m_LibIdxs);
	double NENMed = GetNENMed(m_Chain);
	fprintf(f, "\t%.3f", m_MDL);
	fprintf(f, "\t%.3f", NENMed);
	fprintf(f, "\t%u", N);
	for (uint i = 0; i < N; ++i)
		{
		uint LibIdx = m_LibIdxs[i];
		fprintf(f, "\t%s", (*m_Library)[LibIdx]->m_Label.c_str());
		const coords_t ac = m_AppendCoordsVec[i];
		fprintf(f, "\t%.1f", ac.x);
		fprintf(f, "\t%.1f", ac.y);
		fprintf(f, "\t%.1f", ac.z);
		fprintf(f, "\t%.4f", m_Alphas[i]);
		fprintf(f, "\t%.4f", m_Betas[i]);
		fprintf(f, "\t%.4f", m_Gammas[i]);
		}
	fprintf(f, "\n");
	}

// To color by chain at pymol command line or pml script:
// util.cbc()
void FakeChain::LoadFrag(const string &LoadDir, uint FragIdx, PDBChain &Frag) const
	{
	static const char ChainChars[] =
		"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwzyz";
	const uint NrChainChars = sizeof(ChainChars);
	asserta(FragIdx < SIZE(m_Frags));
	const char ChainChar = ChainChars[FragIdx%NrChainChars];
	const string &Label = m_Frags[FragIdx]->m_Label;
	vector<string> Fields;
	Split(Label, Fields, '|');
//        0       1   2   3       4
// frag6023|d1ijba_|107|130|H2S8L14
	asserta(SIZE(Fields) == 5);
	const string &Name = Fields[1];
	uint Lo = StrToUint(Fields[2]);
	uint Hi = StrToUint(Fields[3]);
	string FN = LoadDir + Name;
	PDBChain FullChain;
	vector<string> Lines;
	ReadLinesFromFile(FN, Lines);
	FullChain.FromPDBLines(Name, Lines);
	const uint FullLength = FullChain.GetSeqLength();
	asserta(SIZE(FullChain.m_ATOMs) == FullLength);
	asserta(Lo < Hi);
	asserta(Hi < FullLength);
	Frag.Clear();
	Frag.m_Label = Label;
	const uint FragL = Hi - Lo + 1;
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		char c = FullChain.m_Seq[Pos];
		double X = FullChain.m_Xs[Pos];
		double Y = FullChain.m_Ys[Pos];
		double Z = FullChain.m_Zs[Pos];

		Frag.m_Seq.push_back(c);
		Frag.m_Xs.push_back(X);
		Frag.m_Ys.push_back(Y);
		Frag.m_Zs.push_back(Z);

		Frag.m_ATOMs.push_back(FullChain.m_ATOMs[Pos]);
		}
	Frag.ZeroOrigin();
	const coords_t AppendCoords = m_AppendCoordsVec[FragIdx];
	RotateChain(Frag, m_Alphas[FragIdx], m_Betas[FragIdx], m_Gammas[FragIdx]);
	asserta(Frag.GetSeqLength() == FragL);
	for (uint FragPos = 0; FragPos < FragL; ++FragPos)
		{
		Frag.m_Xs[FragPos] += AppendCoords.x;
		Frag.m_Ys[FragPos] += AppendCoords.y;
		Frag.m_Zs[FragPos] += AppendCoords.z;

		vector<string> ResATOMs;
		const uint n = SIZE(Frag.m_ATOMs[FragPos]);
		for (uint i = 0; i < n; ++i)
			{
			double x, y, z;
			const string &Line = Frag.m_ATOMs[FragPos][i];
			PDBChain::GetXYZFromATOMLine(Line, x, y, z);
			x += AppendCoords.x;
			y += AppendCoords.y;
			z += AppendCoords.z;
			string UpdatedLine;
			PDBChain::SetXYZInATOMLine(Line, x, y, z, UpdatedLine);
			UpdatedLine[21] = ChainChar;
			ResATOMs.push_back(UpdatedLine);
			}
		Frag.m_ATOMs[FragPos] = ResATOMs;
		}
	}

void FakeChain::BuildPDB(const string &LoadDir,
						 PDBChain &Chain, PDBChain &ChainShuffledCA) const
	{
	Chain.Clear();
	ChainShuffledCA.Clear();

	Chain.m_Label = m_Chain.m_Label;
	ChainShuffledCA.m_Label = m_Chain.m_Label;

	const uint FragCount = SIZE(m_LibIdxs);
	for (uint FragIdx = 0; FragIdx < FragCount; ++FragIdx)
		{
		PDBChain Frag;
		LoadFrag(LoadDir, FragIdx, Frag);

		const uint FragL = Frag.GetSeqLength();
		asserta(SIZE(Frag.m_Seq) == FragL);
		asserta(SIZE(Frag.m_Xs) == FragL);
		asserta(SIZE(Frag.m_Ys) == FragL);
		asserta(SIZE(Frag.m_Zs) == FragL);
		asserta(SIZE(Frag.m_ATOMs) == FragL);

		string FragSeq = Frag.m_Seq;
		Chain.m_Seq += FragSeq;

		random_shuffle(FragSeq.begin(), FragSeq.end());
		ChainShuffledCA.m_Seq += FragSeq;
		
		for (uint i = 0; i < FragL; ++i)
			{
			Chain.m_Xs.push_back(Frag.m_Xs[i]);
			Chain.m_Ys.push_back(Frag.m_Ys[i]);
			Chain.m_Zs.push_back(Frag.m_Zs[i]);
			Chain.m_ATOMs.push_back(Frag.m_ATOMs[i]);

			ChainShuffledCA.m_Xs.push_back(Frag.m_Xs[i]);
			ChainShuffledCA.m_Ys.push_back(Frag.m_Ys[i]);
			ChainShuffledCA.m_Zs.push_back(Frag.m_Zs[i]);
			}
		}

	Chain.RenumberResidues(1);
	}

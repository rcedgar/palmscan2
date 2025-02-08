#include "myutils.h"
#include "chainfaker.h"
#include "abcxyz.h"

void GetThreeFromOne(char aa, string &AAA);
void LogCoords(const char *Name, coords_t c);
coords_t subtract(const coords_t &a, const coords_t &b);
coords_t add(const coords_t &a, const coords_t &b);
double get_angle(const coords_t &a, const coords_t &b);
coords_t cross_product(const coords_t &a, const coords_t &b);
coords_t normalize(const coords_t &a);
double get_dist(const coords_t &a, const coords_t &b);
coords_t rotate_around_vector(const coords_t &v,
							  const coords_t &axis, double theta_rad);

// Fisher-Yates shuffle:
// To shuffle an array a of n elements (indices 0 .. n-1):
//  for i from n - 1 downto 1 do
//       j := random integer with 0 <= j <= i
//       exchange a[j] and a[i]
static void Shuffle(string &Seq, uint Lo, uint Hi)
	{
	const unsigned N = Hi - Lo + 1;
	asserta(N > 0);
	for (unsigned i = N - 1; i >= 1; --i)
		{
		unsigned j = randu32()%(i + 1);
		
		char ci = Seq[Lo+i];
		char cj = Seq[Lo+j];

		Seq[Lo+i] = cj;
		Seq[Lo+j] = ci;
		}
	}

void GetSCOPSFFromLabel(const string &Label, string &SF)
	{
	vector<string> Fields;
	Split(Label, Fields, '/');
	asserta(SIZE(Fields) >= 2);
	const string &Fam = Fields[1];

	vector<string> Fields2;
	Split(Fam, Fields2, '.');
	SF = Fields2[0] + "." + Fields2[1] + "." + Fields2[2];
	}

vector<PDBChain *> ChainFaker::m_SCOP40;
vector<vector<uint> > ChainFaker::m_SFIdxToChainIdxs;
vector<string> ChainFaker::m_SFs;
vector<uint> ChainFaker::m_ChainIdxToSFIdx;
map<string, uint> ChainFaker::m_SFToIdx;

void ChainFaker::ToPDB(const string &FN) const
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);

	const string &Seq = m_FakeChain->m_Seq;
	const vector<double> &Xs = m_FakeChain->m_Xs;
	const vector<double> &Ys = m_FakeChain->m_Ys;
	const vector<double> &Zs = m_FakeChain->m_Zs;

	const size_t L = Xs.size();
	asserta(Ys.size() == L);
	asserta(Zs.size() == L);
	asserta(Seq.size() == L);
	asserta(m_PosToInsertIdx.size() == L);
	const string &Template = m_SCOP40[m_RealChainIdx]->m_Label;

	fprintf(f, "TITLE     %s\n", m_FakeChain->m_Label.c_str());
	fprintf(f, "REMARK    ");
	fprintf(f, "TEMPLATE %s", Template.c_str());
	fprintf(f, "  REPLACED %.1f%%", GetPctReplaced());
	fprintf(f, "\n");

	const uint InsertCount = SIZE(m_InsertedChainIdxs);
	for (uint Idx = 0; Idx < InsertCount; ++Idx)
		{

		uint TakeoutLo = m_InsertedTakeoutLos[Idx];
		uint TakeoutHi = m_InsertedTakeoutHis[Idx];
		uint ChainIdx = m_InsertedChainIdxs[Idx];
		uint Lo = m_InsertedLos[Idx];
		uint Hi = m_InsertedHis[Idx];
		double theta_rad = m_InsertedTheta2_rads[Idx];
		double theta_deg = degrees(theta_rad);

		uint PL = Hi - Lo + 1;
		uint FL = TakeoutHi - TakeoutLo + 1;

		const string &Dom = m_SCOP40[ChainIdx]->m_Label;
		string Tmp;
		Ps(Tmp, "(%u-%u,%u)", TakeoutLo+1, TakeoutHi+1, FL);

		fprintf(f, "REMARK    ");
		fprintf(f, "REPLACE %3u  %16s   =>",
				Idx+1, Tmp.c_str());
		fprintf(f, " %s(%u-%u,%u)%.1f",
				Dom.c_str(),
				Lo+1, Hi+1, PL,
				theta_deg);
		fprintf(f, "\n");
		}

	char Chain = 'A';

	uint LastInsertIdx = UINT_MAX;
	for (uint i = 0; i < L; ++i)
		{
		char aa = Seq[i];
		string sAAA;
		GetThreeFromOne(aa, sAAA);
		const char *AAA = sAAA.c_str();
		uint InsertIdx = m_PosToInsertIdx[i];
		if (InsertIdx != LastInsertIdx)
			{
			LastInsertIdx = InsertIdx;
			if (i > 0)
				++Chain;
			}

		fprintf(f, "ATOM  ");			//  1 -  6        Record name   "ATOM  "
		fprintf(f, "%5u", i+1);			//  7 - 11        Integer       serial       Atom serial number.
		fprintf(f, " ");				// 12
		fprintf(f, " CA ");				// 13 - 16        Atom          name         Atom name.
		fprintf(f, " ");				// 17             Character     altLoc       Alternate location indicator.
		fprintf(f, "%3.3s", AAA);		// 18 - 20        Residue name  resName      Residue name.
		fprintf(f, " ");				// 21
		fprintf(f, "%c", Chain);		// 22             Character     chainID      Chain identifier.
		fprintf(f, "%4u", i+1);			// 23 - 26        Integer       resSeq       Residue sequence number.
		fprintf(f, " ");				// 27             AChar         iCode        Code for insertion of residues.
		fprintf(f, "   ");				// 28 - 30
		fprintf(f, "%8.3f", Xs[i]);	// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		fprintf(f, "%8.3f", Ys[i]);	// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		fprintf(f, "%8.3f", Zs[i]);	// 47 - 57        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		fprintf(f, "%6.2f", 1.0);		// 55 - 60        Real(6.2)     occupancy    Occupancy.
		fprintf(f, "%6.2f", 0.0);		// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		fprintf(f, "          ");		// 67 - 76
		fprintf(f, " C");				// 77 - 78        LString(2)    element      Element symbol, right-justified.
		fprintf(f, "  ");				// 79 - 80        LString(2)    charge       Charge on the atom.

		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}

void ChainFaker::GetDistMetrics(double &MinNENDist,
	double &MinCADist, double &MaxCADist) const
	{
	const uint FL = m_FakeChain->GetSeqLength();
	MinNENDist = DBL_MAX;
	MinCADist = DBL_MAX;
	MaxCADist = 0;
	for (uint i = 0; i < FL; ++i)
		{
		if (i > 0)
			{
			double cad = m_FakeChain->GetDist(i, i-1);
			MinCADist = min(cad, MinCADist);
			MaxCADist = max(cad, MaxCADist);
			}
		double NENDist;
		uint NEN = m_FakeChain->GetNEN(i, NENDist);
		MinNENDist = min(NENDist, MinNENDist);
		}
	}

void ChainFaker::LogFake() const
	{
	const uint FL = m_FakeChain->GetSeqLength();
	Log("LogFake() L=%u\n", FL);
	for (uint i = 0; i < FL; ++i)
		{
		char c = m_FakeChain->m_Seq[i];

		double NENd;
		uint NEN = m_FakeChain->GetNEN(i, NENd);
		uint InsertIdx = m_PosToInsertIdx[i];
		Log("%5u", i);
		Log(" |%c|", c);
		if (InsertIdx == 0)
			Log("   ");
		else
			Log(" >%u", InsertIdx);
		Log("  %8.1f", m_FakeChain->m_Xs[i]);
		Log("  %8.1f", m_FakeChain->m_Ys[i]);
		Log("  %8.1f", m_FakeChain->m_Zs[i]);
		Log("  <%u,%.3g>", NEN, NENd);

		if (i > 0)
			{
			double cad = m_FakeChain->GetDist(i, i-1);
			if (cad < 3.7 || cad > 3.9)
				Log(" **CAD=%.3g", cad);
			}

		Log("\n");
		}
	}

bool ChainFaker::SetTakeout()
	{
	m_TakeoutLo = UINT_MAX;
	m_TakeoutHi = UINT_MAX;
	m_PlugChainIdxs.clear();
	m_PlugLos.clear();
	m_PlugHis.clear();

	const uint FL = m_FakeChain->GetSeqLength();
	for (uint Try = 0; Try < 10; ++Try)
		{
		double zero_to_one = (randu32()%1001)/1000.0;
		uint Lo = randu32()%FL;
		uint FragL = m_MinTakeoutLen +
			uint(zero_to_one*(m_MaxTakeoutLen - m_MinTakeoutLen));
		uint Hi = Lo + FragL -1;
		if (Hi >= FL)
			continue;
		m_TakeoutLo = Lo;
		m_TakeoutHi = Hi;
		uint TL = m_TakeoutHi - m_TakeoutLo + 1;
		if (m_Trace)
			Log("	Takeout %u .. %u (%u)\n",
				m_TakeoutLo, m_TakeoutHi, TL);
		return true;
		}
	return false;
	}

void ChainFaker::ReadSCOP40(const string &FN)
	{
	ReadChains(FN, m_SCOP40);
	const uint ChainCount = SIZE(m_SCOP40);

	map<string, vector<PDBChain *> > SFToChains;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		PDBChain *Chain = m_SCOP40[ChainIdx];
		const string &Label = Chain->m_Label;
		string SF;
		GetSCOPSFFromLabel(Label, SF);
		uint SFIdx = UINT_MAX;
		if (m_SFToIdx.find(SF) == m_SFToIdx.end())
			{
			SFIdx = SIZE(m_SFs);
			m_SFToIdx[SF] = SFIdx;
			m_SFs.push_back(SF);

			vector<uint> Empty;
			m_SFIdxToChainIdxs.push_back(Empty);
			}
		else
			SFIdx = m_SFToIdx[SF];
		m_SFIdxToChainIdxs[SFIdx].push_back(ChainIdx);
		m_ChainIdxToSFIdx.push_back(SFIdx);
		}
	}

uint ChainFaker::GetSFIdx(uint ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_ChainIdxToSFIdx));
	uint SFIdx = m_ChainIdxToSFIdx[ChainIdx];
	return SFIdx;
	}

const vector<uint> &ChainFaker::GetChainIdxsBySFIdx(uint SFIdx) const
	{
	asserta(SFIdx < SIZE(m_SFIdxToChainIdxs));
	const vector<uint> &ChainIdxs = m_SFIdxToChainIdxs[SFIdx];
	return ChainIdxs;
	}

uint ChainFaker::GetSFSize(uint SFIdx) const
	{
	const vector<uint> &ChainIdxs = GetChainIdxsBySFIdx(SFIdx);
	uint Size = SIZE(ChainIdxs);
	return Size;
	}

const PDBChain &ChainFaker::GetRealChain(uint ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_SCOP40));
	return *m_SCOP40[ChainIdx];
	}

const char *ChainFaker::GetChainLabel(uint ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_SCOP40));
	return m_SCOP40[ChainIdx]->m_Label.c_str();
	}

const char *ChainFaker::GetSFName(uint SFIdx) const
	{
	asserta(SFIdx < SIZE(m_SFs));
	return m_SFs[SFIdx].c_str();
	}

const PDBChain &ChainFaker::GetChain(uint ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_SCOP40));
	return *m_SCOP40[ChainIdx];
	}

double ChainFaker::GetCandidateDiameter(uint Idx) const
	{
	asserta(Idx < SIZE(m_PlugChainIdxs));
	uint ChainIdx = m_PlugChainIdxs[Idx];
	uint Lo = m_PlugLos[Idx];
	uint Hi = m_PlugHis[Idx];
	const PDBChain &Chain = GetRealChain(ChainIdx);
	asserta(Lo < Hi && Hi < Chain.GetSeqLength());
	double diameter_p = Chain.GetDist(Lo, Hi);
	return diameter_p;
	}

void ChainFaker::ValidateCandidatePlugDiameter(uint Idx) const
	{
	double diameter_t = GetTakeoutDiameter();
	double diameter_p = GetCandidateDiameter(Idx);
	double diameter_diff = fabs(diameter_p - diameter_t);
	if (diameter_diff > m_MaxFitError)
		{
		uint ChainIdx = m_PlugChainIdxs[Idx];
		uint Lo = m_PlugLos[Idx];
		uint Hi = m_PlugHis[Idx];
		Log("\nChainIdx = %u, %u .. %u diameters=%.3g, %.3g\n",
			ChainIdx, Lo, Hi, diameter_t, diameter_p);
		Die("ValidateCandidatePlugDiameter");
		}
	}

void ChainFaker::ValidateCandidatePlugDiameters() const
	{
#if DEBUG
	const uint N = SIZE(m_PlugChainIdxs);
	asserta(SIZE(m_PlugLos) == N);
	asserta(SIZE(m_PlugHis) == N);
	double diameter_t = GetTakeoutDiameter();

	for (uint Idx = 0; Idx < N; ++Idx)
		{
		double diameter_p = GetCandidateDiameter(Idx);
		double diameter_diff = fabs(diameter_p - diameter_t);
		if (diameter_diff > m_MaxFitError)
			{
			uint ChainIdx = m_PlugChainIdxs[Idx];
			uint Lo = m_PlugLos[Idx];
			uint Hi = m_PlugHis[Idx];
			Log("\nChainIdx = %u, %u .. %u diameters=%.3g, %.3g\n",
				ChainIdx, Lo, Hi, diameter_t, diameter_p);
			Die("ValidateCandidatePlugDiameters");
			}
		}
#endif
	}

double ChainFaker::GetTakeoutDiameter() const
	{
	coords_t tlo, thi;
	tlo.x = m_FakeChain->m_Xs[m_TakeoutLo];
	tlo.y = m_FakeChain->m_Ys[m_TakeoutLo];
	tlo.z = m_FakeChain->m_Zs[m_TakeoutLo];

	thi.x = m_FakeChain->m_Xs[m_TakeoutHi];
	thi.y = m_FakeChain->m_Ys[m_TakeoutHi];
	thi.z = m_FakeChain->m_Zs[m_TakeoutHi];

	double diameter_t = get_dist(thi, tlo);
	return diameter_t;
	}

void ChainFaker::FindCandidatePlugs1(uint ChainIdx,
									 uint PL_lo, uint PL_hi)
	{
	if (ChainIdx == m_RealChainIdx)
		return;
	const double TakeoutTermDist = GetTakeoutDiameter();
	const PDBChain &Chain = GetRealChain(ChainIdx);
	const uint L = Chain.GetSeqLength();
	const uint TL = m_TakeoutHi - m_TakeoutLo + 1;
	for (uint Lo = 0; Lo + PL_lo <= L; ++Lo)
		{
		for (uint Hi = Lo + PL_lo-1; Hi <= Lo + PL_hi-1; ++Hi)
			{
			if (Hi >= L)
				break;
			const uint PL = Hi - Lo + 1;
			if (PL == TL)
				continue;
			double d = Chain.GetDist(Lo, Hi);
			double Error = fabs(d - TakeoutTermDist);
			if (Error <= m_MaxFitError)
				{
				uint Idx = SIZE(m_PlugChainIdxs);
				m_PlugChainIdxs.push_back(ChainIdx);
				m_PlugLos.push_back(Lo);
				m_PlugHis.push_back(Hi);

				if (0 && m_Trace)
					{
					Log(" candidate plug[%u]", Idx);
					Log("  %4u .. %4u", Lo, Hi);
					Log("  d=%.2f, err=%.2f", d, Error);
					Log("  >%s", GetChainLabel(ChainIdx));
					Log("\n");
					}
				}
			}
		}
	}

void ChainFaker::GetTakeoutCoords(coords_t &tlo, coords_t &thi,
								  coords_t &tcentroid) const
	{
	tlo = m_FakeChain->GetCoords(m_TakeoutLo);
	thi = m_FakeChain->GetCoords(m_TakeoutHi);
	tcentroid.x = (tlo.x + thi.x)/2;
	tcentroid.y = (tlo.y + thi.y)/2;
	tcentroid.z = (tlo.z + thi.z)/2;
	}

void ChainFaker::FindCandidatePlugs()
	{
	m_PlugChainIdxs.clear();
	m_PlugLos.clear();
	m_PlugHis.clear();

	asserta(m_TakeoutLo < m_TakeoutHi);
	uint TL = m_TakeoutHi - m_TakeoutLo + 1;
	uint PL_min = TL - 2;
	uint PL_max = TL + 4;

	const uint N = SIZE(*m_RealSFChainIdxs);
	for (uint i = 0; i < N; ++i)
		{
		uint RealChainIdx = (*m_RealSFChainIdxs)[i];
		FindCandidatePlugs1(RealChainIdx, PL_min, PL_max);
		}
	}

void ChainFaker::AlignCentroid(PDBChain &Plug) const
	{
	uint PL = Plug.GetSeqLength();
	asserta(SIZE(Plug.m_Xs) == PL);
	asserta(SIZE(Plug.m_Ys) == PL);
	asserta(SIZE(Plug.m_Zs) == PL);

	asserta(m_FakeChain != 0);

	coords_t plo, phi, centroid_p;
	GetPlugCoords(Plug, plo, phi, centroid_p);

	coords_t tlo, thi, centroid_t;
	GetTakeoutCoords(tlo, thi, centroid_t);

	double diameter_t = get_dist(thi, tlo);
	double diameter_p = get_dist(phi, plo);
	double diameter_diff = fabs(diameter_p - diameter_t);
	if (diameter_diff > m_MaxFitError)
		Warning("diameter_diff %.3g > m_MaxFitErr #1",
				diameter_diff);
	
	coords_t t;
	t.x = centroid_t.x - centroid_p.x;
	t.y = centroid_t.y - centroid_p.y;
	t.z = centroid_t.z - centroid_p.z;

// Translate Plug so that its centroid coincides
// with Takeout centroid.
	for (uint i = 0; i < PL; ++i)
		{
		Plug.m_Xs[i] += t.x;
		Plug.m_Ys[i] += t.y;
		Plug.m_Zs[i] += t.z;
		}
	}

void ChainFaker::GetRealSubchain(uint ChainIdx, uint Lo, uint Hi,
						 PDBChain &Plug) const
	{
	Plug.Clear();

	const PDBChain &Chain = GetChain(ChainIdx);
	Ps(Plug.m_Label, "%s[%u-%u]", Chain.m_Label.c_str(), Lo, Hi);

	for (uint i = Lo; i <= Hi; ++i)
		{
		char c = Chain.m_Seq[i];
		Plug.m_Seq.push_back(c);

		double X = Chain.m_Xs[i];
		double Y = Chain.m_Ys[i];
		double Z = Chain.m_Zs[i];

		Plug.m_Xs.push_back(X);
		Plug.m_Ys.push_back(Y);
		Plug.m_Zs.push_back(Z);
		}
	}

double ChainFaker::RotatePlug1(const PDBChain &Plug,
							 PDBChain &RotatedPlug) const
	{
	RotatedPlug.Clear();
	RotatedPlug.m_Label = Plug.m_Label;
	RotatedPlug.m_Seq = Plug.m_Seq;

	coords_t plo, phi, centroid_p;
	GetPlugCoords(Plug, plo, phi, centroid_p);

	coords_t tlo, thi, centroid_t;
	GetTakeoutCoords(tlo, thi, centroid_t);

// Find angle between vectors pv=PlugLo->Hi and tv=TakeoutLo->Hi
	coords_t pv = subtract(phi, plo);
	coords_t tv = subtract(thi, tlo);
	double theta1_rad = get_angle(tv, pv);

// Take cross product to get direction perpendicular to plane
// defined by pv and tv
	coords_t axis1 = normalize(cross_product(tv, pv));

// Rotate by angle theta around this direction, this brings
// PlugLo close to TakeoutLo and
// PlugHi close to TakeoutHi
	const uint PL = Plug.GetSeqLength();
	for (uint i = 0; i < PL; ++i)
		{
		coords_t pt = Plug.GetCoords(i);
		coords_t pt2 = subtract(pt, centroid_t);
		coords_t pt3 = rotate_around_vector(pt2, axis1, theta1_rad);
		coords_t pt4 = add(pt3, centroid_t);

		RotatedPlug.m_Xs.push_back(pt4.x);
		RotatedPlug.m_Ys.push_back(pt4.y);
		RotatedPlug.m_Zs.push_back(pt4.z);
		}

	coords_t final_plo, final_phi, final_centroid_p;
	GetPlugCoords(RotatedPlug, final_plo, final_phi, final_centroid_p);

	double d = get_dist(final_centroid_p, centroid_t);
	asserta(d < 1);

// These distances should be small
	double dlo_final = get_dist(tlo, final_plo);
	double dhi_final = get_dist(thi, final_phi);
	double dtot = dlo_final + dhi_final;
	if (dtot > 1)
		{
		Warning("dtot=%.3g", dtot);
		return DBL_MAX;
		}
	if (0 && m_Trace)
		Log("TryFitPlug1() dtot=%.3g\n", dtot);
	return dtot;
	}

void ChainFaker::GetPlugCoords(
	const PDBChain &Plug, coords_t &plo, coords_t &phi,
	coords_t &pcentroid) const
	{
	const uint PL = Plug.GetSeqLength();
	plo = Plug.GetCoords(0);
	phi = Plug.GetCoords(PL-1);
	pcentroid.x = (plo.x + phi.x)/2;
	pcentroid.y = (plo.y + phi.y)/2;
	pcentroid.z = (plo.z + phi.z)/2;
	}

double ChainFaker::TryFitPlug1(uint PlugIdx,
	coords_t &t,
	double &theta1_rad,
	PDBChain &Plug) const
	{
	uint ChainIdx = m_PlugChainIdxs[PlugIdx];
	uint Lo = m_PlugLos[PlugIdx];
	uint Hi = m_PlugHis[PlugIdx];
	const uint PL = Hi - Lo + 1;

	PDBChain Plug1;
	GetRealSubchain(ChainIdx, Lo, Hi, Plug1);
	AlignCentroid(Plug1);

	coords_t tlo, thi, centroid_t;
	GetTakeoutCoords(tlo, thi, centroid_t);

	coords_t translated_plo, translated_phi, translated_centroid_p;
	GetPlugCoords(Plug1, translated_plo, translated_phi,
				   translated_centroid_p);

// Validate that centroids coincide after translation
	asserta(feq(translated_centroid_p.x, centroid_t.x));
	asserta(feq(translated_centroid_p.y, centroid_t.y));
	asserta(feq(translated_centroid_p.z, centroid_t.z));

//// Find angle between vectors pv=PlugLo->Hi and tv=TakeoutLo->Hi
//	coords_t pv = subtract(translated_phi, translated_plo);
//	coords_t tv = subtract(thi, tlo);
//	theta1_rad = get_angle(tv, pv);
//
//// Take cross product to get direction perpendicular to plane
//// defined by pv and tv
//	coords_t axis1 = normalize(cross_product(tv, pv));
//
//// Rotate by angle theta around this direction, this brings
//// PlugLo close to TakeoutLo and
//// PlugHi close to TakeoutHi
//	for (uint i = 0; i < PL; ++i)
//		{
//		coords_t pt = Plug.GetCoords(i);
//		coords_t pt2 = subtract(pt, centroid_t);
//		coords_t pt3 = rotate_around_vector(pt2, axis1, theta1_rad);
//		coords_t pt4 = add(pt3, centroid_t);
//
//		Plug.m_Xs[i] = pt4.x;
//		Plug.m_Ys[i] = pt4.y;
//		Plug.m_Zs[i] = pt4.z;
//		}
//
//	coords_t final_plo, final_phi, final_centroid_p;
//	GetPlugCoords(Plug, final_plo, final_phi, final_centroid_p);
//
//	double d = get_dist(final_centroid_p, centroid_t);
//	asserta(d < 1);

//// These distances should be small
//	double dlo_final = get_dist(tlo, final_plo);
//	double dhi_final = get_dist(thi, final_phi);
//	double dtot = dlo_final + dhi_final;
//	if (dtot > 1)
//		{
//		Warning("dtot=%.3g", dtot);
//		return DBL_MAX;
//		}
//	if (0 && m_Trace)
//		Log("TryFitPlug1() dtot=%.3g\n", dtot);

	double dtot = RotatePlug1(Plug1, Plug);
	return dtot;
	}

double ChainFaker::FindBadNENDist(const PDBChain &Plug,
	uint &FakePos, uint &PlugPos) const
	{
	FakePos = UINT_MAX;
	PlugPos = UINT_MAX;
	const uint FL = m_FakeChain->GetSeqLength();
	const uint PL = Plug.GetSeqLength();
	for (uint i = 0; i < PL; ++i)
		{
		const coords_t ppt = Plug.GetCoords(i);
		for (uint j = 0; j < FL; ++j)
			{
			if (j >= m_TakeoutLo && j <= m_TakeoutHi)
				continue;
			if (i == 0 && j == m_TakeoutLo - 1)
				continue;
			if (i + 1 == PL && j == m_TakeoutHi + 1)
				continue;
			const coords_t fpt = m_FakeChain->GetCoords(j);
			double d = get_dist(ppt, fpt);
			if (d < m_MinNENDist)
				{
				PlugPos = i;
				FakePos = j;
				return d;
				}
			}
		}
	return DBL_MAX;
	}

void ChainFaker::RotatePlug2(const PDBChain &Plug, double theta2_rad,
	PDBChain &RotatedPlug) const
	{
	const uint PL = Plug.GetSeqLength();
	RotatedPlug.Clear();
	RotatedPlug.m_Seq = Plug.m_Seq;

	coords_t plo, phi, centroid_p;
	GetPlugCoords(Plug, plo, phi, centroid_p);

// Centroids of takeout and plug should be identical
// (withing rounding error)
// Verify centroids are the same
	coords_t tlo, thi, centroid_t;
	GetTakeoutCoords(tlo, thi, centroid_t);
	asserta(fabs(centroid_t.x - centroid_p.x) < 0.5);
	asserta(fabs(centroid_t.y - centroid_p.y) < 0.5);
	asserta(fabs(centroid_t.z - centroid_p.z) < 0.5);

// Axis of rotation is start-end of takeout (same as 
// start-end of plug)
	coords_t axis = normalize(subtract(phi, plo));

// Rotate by angle theta2 around this direction, this leaves
// coordinates of start, end and centroid unchanged
	for (uint i = 0; i < PL; ++i)
		{
		coords_t pt = Plug.GetCoords(i);
		coords_t pt2 = subtract(pt, centroid_p);
		coords_t pt3 = rotate_around_vector(pt2, axis, theta2_rad);
		coords_t pt4 = add(pt3, centroid_p);

		RotatedPlug.m_Xs.push_back(pt4.x);
		RotatedPlug.m_Ys.push_back(pt4.y);
		RotatedPlug.m_Zs.push_back(pt4.z);
		}

#if DEBUG
// Verify unchanged
	coords_t centroid_r; // _r for rotated
	{
	coords_t rlo, rhi, centroid_r;
	GetPlugCoords(RotatedPlug, rlo, rhi, centroid_r);

	asserta(feq(rlo.x, plo.x));
	asserta(feq(rlo.y, plo.y));
	asserta(feq(rlo.z, plo.z));

	asserta(feq(rhi.x, phi.x));
	asserta(feq(rhi.y, phi.y));
	asserta(feq(rhi.z, phi.z));

	asserta(fabs(centroid_r.x - centroid_p.x) < 0.5);
	asserta(fabs(centroid_r.y - centroid_p.y) < 0.5);
	asserta(fabs(centroid_r.z - centroid_p.z) < 0.5);
	}
#endif
	}

bool ChainFaker::TryFitPlug2(PDBChain &Plug,
	double &theta2_rad) const
	{
	const uint PL = Plug.GetSeqLength();
	asserta(SIZE(Plug.m_Xs) == PL);
	asserta(SIZE(Plug.m_Ys) == PL);
	asserta(SIZE(Plug.m_Zs) == PL);

	coords_t plo, phi;
	plo.x = Plug.m_Xs[0];
	plo.y = Plug.m_Ys[0];
	plo.z = Plug.m_Zs[0];

	phi.x = Plug.m_Xs[PL-1];
	phi.y = Plug.m_Ys[PL-1];
	phi.z = Plug.m_Zs[PL-1];

#if DEBUG
	{
	coords_t tlo, thi;
	tlo.x = m_FakeChain->m_Xs[m_TakeoutLo];
	tlo.y = m_FakeChain->m_Ys[m_TakeoutLo];
	tlo.z = m_FakeChain->m_Zs[m_TakeoutLo];

	thi.x = m_FakeChain->m_Xs[m_TakeoutHi];
	thi.y = m_FakeChain->m_Ys[m_TakeoutHi];
	thi.z = m_FakeChain->m_Zs[m_TakeoutHi];
	double dlo = get_dist(tlo, plo);
	double dhi = get_dist(thi, phi);
	double dtot = dlo + dhi;
	asserta(dtot < 1);
	}
#endif

	const uint TICKS = 10;
	const double Tick = 2*PI/TICKS;
	PDBChain RotatedPlug;
	uint Counter = 0;
	for (double theta_rad = 0; theta_rad < 2*PI; theta_rad += Tick)
		{
		RotatePlug2(Plug, theta_rad, RotatedPlug);

		uint BadFakePos = UINT_MAX;
		uint BadPlugPos = UINT_MAX;
		double d = FindBadNENDist(RotatedPlug, BadFakePos, BadPlugPos);
		if (d == DBL_MAX)
			{
			theta2_rad = theta_rad;
			Plug = RotatedPlug;
			return true;
			}
		}
	return false;
	}

void ChainFaker::ReplaceTakeoutWithPlug(const PDBChain &Plug)
	{
	const uint FL = m_FakeChain->GetSeqLength();
	const uint PL = Plug.GetSeqLength();
	++m_InsertIdx;

	vector<uint> NewPosToInsertIdx;
	PDBChain NewChain;
	for (uint i = 0; i < m_TakeoutLo; ++i)
		{
		NewChain.m_Seq += m_FakeChain->m_Seq[i];
		NewChain.m_Xs.push_back(m_FakeChain->m_Xs[i]);
		NewChain.m_Ys.push_back(m_FakeChain->m_Ys[i]);
		NewChain.m_Zs.push_back(m_FakeChain->m_Zs[i]);
		NewPosToInsertIdx.push_back(m_PosToInsertIdx[i]);
		}

	for (uint i = 0; i < PL; ++i)
		{
		NewChain.m_Seq += Plug.m_Seq[i];
		NewChain.m_Xs.push_back(Plug.m_Xs[i]);
		NewChain.m_Ys.push_back(Plug.m_Ys[i]);
		NewChain.m_Zs.push_back(Plug.m_Zs[i]);
		NewPosToInsertIdx.push_back(m_InsertIdx);
		}

	for (uint i = m_TakeoutHi+1; i < FL; ++i)
		{
		NewChain.m_Seq += m_FakeChain->m_Seq[i];
		NewChain.m_Xs.push_back(m_FakeChain->m_Xs[i]);
		NewChain.m_Ys.push_back(m_FakeChain->m_Ys[i]);
		NewChain.m_Zs.push_back(m_FakeChain->m_Zs[i]);
		NewPosToInsertIdx.push_back(m_PosToInsertIdx[i]);
		}
	*m_FakeChain = NewChain;
	m_PosToInsertIdx = NewPosToInsertIdx;
	}

bool ChainFaker::TryFitPlug(uint PlugIdx,
	coords_t &t,
	double &theta1_rad,
	double &theta2_rad,
	PDBChain &Plug) const
	{
	double TakeoutDiameter = GetTakeoutDiameter();
	double CandidateDiameter = GetCandidateDiameter(PlugIdx);
	double Diff = fabs(CandidateDiameter - TakeoutDiameter);
	asserta(Diff <= m_MaxFitError);
	double dtot = TryFitPlug1(PlugIdx, t, theta1_rad, Plug);
	if (dtot > 1)
		return false;
	bool Ok = TryFitPlug2(Plug, theta2_rad);
	return Ok;
	}

void ChainFaker::MakePlug(uint ChainIdx, uint Lo, uint Hi,
				double theta2_rad, PDBChain &Plug) const
	{
	Plug.Clear();

	coords_t tlo, thi, centroid_t;
	GetTakeoutCoords(tlo, thi, centroid_t);

	const uint PL = Hi - Lo + 1;
	PDBChain Plug1;
	GetRealSubchain(ChainIdx, Lo, Hi, Plug1);

	coords_t plo1, phi1, centroid_p1;
	GetPlugCoords(Plug1, plo1, phi1, centroid_p1);

	double diameter_t = get_dist(thi, tlo);
	double diameter_p1 = get_dist(phi1, plo1);
	double diameter_diff = fabs(diameter_p1 - diameter_t);
	if (diameter_diff > m_MaxFitError)
		Warning("diameter_diff %.3g-%.3g%.3g > m_MaxFitErr #2",
				diameter_t, diameter_p1, diameter_diff);

	AlignCentroid(Plug1);

// Validate that centroids coincide after translation
	coords_t translated_plo, translated_phi, translated_centroid_p;
	GetPlugCoords(Plug1, translated_plo, translated_phi,
				  translated_centroid_p);

	asserta(feq(translated_centroid_p.x, centroid_t.x));
	asserta(feq(translated_centroid_p.y, centroid_t.y));
	asserta(feq(translated_centroid_p.z, centroid_t.z));

	PDBChain Plug2;
	RotatePlug1(Plug1, Plug2);

// Validate that centroids coincide after first rotation
	coords_t rotated_plo, rotated_phi, rotated_centroid_p;
	GetPlugCoords(Plug1, rotated_plo, rotated_phi,
				  rotated_centroid_p);

	asserta(feq(rotated_centroid_p.x, centroid_t.x));
	asserta(feq(rotated_centroid_p.y, centroid_t.y));
	asserta(feq(rotated_centroid_p.z, centroid_t.z));

	RotatePlug2(Plug2, theta2_rad, Plug);
	}

void ChainFaker::AssertPlugsEq(const PDBChain &Plug1, const PDBChain &Plug2) const
	{
	const uint L = Plug1.GetSeqLength();
	asserta(Plug2.GetSeqLength() == L);
	asserta(Plug1.m_Seq == Plug2.m_Seq);
	for (uint i = 0; i < L; ++i)
		{
		asserta(feq(Plug1.m_Xs[i], Plug2.m_Xs[i]));
		asserta(feq(Plug1.m_Ys[i], Plug2.m_Ys[i]));
		asserta(feq(Plug1.m_Zs[i], Plug2.m_Zs[i]));
		}
	}

bool ChainFaker::Mutate()
	{
	bool Ok = SetTakeout();
	if (!Ok)
		{
		if (m_Trace)
			Log("	**FAIL** SetTakeout()\n");
		return false;
		}
	if (m_Trace)
		{
		double TakeoutDiameter = GetTakeoutDiameter();
		Log("	diameter_t = %.3g\n", TakeoutDiameter);
		}

	FindCandidatePlugs();
	ValidateCandidatePlugDiameters();
	const uint NrCandidatePlugs = SIZE(m_PlugChainIdxs);
	asserta(SIZE(m_PlugLos) == NrCandidatePlugs);
	asserta(SIZE(m_PlugHis) == NrCandidatePlugs);
	if (m_Trace)
		Log("	%u candidate plugs\n", NrCandidatePlugs);
	if (NrCandidatePlugs == 0)
		return false;

	vector<uint> PlugIdxs;
	for (uint i = 0; i < NrCandidatePlugs; ++i)
		PlugIdxs.push_back(i);
	Shuffle(PlugIdxs);

	uint FailedCount = 0;
	for (uint i = 0; i < NrCandidatePlugs; ++i)
		{
		double theta1_rad, theta2_rad;
		PDBChain Plug;
		uint PlugIdx = PlugIdxs[i];
		coords_t t;
		bool Ok = TryFitPlug(PlugIdx, t, theta1_rad, theta2_rad, Plug);
		if (!Ok)
			{
			++FailedCount;
			continue;
			}
		uint ChainIdx = m_PlugChainIdxs[PlugIdx];
		uint Lo = m_PlugLos[PlugIdx];
		uint Hi = m_PlugHis[PlugIdx];

#if DEBUG
		{
		PDBChain CheckPlug;
		MakePlug(ChainIdx, Lo, Hi, theta2_rad, CheckPlug);
		AssertPlugsEq(CheckPlug, Plug);
		}
#endif
		ReplaceTakeoutWithPlug(Plug);
		if (m_Trace)
			{
			m_InsertedChainIdxs.push_back(ChainIdx);
			m_InsertedTakeoutLos.push_back(m_TakeoutLo);
			m_InsertedTakeoutHis.push_back(m_TakeoutHi);
			m_InsertedLos.push_back(Lo);
			m_InsertedHis.push_back(Hi);
			m_InsertedTheta2_rads.push_back(theta2_rad);

			Log("   Inserted plug %u '%s', fails=%u\n",
				PlugIdx, Plug.m_Seq.c_str(), FailedCount);
			//LogFake();
			double MinNENDist, MinCADist, MaxCADist;
			GetDistMetrics(MinNENDist, MinCADist, MaxCADist);
			Log("	ninnen=%.3g cadist %.2f .. %.2f\n",
				MinNENDist, MinCADist, MaxCADist);
			}
		return true;
		}
	if (m_Trace)
		Log("	Mutate() failed\n");
	return false;
	}

void ChainFaker::ShuffleFakeSequence(uint Pct)
	{
	string &Seq = m_FakeChain->m_Seq;
	const uint L = SIZE(Seq);

// Each swap creates two new diffs
	double Fract = Pct/200.0;
	const uint N = uint(L*Fract+0.5);
	for (uint i = 0; i < N; ++i)
		{
		uint Pos1 = randu32()%L;
		uint Pos2 = randu32()%L;
		std::swap(Seq[Pos1], Seq[Pos2]);
		}
	}

void ChainFaker::GetSFStr(string &SF) const
	{
	asserta(m_RealSFIdx < SIZE(m_SFs));
	SF = m_SFs[m_RealSFIdx];
	}

double ChainFaker::GetPctReplaced() const
	{
	const uint L = m_FakeChain->GetSeqLength();
	asserta(SIZE(m_PosToInsertIdx) == L);
	uint n = 0;
	for (uint i = 0; i < L; ++i)
		if (m_PosToInsertIdx[i] != 0)
			++n;
	return n*100.0/L;
	}

void ChainFaker::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, "%s", m_FakeChain->m_Label.c_str());
	const string &Template = m_SCOP40[m_RealChainIdx]->m_Label;
	fprintf(f, "\t%s", Template.c_str());
	const uint InsertCount = SIZE(m_InsertedChainIdxs);
	double Pct = GetPctReplaced();
	fprintf(f, "\t%.1f", Pct);
	fprintf(f, "\t%u", InsertCount);
	for (uint Idx = 0; Idx < InsertCount; ++Idx)
		{
		uint TakeoutLo = m_InsertedTakeoutLos[Idx];
		uint TakeoutHi = m_InsertedTakeoutHis[Idx];
		uint ChainIdx = m_InsertedChainIdxs[Idx];
		uint Lo = m_InsertedLos[Idx];
		uint Hi = m_InsertedHis[Idx];
		double theta_rad = m_InsertedTheta2_rads[Idx];

		const string &Dom = m_SCOP40[ChainIdx]->m_Label;
		fprintf(f, "\t%u\t%u", TakeoutLo+1, TakeoutHi+1);
		fprintf(f, " \t%s\t%u\t%u\t%.1f", Dom.c_str(),
				Lo+1, Hi+1, theta_rad);
		}
	fprintf(f, "\n");
	}

bool ChainFaker::MakeFake(uint ChainIdx, PDBChain &FakeChain)
	{
	Reset();
	m_FakeChain = &FakeChain;
	m_RealChainIdx = ChainIdx;
	m_RealChain = &GetChain(m_RealChainIdx);
	m_RealSFIdx = GetSFIdx(ChainIdx);
	uint SFSize = GetSFSize(m_RealSFIdx);
	if (m_Trace)
		{
		Log("\n");
		Log("MakeFake(ChainIdx=%u)", m_RealChainIdx);
		Log(" fold=%s", GetSFName(m_RealSFIdx));
		Log(" size=%u", SFSize);
		Log("\n");
		}

	if (SFSize < 3)
		{
		if (m_Trace)
			Log("	**FAIL** size too small\n");
		return false;
		}

	m_RealSFChainIdxs =
		&GetChainIdxsBySFIdx(m_RealSFIdx);
	const PDBChain &RealChain = GetRealChain(m_RealChainIdx);
	*m_FakeChain = RealChain;
	const uint L = m_FakeChain->GetSeqLength();
	m_InsertIdx = 0;
	for (uint i = 0; i < L; ++i)
		m_PosToInsertIdx.push_back(0);

	double ReplacedPct = 0;
	for (uint k = 0; k < 100; ++k)
		{
		ReplacedPct = GetPctReplaced();
		if (ReplacedPct >= m_MinReplacedPct)
			break;
		Mutate();
		}

	if (m_Trace)
		Log("	%u inserts, replaced %.1f%%\n",
			SIZE(m_InsertedChainIdxs), ReplacedPct);
	if (ReplacedPct < m_MinReplacedPct)
		return false;
	ShuffleFakeSequence(m_ShufflePct);
	if (m_Reverse)
		{
		PDBChain ReversedChain;
		m_FakeChain->GetReverse(ReversedChain);
		*m_FakeChain = ReversedChain;
		}
	return true;
	}

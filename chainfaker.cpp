#include "myutils.h"
#include "chainfaker.h"

static const double PI = 3.1415926535;

void LogCoords(const char *Name, coords_t c);
coords_t subtract(const coords_t &a, const coords_t &b);
coords_t add(const coords_t &a, const coords_t &b);
double get_angle(const coords_t &a, const coords_t &b);
coords_t cross_product(const coords_t &a, const coords_t &b);
coords_t normalize(const coords_t &a);
double get_dist(const coords_t &a, const coords_t &b);
coords_t rotate_around_vector(const coords_t &v,
							  const coords_t &axis, double theta_rad);

void GetSCOPFoldFromLabel(const string &Label, string &Fold)
	{
	vector<string> Fields;
	Split(Label, Fields, '/');
	asserta(SIZE(Fields) >= 2);
	const string &Fam = Fields[1];

	vector<string> Fields2;
	Split(Fam, Fields2, '.');
	Fold = Fields2[0] + "." + Fields2[1];
	}

vector<PDBChain *> ChainFaker::m_SCOP40;
vector<vector<uint> > ChainFaker::m_FoldIdxToChainIdxs;
vector<string> ChainFaker::m_Folds;
vector<uint> ChainFaker::m_ChainIdxToFoldIdx;
map<string, uint> ChainFaker::m_FoldToIdx;

bool ChainFaker::SetTakeout()
	{
	for (uint Try = 0; Try < 10; ++Try)
		{
		uint Pos1 = randu32()%m_RL;
		uint Pos2 = randu32()%m_RL;
		uint Lo = min(Pos1, Pos2);
		uint Hi = max(Pos1, Pos2);
		uint FragL = Hi - Lo + 1;
		if (FragL >= m_MinTakeoutLen && FragL <= m_MaxTakeoutLen)
			{
			m_TakeoutLo = Lo;
			m_TakeoutHi = Hi;
			uint TL = m_TakeoutHi - m_TakeoutLo + 1;
			if (m_Trace)
				Log("	Takeout %u .. %u (%u)\n",
					m_TakeoutLo, m_TakeoutHi, TL);
			return true;
			}
		}
	m_TakeoutLo = UINT_MAX;
	m_TakeoutHi = UINT_MAX;
	return false;
	}

void ChainFaker::ReadSCOP40(const string &FN)
	{
	ReadChains(FN, m_SCOP40);
	const uint ChainCount = SIZE(m_SCOP40);

	map<string, vector<PDBChain *> > FoldToChains;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		PDBChain *Chain = m_SCOP40[ChainIdx];
		const string &Label = Chain->m_Label;
		string Fold;
		GetSCOPFoldFromLabel(Label, Fold);
		uint FoldIdx = UINT_MAX;
		if (m_FoldToIdx.find(Fold) == m_FoldToIdx.end())
			{
			FoldIdx = SIZE(m_Folds);
			m_FoldToIdx[Fold] = FoldIdx;
			m_Folds.push_back(Fold);

			vector<uint> Empty;
			m_FoldIdxToChainIdxs.push_back(Empty);
			}
		else
			FoldIdx = m_FoldToIdx[Fold];
		m_FoldIdxToChainIdxs[FoldIdx].push_back(ChainIdx);
		m_ChainIdxToFoldIdx.push_back(FoldIdx);
		}
	}

uint ChainFaker::GetFoldIdx(uint ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_ChainIdxToFoldIdx));
	uint FoldIdx = m_ChainIdxToFoldIdx[ChainIdx];
	return FoldIdx;
	}

const vector<uint> &ChainFaker::GetChainIdxsByFoldIdx(uint FoldIdx) const
	{
	asserta(FoldIdx < SIZE(m_FoldIdxToChainIdxs));
	const vector<uint> &ChainIdxs = m_FoldIdxToChainIdxs[FoldIdx];
	return ChainIdxs;
	}

uint ChainFaker::GetFoldSize(uint FoldIdx) const
	{
	const vector<uint> &ChainIdxs = GetChainIdxsByFoldIdx(FoldIdx);
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

const char *ChainFaker::GetFoldName(uint FoldIdx) const
	{
	asserta(FoldIdx < SIZE(m_Folds));
	return m_Folds[FoldIdx].c_str();
	}

const PDBChain &ChainFaker::GetChain(uint ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_SCOP40));
	return *m_SCOP40[ChainIdx];
	}

void ChainFaker::FindCandidatePlugs1(uint ChainIdx,
									 uint PL_lo, uint PL_hi)
	{
	if (ChainIdx == m_RealChainIdx)
		return;

	const PDBChain &Chain = GetRealChain(ChainIdx);
	uint L = Chain.GetSeqLength();
	for (uint Lo = 0; Lo + PL_lo <= L; ++Lo)
		{
		for (uint Hi = Lo + PL_lo-1; Hi <= Lo + PL_hi-1; ++Hi)
			{
			if (Hi >= L)
				break;
			double d = Chain.GetDist(Lo, Hi);
			double Error = fabs(d - m_TakeoutTermDist);
			if (Error <= m_MaxFitError)
				{
				m_PlugChainIdxs.push_back(ChainIdx);
				m_PlugLos.push_back(Lo);
				m_PlugHis.push_back(Hi);

				if (false && m_Trace)
					{
					Log(" candidate plug");
					Log("  %4u .. %4u", Lo, Hi);
					Log("  d=%.1f, err=%.1f", d, Error);
					Log("  >%s", GetChainLabel(ChainIdx));
					Log("\n");
					}
				}
			}
		}
	}

void ChainFaker::FindCandidatePlugs()
	{
	asserta(m_TakeoutLo < m_TakeoutHi);
	uint TL = m_TakeoutHi - m_TakeoutLo + 1;
	uint PL_min = TL - 2;
	uint PL_max = TL + 4;

	const uint N = SIZE(*m_RealFoldChainIdxs);
	for (uint i = 0; i < N; ++i)
		{
		uint RealChainIdx = (*m_RealFoldChainIdxs)[i];
		FindCandidatePlugs1(RealChainIdx, PL_min, PL_max);
		}
	}

void ChainFaker::MakePlug(uint PlugIdx, PDBChain &Plug) const
	{
	Plug.Clear();
	Ps(Plug.m_Label, "Plug%u", PlugIdx);

	uint PlugChainIdx = m_PlugChainIdxs[PlugIdx];
	uint PlugLo = m_PlugLos[PlugIdx];
	uint PlugHi = m_PlugHis[PlugIdx];

	const PDBChain &Chain = GetChain(PlugChainIdx);
	for (uint i = PlugLo; i <= PlugHi; ++i)
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

double ChainFaker::TryFitPlug1(uint PlugIdx,
	coords_t &t,
	coords_t &axis1,
	double &theta1_rad,
	PDBChain &Plug) const
	{
	MakePlug(PlugIdx, Plug);
	uint PL = Plug.GetSeqLength();
	asserta(SIZE(Plug.m_Xs) == PL);
	asserta(SIZE(Plug.m_Ys) == PL);
	asserta(SIZE(Plug.m_Zs) == PL);

	asserta(m_FakeChain != 0);

	coords_t plo, phi;
	plo.x = Plug.m_Xs[0];
	plo.y = Plug.m_Ys[0];
	plo.z = Plug.m_Zs[0];

	phi.x = Plug.m_Xs[PL-1];
	phi.y = Plug.m_Ys[PL-1];
	phi.z = Plug.m_Zs[PL-1];

	coords_t tlo, thi;
	tlo.x = m_FakeChain->m_Xs[m_TakeoutLo];
	tlo.y = m_FakeChain->m_Ys[m_TakeoutLo];
	tlo.z = m_FakeChain->m_Zs[m_TakeoutLo];

	thi.x = m_FakeChain->m_Xs[m_TakeoutHi];
	thi.y = m_FakeChain->m_Ys[m_TakeoutHi];
	thi.z = m_FakeChain->m_Zs[m_TakeoutHi];

	coords_t centroid_p, centroid_t;
	centroid_p.x = (plo.x + phi.x)/2;
	centroid_p.y = (plo.y + phi.y)/2;
	centroid_p.z = (plo.z + phi.z)/2;

	centroid_t.x = (tlo.x + thi.x)/2;
	centroid_t.y = (tlo.y + thi.y)/2;
	centroid_t.z = (tlo.z + thi.z)/2;
	
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

	coords_t translated_plo, translated_phi;
	translated_plo.x = Plug.m_Xs[0];
	translated_plo.y = Plug.m_Ys[0];
	translated_plo.z = Plug.m_Zs[0];

	translated_phi.x = Plug.m_Xs[PL-1];
	translated_phi.y = Plug.m_Ys[PL-1];
	translated_phi.z = Plug.m_Zs[PL-1];

#if DEBUG
	{
	coords_t translated_centroid_p;
	translated_centroid_p.x += (translated_plo.x + translated_phi.x)/2;
	translated_centroid_p.y += (translated_plo.y + translated_phi.y)/2;
	translated_centroid_p.z += (translated_plo.z + translated_phi.z)/2;

// Validate that centroids coincide after translation
	asserta(feq(translated_centroid_p.x, centroid_t.x));
	asserta(feq(translated_centroid_p.y, centroid_t.y));
	asserta(feq(translated_centroid_p.z, centroid_t.z));
	}
#endif

// Find angle between vectors pv=PlugLo->Hi and tv=TakeoutLo->Hi
	coords_t pv = subtract(translated_phi, translated_plo);
	coords_t tv = subtract(thi, tlo);
	theta1_rad = get_angle(tv, pv);

// Take cross product to get direction perpendicular to plane
// defined by pv and tv
	axis1 = normalize(cross_product(tv, pv));

// Rotate by angle theta around this direction, this brings
// PlugLo close to TakeoutLo and
// PlugHi close to TakeoutHi
	for (uint i = 0; i < PL; ++i)
		{
		coords_t pt = Plug.GetCoords(i);
		coords_t pt2 = subtract(pt, centroid_t);
		coords_t pt3 = rotate_around_vector(pt2, axis1, theta1_rad);
		coords_t pt4 = add(pt3, centroid_t);

		Plug.m_Xs[i] = pt4.x;
		Plug.m_Ys[i] = pt4.y;
		Plug.m_Zs[i] = pt4.z;
		}

	coords_t final_plo, final_phi;
	final_plo.x = Plug.m_Xs[0];
	final_plo.y = Plug.m_Ys[0];
	final_plo.z = Plug.m_Zs[0];

	final_phi.x = Plug.m_Xs[PL-1];
	final_phi.y = Plug.m_Ys[PL-1];
	final_phi.z = Plug.m_Zs[PL-1];

#if DEBUG
	{
	coords_t final_centroid_p;
	final_centroid_p.x += (final_plo.x + final_phi.x)/2;
	final_centroid_p.y += (final_plo.y + final_phi.y)/2;
	final_centroid_p.z += (final_plo.z + final_phi.z)/2;
	double d = get_dist(final_centroid_p, centroid_t);
	asserta(d < 1);
	}
#endif

// These distances should be small
	double dlo = get_dist(tlo, final_plo);
	double dhi = get_dist(thi, final_phi);
	double dtot = dlo + dhi;
	if (dtot > 1)
		{
		Warning("dtot=%.3g", dtot);
		return DBL_MAX;
		}
	if (m_Trace)
		Log("TryFitPlug1() dtot=%.3g\n", dtot);
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

	coords_t plo, phi;
	plo.x = Plug.m_Xs[0];
	plo.y = Plug.m_Ys[0];
	plo.z = Plug.m_Zs[0];

	phi.x = Plug.m_Xs[PL-1];
	phi.y = Plug.m_Ys[PL-1];
	phi.z = Plug.m_Zs[PL-1];

// Centroids of takeout and plug should be identical
// (withing rounding error)
	coords_t centroid;
	centroid.x = (plo.x + phi.x)/2;
	centroid.y = (plo.y + phi.y)/2;
	centroid.z = (plo.z + phi.z)/2;

#if DEBUG
// Verify centroids are the same
	coords_t centroid_t;
	{
	coords_t tlo, thi;
	tlo.x = m_FakeChain->m_Xs[m_TakeoutLo];
	tlo.y = m_FakeChain->m_Ys[m_TakeoutLo];
	tlo.z = m_FakeChain->m_Zs[m_TakeoutLo];

	thi.x = m_FakeChain->m_Xs[m_TakeoutHi];
	thi.y = m_FakeChain->m_Ys[m_TakeoutHi];
	thi.z = m_FakeChain->m_Zs[m_TakeoutHi];

	centroid_t.x = (tlo.x + thi.x)/2;
	centroid_t.y = (tlo.y + thi.y)/2;
	centroid_t.z = (tlo.z + thi.z)/2;
	asserta(fabs(centroid_t.x - centroid.x) < 0.5);
	asserta(fabs(centroid_t.y - centroid.y) < 0.5);
	asserta(fabs(centroid_t.z - centroid.z) < 0.5);
	}
#endif

// Axis of rotation is start-end of takeout (same as 
// start-end of plug)
	coords_t axis = normalize(subtract(phi, plo));

// Rotate by angle theta2 around this direction, this leaves
// coordinates of start, end and centroid unchanged
	for (uint i = 0; i < PL; ++i)
		{
		coords_t pt = Plug.GetCoords(i);
		coords_t pt2 = subtract(pt, centroid_t);
		coords_t pt3 = rotate_around_vector(pt2, axis, theta2_rad);
		coords_t pt4 = add(pt3, centroid_t);

		RotatedPlug.m_Xs.push_back(pt4.x);
		RotatedPlug.m_Ys.push_back(pt4.y);
		RotatedPlug.m_Zs.push_back(pt4.z);
		}

#if DEBUG
// Verify unchanged
	coords_t centroid_r; // _r for rotated
	{
	coords_t rlo, rhi;
	rlo.x = RotatedPlug.m_Xs[0];
	rlo.y = RotatedPlug.m_Ys[0];
	rlo.z = RotatedPlug.m_Zs[0];

	rhi.x = RotatedPlug.m_Xs[PL-1];
	rhi.y = RotatedPlug.m_Ys[PL-1];
	rhi.z = RotatedPlug.m_Zs[PL-1];

	asserta(feq(rlo.x, plo.x));
	asserta(feq(rlo.y, plo.y));
	asserta(feq(rlo.z, plo.z));

	asserta(feq(rhi.x, phi.x));
	asserta(feq(rhi.y, phi.y));
	asserta(feq(rhi.z, phi.z));

	centroid_r.x = (rlo.x + rhi.x)/2;
	centroid_r.y = (rlo.y + rhi.y)/2;
	centroid_r.z = (rlo.z + rhi.z)/2;
	asserta(fabs(centroid_r.x - centroid.x) < 0.5);
	asserta(fabs(centroid_r.y - centroid.y) < 0.5);
	asserta(fabs(centroid_r.z - centroid.z) < 0.5);
	}
#endif
	return;
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

	PDBChain NewChain;
	for (uint i = 0; i < m_TakeoutLo; ++i)
		{
		NewChain.m_Seq += m_FakeChain->m_Seq[i];
		NewChain.m_Xs.push_back(m_FakeChain->m_Xs[i]);
		NewChain.m_Ys.push_back(m_FakeChain->m_Ys[i]);
		NewChain.m_Zs.push_back(m_FakeChain->m_Zs[i]);
		}

	for (uint i = 0; i < PL; ++i)
		{
		NewChain.m_Seq += Plug.m_Seq[i];
		NewChain.m_Xs.push_back(Plug.m_Xs[i]);
		NewChain.m_Ys.push_back(Plug.m_Ys[i]);
		NewChain.m_Zs.push_back(Plug.m_Zs[i]);
		}

	for (uint i = m_TakeoutHi+1; i < FL; ++i)
		{
		NewChain.m_Seq += m_FakeChain->m_Seq[i];
		NewChain.m_Xs.push_back(m_FakeChain->m_Xs[i]);
		NewChain.m_Ys.push_back(m_FakeChain->m_Ys[i]);
		NewChain.m_Zs.push_back(m_FakeChain->m_Zs[i]);
		}
	*m_FakeChain = NewChain;
	}

bool ChainFaker::TryFitPlug(uint PlugIdx,
	coords_t &t,
	coords_t &axis1,
	double &theta1_rad,
	double &theta2_rad,
	PDBChain &Plug) const
	{
	double dtot = TryFitPlug1(PlugIdx, t, axis1, theta1_rad, Plug);
	if (dtot > 1)
		return false;
	bool Ok = TryFitPlug2(Plug, theta2_rad);
	return Ok;
	}

bool ChainFaker::MakeFake(uint ChainIdx, PDBChain &FakeChain)
	{
	Reset();
	m_FakeChain = &FakeChain;
	m_RealChainIdx = ChainIdx;
	m_RealChain = &GetChain(m_RealChainIdx);
	m_RL = m_RealChain->GetSeqLength();
	m_RealFoldIdx = GetFoldIdx(ChainIdx);
	uint FoldSize = GetFoldSize(m_RealFoldIdx);
	if (m_Trace)
		{
		Log("\n");
		Log("MakeFake(ChainIdx=%u)", m_RealChainIdx);
		Log(" fold=%s", GetFoldName(m_RealFoldIdx));
		Log(" size=%u", FoldSize);
		Log("\n");
		}

	if (FoldSize < 3)
		{
		if (m_Trace)
			Log("	**FAIL** size too small\n");
		return false;
		}

	m_RealFoldChainIdxs =
		&GetChainIdxsByFoldIdx(m_RealFoldIdx);
	const PDBChain &RealChain = GetRealChain(m_RealChainIdx);
	*m_FakeChain = RealChain;

	bool Ok = SetTakeout();
	if (!Ok)
		{
		if (m_Trace)
			Log("	**FAIL** SetTakeout()\n");
		return false;
		}
	m_TakeoutTermDist =
		m_FakeChain->GetDist(m_TakeoutLo, m_TakeoutHi);
	if (m_Trace)
		Log("	m_TakeoutTermDist = %.3g\n", m_TakeoutTermDist);

	FindCandidatePlugs();
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

	for (uint i = 0; i < NrCandidatePlugs; ++i)
		{
		coords_t t, axis1;
		double theta1_rad, theta2_rad;
		PDBChain Plug;
		uint PlugIdx = PlugIdxs[i];
		bool Ok = TryFitPlug(PlugIdx, t, axis1, theta1_rad, theta2_rad, Plug);
		if (Ok)
			{
			m_FakeChain->ToPDB("before.pdb");
			ReplaceTakeoutWithPlug(Plug);
			m_FakeChain->ToPDB("after.pdb");
			Die("TODO");
			}
		}

	return true;
	}

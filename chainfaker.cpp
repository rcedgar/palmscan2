#include "myutils.h"
#include "chainfaker.h"

void LogCoords(const char *Name, coords_t c);
coords_t subtract(const coords_t &a, const coords_t &b);
double get_angle(const coords_t &a, const coords_t &b);
coords_t cross_product(const coords_t &a, const coords_t &b);
coords_t normalize(const coords_t &a);
coords_t rotate_around_vector(const coords_t &v,
							  const coords_t &axis, double theta);

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
	uint NrCandidatePlugs = SIZE(m_PlugChainIdxs);
	asserta(SIZE(m_PlugLos) == NrCandidatePlugs);
	asserta(SIZE(m_PlugHis) == NrCandidatePlugs);
	if (m_Trace)
		Log("	%u candidate plugs\n", NrCandidatePlugs);
	if (NrCandidatePlugs == 0)
		return false;

	coords_t t, axis1, axis2;
	double theta1_rad, theta2_rad;
	TryFitPlug(0, t, axis1, theta1_rad, axis2, theta2_rad);
	return true;
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

double ChainFaker::TryFitPlug(uint PlugIdx,
	coords_t &t,
	coords_t &axis1,
	double &theta1_rad,
	coords_t &axis2,
	double &theta2_rad) const
	{
	PDBChain Plug;
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

// Translate Plug so that its centroid coincides
// with Takeout centroid.
	for (uint i = 0; i < PL; ++i)
		{
		Plug.m_Xs[i] += centroid_t.x - centroid_p.x;
		Plug.m_Ys[i] += centroid_t.y - centroid_p.y;
		Plug.m_Zs[i] += centroid_t.z - centroid_p.z;
		}

	coords_t translated_plo, translated_phi;
	translated_plo.x = Plug.m_Xs[0];
	translated_plo.y = Plug.m_Ys[0];
	translated_plo.z = Plug.m_Zs[0];

	translated_phi.x = Plug.m_Xs[PL-1];
	translated_phi.y = Plug.m_Ys[PL-1];
	translated_phi.z = Plug.m_Zs[PL-1];

// Validate that centroids coincide after translation
	coords_t translated_centroid_p;
	translated_centroid_p.x += (translated_plo.x + translated_phi.x)/2;
	translated_centroid_p.y += (translated_plo.y + translated_phi.y)/2;
	translated_centroid_p.z += (translated_plo.z + translated_phi.z)/2;

	asserta(feq(translated_centroid_p.x, centroid_t.x));
	asserta(feq(translated_centroid_p.y, centroid_t.y));
	asserta(feq(translated_centroid_p.z, centroid_t.z));

// Find angle between vectors pv=PlugLo->Hi and tv=TakeoutLo->Hi
	coords_t pv = subtract(phi, plo);
	coords_t tv = subtract(thi, tlo);
	theta1_rad = get_angle(tv, pv);

// Take cross product to get axis perpendicular to plane
// defined by pv and tv
	coords_t axis = normalize(cross_product(tv, pv));

// Rotate by angle theta around this axis, this brings
// PlugLo,TakeoutLo very close and similarly PlugHi,TakeoutHi.

	for (uint i = 0; i < PL; ++i)
		{
		coords_t pt = Plug.GetCoords(i);
		coords_t pt_rot = rotate_around_vector(pt, axis, -theta1_rad);

		Plug.m_Xs.push_back(pt.x);
		Plug.m_Ys.push_back(pt.y);
		Plug.m_Zs.push_back(pt.z);
		}

	coords_t final_plo, final_phi;
	final_plo.x = Plug.m_Xs[0];
	final_plo.y = Plug.m_Ys[0];
	final_plo.z = Plug.m_Zs[0];

	final_phi.x = Plug.m_Xs[PL-1];
	final_phi.y = Plug.m_Ys[PL-1];
	final_phi.z = Plug.m_Zs[PL-1];

	coords_t final_centroid_p;
	final_centroid_p.x += (final_plo.x + final_phi.x)/2;
	final_centroid_p.y += (final_plo.y + final_phi.y)/2;
	final_centroid_p.z += (final_plo.z + final_phi.z)/2;

	LogCoords("tlo", tlo);
	LogCoords("thi", thi);
	LogCoords("centroid_t", centroid_t);
	LogCoords("translated_plo", translated_plo);
	LogCoords("translated_phi", translated_phi);
	LogCoords("translated_centroid", translated_centroid_p);
	LogCoords("final_plo", final_plo);
	LogCoords("final_phi", final_phi);
	LogCoords("final_centroid", final_centroid_p);

	return 0;
	}

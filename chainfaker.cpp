#include "myutils.h"
#include "chainfaker.h"

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

double ChainFaker::TryFitPlug(uint PlugIdx,
							  double t[3], double u[3][3]) const
	{
	return 0;
	}

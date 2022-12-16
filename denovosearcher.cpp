#include "myutils.h"
#include "denovosearcher.h"

static const uint min_aadist_AdBg = 10;
static const uint max_aadist_AdBg = 80;

static const uint min_aadist_BgCd = 10;
static const uint max_aadist_BgCd = 80;

static const uint min_aadist_AdCd = 80;
static const uint max_aadist_AdCd = 200;

void DeNovoSearcher::SearchAd(uint AdLo, uint AdHi,
  vector<uint> &PosVec)
	{
	PosVec.clear();
	const uint L = SIZE(m_Seq);
	for (uint Pos = 3; Pos + 13 < L; ++Pos)
		{
		if (m_Seq[Pos] != 'D')
			continue;
		bool Ok = MatchAd_SS(Pos);
		if (Ok)
			PosVec.push_back(Pos);
		}
	}

void DeNovoSearcher::SearchBg(uint BgLo, uint BgHi,
  vector<uint> &PosVec)
	{
	PosVec.clear();
	const uint L = SIZE(m_Seq);
	for (uint Pos = 3; Pos + 13 < L; ++Pos)
		{
		if (m_Seq[Pos] != 'G')
			continue;
		bool Ok = MatchBg_SS(Pos);
		if (Ok)
			PosVec.push_back(Pos);
		}
	}

void DeNovoSearcher::SearchCd(uint CdLo, uint CdHi,
  vector<uint> &PosVec)
	{
	PosVec.clear();
	const uint L = SIZE(m_Seq);
	for (uint Pos = 3; Pos + 13 < L; ++Pos)
		{
		if (m_Seq[Pos] != 'D')
			continue;
		bool Ok = MatchCd_SS(Pos);
		if (Ok)
			PosVec.push_back(Pos);
		}
	}

void DeNovoSearcher::GetBgLoHi(uint Ad, uint &BgLo, uint &BgHi) const
	{
	BgLo = Ad + min_aadist_AdBg;
	BgHi = Ad + max_aadist_AdBg;
	}

void DeNovoSearcher::GetCdLoHi(uint Bg, uint &CdLo, uint &CdHi) const
	{
	CdLo = Bg + min_aadist_BgCd;
	CdHi = Bg + max_aadist_BgCd;
	}

void DeNovoSearcher::Search(PDBChain &Query)
	{
	const uint QL = SIZE(Query.m_Seq);
	asserta(SIZE(Query.m_SS) == QL);

	Clear();
	m_Query = &Query;
	m_Seq = Query.m_Seq;
	m_SS = Query.m_SS;

	vector<uint> AdVec;
	vector<uint> BgVec;
	vector<uint> CdVec;

	SearchAd(0, QL, AdVec);
	const uint nad = SIZE(AdVec);
	for (uint iad = 0; iad < nad; ++iad)
		{
		uint Ad = AdVec[iad];
		uint BgLo, BgHi;
		GetBgLoHi(Ad, BgLo, BgHi);

		SearchBg(BgLo, BgHi, BgVec);
		const uint nbg = SIZE(BgVec);
		for (uint ibg = 0; ibg < nbg; ++ibg)
			{
			uint Bg = BgVec[iad];
			uint CdLo, CdHi;
			GetCdLoHi(Bg, CdLo, CdHi);

			SearchCd(CdLo, CdHi, CdVec);
			const uint ncd = SIZE(CdVec);
			for (uint icd = 0; icd < ncd; ++icd)
				{
				uint Cd = CdVec[icd];
				bool Ok = CheckHit(Ad, Bg, Cd);
				if (Ok)
					{
					m_Ads.push_back(Ad);
					m_Bgs.push_back(Bg);
					m_Cds.push_back(Cd);
					}
				}
			}
		}
	}

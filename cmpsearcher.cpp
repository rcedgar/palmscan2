#include "myutils.h"
#include "cmpsearcher.h"

/***
From D:\src\py\palmscan3d_params.py
-----------------------------------
static const uint min_aadist_AdBg = 21;
static const uint max_aadist_AdBg = 136;

static const uint min_aadist_BgCd = 28;
static const uint max_aadist_BgCd = 140;

static const uint min_aadist_AdCd = 85;
static const uint max_aadist_AdCd = 173;
***/

// Manually extended to allow outliers
static const uint min_aadist_AdBg = 10;
static const uint max_aadist_AdBg = 150;

static const uint min_aadist_BgCd = 10;
static const uint max_aadist_BgCd = 150;

//static const uint min_aadist_AdCd = 70;
//static const uint max_aadist_AdCd = 200;

static const double MINSCORE1 = 0.5;
static const double MINSCORE3 = 0.5;

static const double MAXSCORE1 = 2;
static const double MAXSCORE3 = 2.5;

#define TRACE	0

bool CMPSearcher::GoodScore1(double Score) const
	{
	if (NoStdDevs())
		{
		if (Score <= MAXSCORE1)
			return true;
		}
	else
		{
		if (Score >= MINSCORE1)
			return true;
		}
	return false;
	}

bool CMPSearcher::GoodScore3(double Score) const
	{
	if (NoStdDevs())
		{
		if (Score <= MAXSCORE3)
			return true;
		}
	else
		{
		if (Score >= MINSCORE3)
			return true;
		}
	return false;
	}

bool CMPSearcher::MatchAd(uint Pos) const
	{
	if (Pos < 3)
		return false;
	double Score = m_Prof->GetScoreA(*m_Query, Pos-3);
	if (GoodScore1(Score))
		return true;
	return false;
	}

bool CMPSearcher::MatchBg(uint Pos) const
	{
	if (Pos < 1)
		return false;
	double Score = m_Prof->GetScoreB(*m_Query, Pos-1);
	if (GoodScore1(Score))
		return true;
	return false;
	}

bool CMPSearcher::MatchCd(uint Pos) const
	{
	if (Pos < 3)
		return false;
	double Score = m_Prof->GetScoreC(*m_Query, Pos-3);
	if (GoodScore1(Score))
		return true;
	return false;
	}

void CMPSearcher::SearchAd(uint AdLo, uint AdHi,
  vector<uint> &PosVec)
	{
	PosVec.clear();
	const uint L = SIZE(m_Seq);
	for (uint Pos = 3; Pos + AL+BL+CL < L; ++Pos)
		{
		if (m_Seq[Pos] != 'D')
			continue;
		bool Ok = MatchAd(Pos);
		if (Ok)
			PosVec.push_back(Pos);
		}
	}

void CMPSearcher::SearchBg(uint BgLo, uint BgHi,
  vector<uint> &PosVec)
	{
	PosVec.clear();
	const uint L = SIZE(m_Seq);
	for (uint Pos = BgLo; Pos <= BgHi && Pos + BL+CL < L; ++Pos)
		{
		if (m_Seq[Pos] != 'G')
			continue;
		bool Ok = MatchBg(Pos);
		if (Ok)
			PosVec.push_back(Pos);
		}
	}

void CMPSearcher::SearchCd(uint CdLo, uint CdHi,
  vector<uint> &PosVec)
	{
	PosVec.clear();
	const uint L = SIZE(m_Seq);
	for (uint Pos = CdLo; Pos <= CdHi && Pos + 4 < L; ++Pos)
		{
		if (m_Seq[Pos] != 'D')
			continue;
		bool Ok = MatchCd(Pos);
		if (Ok)
			PosVec.push_back(Pos);
		}
	}

void CMPSearcher::GetBgLoHi(uint Ad, uint &BgLo, uint &BgHi) const
	{
	BgLo = Ad + min_aadist_AdBg;
	BgHi = Ad + max_aadist_AdBg;
	}

void CMPSearcher::GetCdLoHi(uint Bg, uint &CdLo, uint &CdHi) const
	{
	CdLo = Bg + min_aadist_BgCd;
	CdHi = Bg + max_aadist_BgCd;
	}

void CMPSearcher::Search(PDBChain &Query)
	{
	ClearSearch();
	m_Query = &Query;
	m_Seq = Query.m_Seq;

	if (!opt_permuted)
		Search_ABC(Query);

#if 0 // TODO
	if (!opt_notpermuted)
		Search_CAB(Query);
#endif
	}

void CMPSearcher::Search_ABC(PDBChain &Query)
	{
	const uint QL = SIZE(Query.m_Seq);

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
			uint Bg = BgVec[ibg];
			uint CdLo, CdHi;
			GetCdLoHi(Bg, CdLo, CdHi);

			SearchCd(CdLo, CdHi, CdVec);
			const uint ncd = SIZE(CdVec);
			for (uint icd = 0; icd < ncd; ++icd)
				{
				uint Cd = CdVec[icd];
				double Score = GetScore(Ad, Bg, Cd);
#if TRACE
				Log("[Score3] %.4f\n", Score);
#endif
				CheckHit(Ad, Bg, Cd, Score);
				}
			}
		}
	}

double CMPSearcher::GetPSSMStarts(uint &PosA, uint &PosB, uint &PosC) const
	{
	PosA = UINT_MAX;
	PosB = UINT_MAX;
	PosC = UINT_MAX;
	const uint N = SIZE(m_Scores);
	if (N == 0)
		return -1;
	double TopScore = m_Scores[0];
	PosA = m_Ads[0];
	PosB = m_Bgs[0];
	PosC = m_Cds[0];

	const bool NoS = NoStdDevs();

	for (uint i = 1; i < N; ++i)
		{
		bool Better = (NoS ? m_Scores[i] < TopScore : m_Scores[i] > TopScore);
		if (Better)
			{
			TopScore = m_Scores[i];
			PosA = m_Ads[i];
			PosB = m_Bgs[i];
			PosC = m_Cds[i];
			}
		}

	asserta(PosA >= 3);
	asserta(PosB >= 1);
	asserta(PosC >= 3);
	PosA -= 3;
	PosB -= 1;
	PosC -= 3;

	const string &Q = m_Query->m_Seq;
	asserta(Q[PosA+3] == 'D');
	asserta(Q[PosB+1] == 'G');
	asserta(Q[PosC+3] == 'D');

	return TopScore;
	}

double CMPSearcher::GetScore(uint Ad, uint Bg, uint Cd) const
	{
	const uint QL = SIZE(m_Query->m_Seq);
	asserta(Ad >= 3 && Ad < QL);
	asserta(Bg >= 2 && Bg < QL);
	asserta(Cd >= 3 && Cd < QL);
	double Score = m_Prof->GetScore3(*m_Query, Ad-3, Bg-1, Cd-3);
	return Score;
	}

void CMPSearcher::CheckHit(uint Ad, uint Bg, uint Cd, double Score)
	{
	if (!GoodScore3(Score))
		return;

	m_Ads.push_back(Ad);
	m_Bgs.push_back(Bg);
	m_Cds.push_back(Cd);
	m_Scores.push_back(Score);
	}

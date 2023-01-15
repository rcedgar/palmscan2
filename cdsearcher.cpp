#include "myutils.h"
#include "cdsearcher.h"

#define TRACE	0

double GetNormal(double Mu, double Sigma, double x);

void CDSearcher::LogHits() const
	{
	Log("\n");
	Log(">%s\n", m_Query->m_Label.c_str());
	const uint MotifCount = m_Info->GetMotifCount();
	asserta(SIZE(m_MotifIndexToHits) == MotifCount);
	asserta(SIZE(m_MotifIndexToScores) == MotifCount);
	for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
		{
		const vector<uint> &Hits = m_MotifIndexToHits[MotifIndex];
		const vector<double> &Scores = m_MotifIndexToScores[MotifIndex];
		uint n = SIZE(Hits);
		asserta(SIZE(Scores) == n);
		Log("  %2s", m_Info->GetMotifName(MotifIndex));
		Log("  [%3u]", n);
		for (uint i = 0; i < n; ++i)
			Log("  %u=%.4f", Hits[i], Scores[i]);
		Log("\n");
		}
	}

void CDSearcher::ClearSearch()
	{
	m_Query = 0;
	m_MotifStarts.clear();
	m_MotifScores.clear();
	m_MotifIndexToHits.clear();
	m_MotifIndexToScores.clear();

	const uint MotifCount = m_Info->GetMotifCount();
	m_MotifStarts.resize(MotifCount, UINT_MAX);
	m_MotifScores.resize(MotifCount, 0);
	m_MotifIndexToHits.resize(MotifCount);
	m_MotifIndexToScores.resize(MotifCount);
	}

void CDSearcher::Init(const CDInfo &Info, const CDData &Dists,
  const CDData &StdDevs)
	{
	Clear();
	m_Info = &Info;
	m_Dists = &Dists;
	m_StdDevs = &StdDevs;
	ClearSearch();
	}

double CDSearcher::GetScore(
  uint SeqPos1, uint SeqPos2,
  uint Ix1, uint Ix2,
  uint L1, uint L2) const
	{
	bool Diag = (Ix1 == Ix2);
	double XS = (Diag ? 1.5 : 2);
	double Sum = 0;
	uint n = 0;
	asserta(m_Dists != 0);
	asserta(m_StdDevs != 0);
	const CDData &Dists = *m_Dists;
	const CDData &StdDevs = *m_StdDevs;
	for (uint i = 0; i < L1; ++i)
		{
		uint jhi = (Diag ? i : L2);
		for (uint j = 0; j < jhi; ++j)
			{
			double Observed_d = m_Query->GetDist(SeqPos1+i, SeqPos2+j);
			double Mu = Dists.GetByIx(Ix1 + i, Ix2 + j);
			double Sigma = StdDevs.GetByIx(Ix1 + i, Ix2 + j);
			double y = GetNormal(Mu, XS*Sigma, Observed_d);
			double Max = GetNormal(Mu, XS*Sigma, Mu);
			double Ratio = y/Max;
			Sum += Ratio;
			++n;
			}
		}
	asserta(n > 0);
	double Score = Sum/n;
	return Score;
	}

/***
Minimum start position of next motif is:
	{Position of lowest hit to previous motif} +
	  {minimum number of AAs to next motif}

Maximum start position of next motif is:
	{Position of highest hit to previous motif} +
	  {maximum number of AAs to next motif}
***/
void CDSearcher::GetRangeNext(uint NextMotifIndex, uint &Lo, uint &Hi) const
	{
	Lo = UINT_MAX;
	Hi = UINT_MAX;

	uint QL = m_Query->GetSeqLength();
	if (NextMotifIndex == 0)
		{
		uint n = m_Template->GetMinAAsSeqEnd(0);
		if (n > QL)
			return;
		Lo = 0;
		Hi = QL - n;
#if TRACE
		Log("  GetRange(0) 0 .. QL-n=%u\n", Hi);
#endif
		return;
		}

	uint PrevMotifIndex = NextMotifIndex-1;
	const vector<uint> &Hits = m_MotifIndexToHits[PrevMotifIndex];
#if TRACE
	Log("  GetRange(%u), %u prev hits\n", NextMotifIndex, SIZE(Hits));
#endif
	if (Hits.empty())
		return;

	uint FirstHit = Hits.front();
	uint LastHit = Hits.back();
#if TRACE
	Log("  First %u, last %u\n", FirstHit, LastHit);
#endif

	uint MinAAsNext = 
	  m_Template->GetMinAAsNext(PrevMotifIndex);
#if TRACE
	Log("  MinAAsNext %u\n", MinAAsNext);
#endif
	if (MinAAsNext == UINT_MAX)
		return;

	uint ML = m_Info->GetMotifLength(NextMotifIndex);
	Lo = FirstHit + MinAAsNext;
	if (Lo + ML >= QL)
		{
		Lo = UINT_MAX;
		return;
		}

	Hi = LastHit + m_Template->GetMaxAAsNext(PrevMotifIndex);
	if (Hi + ML >= QL)
		Hi = QL - ML;
#if TRACE
	Log("  Hi %u\n", Hi);
#endif
	}

void CDSearcher::Search1(uint MotifIndex, uint Lo, uint Hi,
  vector<uint> &Hits, vector<double> &Scores)
	{
	const char AnchorAA = m_Template->m_AnchorAAs[MotifIndex];
	const uint AnchorAAOffset = m_Template->m_AnchorAAOffsets[MotifIndex];
	double MinScore = m_Template->m_MinScores[MotifIndex];
#if TRACE
	Log("Search1(Mf=%u, Lo=%u, Hi=%u) aa=%c(%u)\n",
	  MotifIndex, Lo, Hi, AnchorAA, AnchorAAOffset);
#endif
	const string &Seq = m_Query->m_Seq;
	const uint QL = SIZE(Seq);
	asserta(Lo <= Hi && Hi < QL);
	const uint Ix = m_Info->GetIx(MotifIndex, 0);
	const uint ML = m_Info->GetMotifLength(MotifIndex);
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		if (AnchorAA != 'x' && Seq[Pos+AnchorAAOffset] != AnchorAA)
			continue;

		double Score = GetScore(Pos, Pos, Ix, Ix, ML, ML);
		if (Score >= MinScore)
			{
#if TRACE
			Log("Mf=%u Pos=%u, score=%.4f hits=%u\n",
			  MotifIndex, Pos, Score, SIZE(Hits));
#endif
			if (!Hits.empty())
				{
				uint LastHit = Hits.back();
				double LastScore = Scores.back();
				if (Pos - LastHit < 10)
					{
					if (Score > LastScore)
						{
						Hits.back() = Pos;
						Scores.back() = Score;
						}
					continue;
					}
				}
			Hits.push_back(Pos);
			Scores.push_back(Score);
			}
		}
	}

void CDSearcher::Search(const PDBChain &Q)
	{
#if TRACE
	Log("\n");
	Log("Search(%s)\n", Q.m_Label.c_str());
#endif
	ClearSearch();

	m_Query = &Q;
	const uint QL = Q.GetSeqLength();
	const uint MotifCount = m_Info->GetMotifCount();

	for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
		{
		uint Lo, Hi;
		GetRangeNext(MotifIndex, Lo, Hi);
#if TRACE
		{
		const char *Name = m_Info->GetMotifName(MotifIndex);
		Log("Motif %s, range %u - %u\n", Name, Lo, Hi);
		}
#endif
		if (Lo == UINT_MAX)
			return;

		vector<uint> &Hits = m_MotifIndexToHits[MotifIndex];
		vector<double> &Scores = m_MotifIndexToScores[MotifIndex];
		Search1(MotifIndex, Lo, Hi, Hits, Scores);
#if TRACE
		{
		Log("  %u hits ", SIZE(Hits));
		for (uint i = 0; i < SIZE(Hits); ++i)
			Log(" %u", Hits[i]);
		Log("\n");
		}
#endif
		if (Hits.empty())
			return;
		}
	}

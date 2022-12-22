#pragma once

#include "pdbchain.h"
#include "ppsp.h"

class DSHit;

class PPSPSearcher
	{
public:
	PDBChain *m_Query = 0;
	string m_Seq;
	const PPSP *m_Prof = 0;

	vector<uint> m_Ads;
	vector<uint> m_Bgs;
	vector<uint> m_Cds;
	vector<double> m_Scores;

public:
	void ClearSearch()
		{
		m_Query = 0;
		m_Seq.clear();
		m_Ads.clear();
		m_Bgs.clear();
		m_Cds.clear();
		m_Scores.clear();
		}

	void Search(PDBChain &Query);
	void Search_ABC(PDBChain &Query);
	void Search_CAB(PDBChain &Query);
	double GetPSSMStarts(uint &PosA, uint &PosB, uint &PosC) const;

	void SearchAd(uint AdLo, uint AdHi, vector<uint> &PosVec);
	void SearchBg(uint BgLo, uint BgHi, vector<uint> &PosVec);
	void SearchCd(uint CdLo, uint CdHi, vector<uint> &PosVec);

	bool MatchAd(uint AdPos) const;
	bool MatchBg(uint BgPos) const;
	bool MatchCd(uint CdPos) const;

	void GetBgLoHi(uint Ad, uint &BgLo, uint &BgHi) const;
	void GetCdLoHi(uint Bg, uint &CdLo, uint &CdHi) const;

	void GetAdLoHi_Permuted(uint Cd, uint &AdLo, uint &AdHi) const;
	void GetBgLoHi_Permuted(uint Ad, uint &BgLo, uint &BgHi) const;

	double GetScore(uint Ad, uint Bg, uint Cd) const;
	void CheckHit(uint Ad, uint Bg, uint Cd, double Score);
	};

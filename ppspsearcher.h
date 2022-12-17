#pragma once
#pragma once

#include "pdbchain.h"
#include "ppsp.h"

class DSHit;

class PPSPSearcher
	{
public:
	PDBChain *m_Query = 0;
	string m_Seq;
	PPSP m_Prof;

	vector<uint> m_Ads;
	vector<uint> m_Bgs;
	vector<uint> m_Cds;
	vector<double> m_Scores;

public:
	void Clear()
		{
		m_Query = 0;
		m_Seq.clear();
		m_Ads.clear();
		m_Bgs.clear();
		m_Cds.clear();
		}

	void Search(PDBChain &Query);

	void SearchAd(uint AdLo, uint AdHi, vector<uint> &PosVec);
	void SearchBg(uint BgLo, uint BgHi, vector<uint> &PosVec);
	void SearchCd(uint CdLo, uint CdHi, vector<uint> &PosVec);

	bool MatchAd(uint AdPos) const;
	bool MatchBg(uint BgPos) const;
	bool MatchCd(uint CdPos) const;

	void GetBgLoHi(uint Ad, uint &BgLo, uint &BgHi) const;
	void GetCdLoHi(uint Bg, uint &CdLo, uint &CdHi) const;

	double GetScore(uint Ad, uint Bg, uint Cd) const;
	void CheckHit(uint Ad, uint Bg, uint Cd, double Score);
	};

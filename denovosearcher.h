#pragma once
#pragma once

#include "pdbchain.h"

class DSHit;

class DeNovoSearcher
	{
public:
	PDBChain *m_Query = 0;
	string m_Seq;
	string m_SS;

	vector<uint> m_Ads;
	vector<uint> m_Bgs;
	vector<uint> m_Cds;

public:
	void Clear()
		{
		m_Query = 0;
		m_Seq.clear();
		m_SS.clear();
		m_Ads.clear();
		m_Bgs.clear();
		m_Cds.clear();
		}

	void LogMe(FILE *f = g_fLog) const;
	void Search(PDBChain &Query);

	void SearchAd(uint AdLo, uint AdHi, vector<uint> &PosVec);
	void SearchBg(uint BgLo, uint BgHi, vector<uint> &PosVec);
	void SearchCd(uint CdLo, uint CdHi, vector<uint> &PosVec);

	bool MatchAd_SS(uint AdPos) const;
	bool MatchBg_SS(uint BgPos) const;
	bool MatchCd_SS(uint CdPos) const;

	void GetBgLoHi(uint Ad, uint &BgLo, uint &BgHi) const;
	void GetCdLoHi(uint Bg, uint &CdLo, uint &CdHi) const;

	bool CheckHit(uint Ad, uint Bg, uint Cd);
	};

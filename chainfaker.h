#pragma once
#pragma once

#include <set>
#include <map>
#include "pt.h"
#include "pdbchain.h"

class ChainFaker
	{
public:
	ChainFaker() {}
	~ChainFaker() {}

public:
	static vector<PDBChain *> m_SCOP40;
	static vector<string> m_Folds;
	static vector<uint> m_ChainIdxToFoldIdx;
	static vector<vector<uint> > m_FoldIdxToChainIdxs;
	static map<string, uint> m_FoldToIdx;
	bool m_Trace = true;

public:
	double m_MinNENDist = 4;
	double m_CADist = 3.81;
	uint m_MinTakeoutLen = 4;
	uint m_MaxTakeoutLen = 16;
	double m_MaxFitError = 0.1;

// Fold data, input to construction
public:
	const PDBChain *m_RealChain = 0;
	uint m_RL = 0;
	uint m_RealChainIdx = UINT_MAX;
	uint m_RealFoldIdx = UINT_MAX;
	const vector<uint> *m_RealFoldChainIdxs = 0;

// Construction state
public:
	PDBChain *m_FakeChain = 0;

// Fragment to be replaced
	uint m_TakeoutLo = UINT_MAX;
	uint m_TakeoutHi = UINT_MAX;
	double m_TakeoutTermDist = DBL_MAX;

// Candidate plugs to replace the fragment
	vector<uint> m_PlugChainIdxs;
	vector<uint> m_PlugLos;
	vector<uint> m_PlugHis;

public:
	void Reset()
		{
		m_RealChain = 0;
		m_FakeChain = 0;
		m_RL = UINT_MAX;
		m_RealChainIdx = UINT_MAX;
		m_RealFoldIdx = UINT_MAX;
		m_RealFoldChainIdxs = 0;
		m_FakeChain = 0;
		m_TakeoutLo = UINT_MAX;
		m_TakeoutHi = UINT_MAX;
		m_TakeoutTermDist = DBL_MAX;
		m_PlugChainIdxs.clear();
		m_PlugLos.clear();
		m_PlugHis.clear();
		}

	void ReadSCOP40(const string &FN);
	const PDBChain &GetChain(uint ChainIdx) const;
	uint GetFoldIdx(uint ChainIdx) const;
	uint GetFoldSize(uint FoldIdx) const;
	const vector<uint> &GetChainIdxsByFoldIdx(uint FoldIdx) const;
	bool MakeFake(uint ChainIdx, PDBChain &Chain);
	const PDBChain &GetRealChain(uint ChainIdx) const;
	const char *GetFoldName(uint FoldIdx) const;
	const char *GetChainLabel(uint ChainIdx) const;
	bool SetTakeout();
	void FindCandidatePlugs();
	void FindCandidatePlugs1(uint ChainIdx, uint PL_min, uint PL_max);
	void MakePlug(uint PlugIdx, PDBChain &Plug) const;
	double TryFitPlug(uint PlugIdx,
					  coords_t &t,
					  coords_t &axis1,
					  double &theta1_rad,
					  coords_t &axis2,
					  double &theta2_rad) const;

	double TryFitPlug1(uint PlugIdx,
					  coords_t &t,
					  coords_t &axis1,
					  double &theta1_rad,
					   PDBChain &Plug) const;
	double TryFitPlug2(PDBChain &Plug,
					  double &theta2_rad) const;
	double FindBadNENDist(const PDBChain &Plug, uint &FakePos, uint &PlugPos) const;
	void RotatePlug2(const PDBChain &Plug, double theta2_rad,
					PDBChain &RotatedPlug) const;

	};

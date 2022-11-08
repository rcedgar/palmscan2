#pragma once

#include "xdpmem.h"
#include "pathinfo.h"
#include "pdb.h"

class TriSearcher
	{
public:
	double MaxTriRMSD2 = 0;
	double Radius = 0;
	uint NABmin = 0;
	uint NABmax = 0;
	uint NBCmin = 0;
	uint NBCmax = 0;
	uint NACmin = 0;
	uint NACmax = 0;

	const PDB *m_Query = 0;
	const PDB *m_Ref = 0;
	vector<uint> m_PosAs;
	vector<uint> m_PosBs;
	vector<uint> m_PosCs;
	vector<double> m_TriRMSD2s;
	vector<double> m_MotifRMSD2s;
	vector<uint> m_HitOrder;

	XDPMem m_XDPMem;

	vector<double> m_TriForm_t;
	vector<vector<double> > m_TriForm_R;

	float **m_DPScoreMx = 0;

public:
	void Clear()
		{
		m_Query = 0;
		m_Ref = 0;
		m_PosAs.clear();
		m_PosBs.clear();
		m_PosCs.clear();
		m_TriRMSD2s.clear();
		m_MotifRMSD2s.clear();
		m_HitOrder.clear();
		}

	void LogMe(FILE *f = g_fLog) const;
	void Search(const PDB &Query, const PDB &Ref);
	void SetHitOrder();
	double GetRMSDMotifs(uint QueryPosA, uint QueryPosB,
	  uint QueryPosC) const;
	double GetRMSD2Segment(uint QueryStartPos, uint RefStartPos, uint n,
	  const vector<double> &t, const vector<vector<double> > &R) const;
	double GetMotifsTM() const;
	double GetTMSum(uint QPos, uint RPos, uint n) const;
	bool GetTopHit(uint &QueryPosA, uint &QueryPosB, uint &QueryPosC,
	  double &MotifRMSD2) const;
	bool AlignPalm(uint QueryPosA, uint QueryPosB, uint QueryPosC,
	  string &Path);
	void AlignSeg(uint QPos, uint QLen, uint RPos, uint RLen,
	  const vector<double> &t, const vector<vector<double> > &R,
	  string &Path);
	void SetTriForm();
	void AllocDPScoreMx();
	};

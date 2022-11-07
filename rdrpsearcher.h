#pragma once

#include "rdrpmodel.h"
#include "palmhit.h"

#define MotLet "ABC"[MotifIndex]

class RdRpSearcher
	{
public:
	const RdRpModel *m_Model = 0;

	string m_QueryLabel;
	string m_QuerySeq;

	PalmHit m_TopPalmHit;
	PalmHit m_SecondPalmHit;

	uint m_MaxX = 10;
	float m_MinCScore = 3.0;
	float m_MinPSSMScore = 2.0;

	uint m_MinPPLength = 90;
	uint m_MaxPPLength = 180;

	uint m_MinInsert = 10;

	bool m_Trace = false;

public:
	RdRpSearcher()
		{
		m_MaxX = 10;
		if (optset_maxx)
			m_MaxX = opt_maxx;
		if (optset_mincscore)
			m_MinCScore = (float) opt_mincscore;
		asserta(m_MinCScore > 0);
		}

	void Clear()
		{
		m_QueryLabel.clear();
		m_QuerySeq.clear();
		m_TopPalmHit.Clear();
		m_SecondPalmHit.Clear();
		}

	void ClearSearch()
		{
		m_QueryLabel.clear();
		m_QuerySeq.clear();
		m_TopPalmHit.Clear();
		m_SecondPalmHit.Clear();
		}

	uint GetGroupCount() const
		{
		asserta(m_Model != 0);
		return m_Model->GetGroupCount();
		}

	void GetGroupName(uint GroupIndex, string &Name) const
		{
		asserta(m_Model != 0);
		return m_Model->GetGroupName(GroupIndex, Name);
		}

public:
	void Init(const RdRpModel &Model);
	void Search(const string &QueryLabel, const string &QuerySeq);
	void SearchGroup(uint GroupIndex);
	void SearchMotif(uint GroupIndex, uint MotifIndex,
	  uint QLo, uint QHi, float MinScore, vector<RPHit> &Hits);
	void SearchMotif_TopHitOnly(uint GroupIndex, uint MotifIndex,
	  uint QLo, uint QHi, float MinScore, RPHit &Hit);
	void SearchAB(uint GroupIndex, const RPHit &Hit_C);
	bool SearchAB_NotPermuted(uint GroupIndex, const RPHit &Hit_C);
	bool SearchAB_Permuted(uint GroupIndex, const RPHit &Hit_C);
	void OnPalmHit(uint GroupIndex, const RPHit &Hit_A, const RPHit &Hit_B,
	  const RPHit &Hit_C, bool Permuted);
	void GetAlnRows(vector<string> &Rows) const;
	void GetMotifsSeq(const string &Sep, string &Seq) const;
	uint GetMotifPos(uint MotifIndex) const;
	void GetTrimmedSeq(string &Seq) const;
	void GetMotifPositions(uint &APos, uint &BPos, uint &CPos) const;
	void GetAln(const RPHit &Hit, string &Q, string &P, string &Annot,
	  uint &QLo, uint &QHi) const;
	uint GetPSSMLength(uint GroupIndex, uint MotifIndex) const;
	void GetSpan(uint &QLo, uint &QHi) const;
	void GetMotifSeq(uint MotifIndex, string &Seq) const;
	void WriteOutput() const;
	void WriteReport(FILE *f) const;
	void WriteTsv(FILE *f) const;

	const PSSM &GetPSSM(uint GroupIndex, uint MotifIndex) const;

public:
	static void InitOutput();
	static void CloseOutput();
	};

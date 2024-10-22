#pragma once

#include "seqdb.h"
#include "rdrpmodel.h"
#include "rdrpsearcher.h"

static const uint LA = 12;
static const uint LB = 14;
static const uint LC = 8;
static const uint NPOS = 3;

class MSAQC2
	{
public:
	RdRpModel m_Mod;
	RdRpSearcher m_RS;
	const SeqDB *m_MSA = 0;
	uint m_ColCount = 0;
	uint m_SeqCount = 0;
	vector<string> m_AlignedSeqs;
	vector<string> m_UnalignedSeqs;
	vector<vector<uint> > m_UnalignedPosToAlignedColVec;
	vector<vector<uint> > m_AlignedColToUnalignedPosVec;

	vector<uint> m_Ls;
	vector<vector<uint> > m_SeqIndexToPosVec;
	vector<vector<uint> > m_SeqIndexToColVec;
	vector<uint> m_ConsensusColVec;

	vector<float> m_ColIndexToGapFract;

public:
	void Init(const SeqDB &MSA);
	uint GetMotifIndex(char c) const;
	uint GetTopIx(const vector<uint> &v, uint MinCount = UINT_MAX) const;
	char mtoc(uint m) const { asserta(m < 3); return "ABC"[m]; }
	bool CheckIncreasingOrder(const vector<uint> &v) const;
	void GetCorrectLetterCount(uint &CorrectCount, uint &LetterCount) const;
	void GetCorrectlyAlignedLetterCount(uint SeqIndex,
	  uint &CorrectCount, uint &LetterCount) const;
	uint GetCorrectSeqCount() const;
	float GetGapFract(uint ColIndex) const;
	void WriteAlnA(FILE *f) const;
	void WriteAlnB(FILE *f) const;
	void WriteAlnC(FILE *f) const;
	bool GetColIsMatchState(uint ColIndex) const;
	uint GetColIsMatchStateVec(vector<bool> &v) const;
	void WriteFev(FILE *f, const string &Name) const;

public:
	static void GetColMaps(const string &UnalignedSeq, const string &AlignedSeq,
	  vector<uint> &UnalignedPosToAlignedCol,
	  vector<uint> &AlignedColToUnalginedPos);
	};

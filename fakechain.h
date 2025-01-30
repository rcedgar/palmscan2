#pragma once

#include <set>
#include <map>
#include "pt.h"
class PDBChain;

class FakeChain
	{
public:
	PDBChain m_Chain;
	double m_Size = 5;
	vector<const PDBChain *> m_Frags;
	vector<coords_t> m_AppendCoordsVec;
	map<int64, vector<uint> > m_HashToPosVec;
	vector<uint> m_PosToFragIdx;
	vector<uint> m_PosToFragPos;

	set<int64> m_HashSet;

	void GetCube(coords_t c, intpt_t &Coords) const;

	bool IsOccupied(coords_t c) const;
	uint FindOverlap(const PDBChain &Chain) const;
	bool TryAppendFrag(const PDBChain &Frag);
	uint64 GetHash(intpt_t p) const;
	uint64 GetHash(int x, int y, int z) const;
	bool GetAppendCoords(coords_t &Coords) const;
	void Validate() const;
	void LogMe() const;

private:
	void AppendFrag(const PDBChain &Frag, coords_t AppendCoords);
	};

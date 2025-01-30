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
	vector<uint> m_PosToFragIdx;
	vector<uint> m_PosToFragPos;
	double m_MinNENDist = 4;
	double m_CADist = 3.81;

	coords_t NextCoords(const coords_t cterm2,
						const coords_t cterm1, 
						const coords_t cterm,
						double theta_bc, double theta_vc) const;
	bool IsOccupied(coords_t c, uint &Pos) const;
	uint FindOverlap(const PDBChain &Chain) const;
	bool TryAppendFrag(const PDBChain &Frag);
	bool GetAppendCoords(coords_t &Coords) const;
	void Validate() const;
	void LogMe() const;
	bool MakeFake(const vector<PDBChain *> &Frags,
				  uint L, PDBChain &Fake);
	bool AppendBest(const vector<PDBChain *> &Frags, uint Iters);
	void GetNEN_Plus(coords_t p, uint &Pos, double &Dist) const;
	void GetNEN_Minus(coords_t p, uint &Pos, double &Dist) const;
	coords_t GetCoords(int Pos) const { return m_Chain.GetCoords(Pos); }

private:
	void AppendFrag(const PDBChain &Frag, coords_t AppendCoords);
	};

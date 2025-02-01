#pragma once

#include <set>
#include <map>
#include "pt.h"
class PDBChain;

class FakeChain
	{
public:
	FakeChain() {}
	~FakeChain();

public:
	double m_MinNENDist = 4;
	double m_CADist = 3.81;
	PDBChain m_Chain;
	const vector<PDBChain *> *m_Library = 0;
	vector<const PDBChain *> m_Frags;
	vector<coords_t> m_AppendCoordsVec;
	vector<uint> m_LibIdxs;
	vector<double> m_Alphas;
	vector<double> m_Betas;
	vector<double> m_Gammas;
	double m_MDL = DBL_MAX;
	double m_NENMed = DBL_MAX;

public:
	void DeleteFrags();
	void Init(uint LibIdx);
	const PDBChain *CreateFrag(uint LibIdx,
							   const coords_t &AppendCoords,
							   double alpha,
							   double beta,
							   double gamma) const;
	void AppendFrag(uint LibIdx,
					const coords_t &AppendCoords,
					double alpha, double beta, double gamma);

	uint FindCollision(coords_t c, uint Lo, uint Hi) const;
	double FitOk(const PDBChain &Frag,
			   uint &CollisionFakePos,
			   uint &CollisionFragPos) const;
	bool GetAppendCoords(coords_t &Coords) const;
	void Validate() const;
	void LogMe() const;
	bool MakeFake(uint L);
	bool BestFit(uint Iters,
				 uint &BestLibIdx,
				 double &BestAlpha,
				 double &BestBeta,
				 double &BestGamma,
				 coords_t &BestAppendCoords) const;
	void GetNEN_Plus(uint Pos, uint &NENPos, double &Dist) const;
	void GetNEN_Minus(uint Pos, uint &NENPos, double &Dist) const;
	coords_t GetCoords(int Pos) const { return m_Chain.GetCoords(Pos); }
	double GetQualityScore(const PDBChain &Chain) const;
	double GetQualityScoreFrag(const PDBChain &Frag) const;
	void ToTsv(FILE *f) const;
	void LoadFrag(const string &LoadDir, uint FragIdx, PDBChain &Frag) const;
	void BuildPDB(const string &LoadDir,
				  PDBChain &Chain, PDBChain &ChainShuffledCA) const;
	};

#pragma once

class PDBChain;

class XProf
	{
public:
	const PDBChain *m_Chain;
	uint m_L = 0;

public:
	void Init(const PDBChain &Chain);
	void ToTsv(FILE *f) const;
	void PosToTsv(FILE *f, uint Pos) const;
	uint GetSeqLength() const { return m_Chain->GetSeqLength(); }
	uint GetFeatureCount() const;
	double GetFeature(uint Pos, uint FeatureIndex) const;
	double GetFeature0(uint Pos) const;
	double GetFeature1(uint Pos) const;
	double GetFeature2(uint Pos) const;
	double GetFeature3(uint Pos) const;
	double GetFeature4(uint Pos) const;
	double GetFeature5(uint Pos) const;
	double GetFeature6(uint Pos) const;
	double GetFeature7(uint Pos) const;

	double GetAngle(
	  uint PosA1, uint PosA2, 
	  uint PosB1, uint PosB2) const;
	uint GetSphereNr(uint Pos, double Radius) const;
	};

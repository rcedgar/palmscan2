#ifndef xprof_h
#define xprof_h

#include "myutils.h"
#include "pdbchain.h"

const uint XBINS = 10;
const uint XFEATS = 6;

class XProf
	{
public:
	const PDBChain *m_Chain = 0;
	uint m_L = 0;

	vector<double> m_NUDX_ScaledValues;
	string m_SS;

public:
	static vector<vector<double> > g_BinLos;
	static vector<vector<double> > g_Scores;

public:
	void Init(const PDBChain &Chain);
	void ToCfv(FILE *f) const;
	void PosToCfv(FILE *f, uint Pos) const;

	uint GetSeqLength() const { return m_Chain->GetSeqLength(); }

	double GetFeatureValue(uint FeatureIndex, uint Pos) const;
	void GetFeature(uint FeatureIndex, uint Pos,
	  double &Value, uint &iValue) const;
	void GetFeatures(uint Pos, vector<double> &Values) const;

	void Get_NU(uint Pos, double &Value, uint &iValue) const;
	void Get_ND(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang_m1_p1(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang_m2_p2(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang_m3_p3(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang01_23(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang01_34(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang01_45(uint Pos, double &Value, uint &iValue) const;
	void Get_ED_p4(uint Pos, double &Value, uint &iValue) const;
	void Get_ED_m4(uint Pos, double &Value, uint &iValue) const;

	uint GetFeatureX(uint FeatureIndex, uint Pos);
	double Get_NUDX(uint Pos);
	char Get_SSX(uint Pos);
	char Get_SSX2(uint Pos);
	double Get_SSD2(uint Pos);
	double Get_SSAngle2(uint Pos);
	void Get_NUDX_Lo(uint Pos, double &NU, double &ND) const;
	void Set_NUDXVec();

	uint IntScale(double Value, double MinVal, double MedVal) const;

	double GetAngle(
	  uint PosA1, uint PosA2, 
	  uint PosB1, uint PosB2) const;
	uint GetSphereNr(uint Pos, double Radius) const;

public:
	static uint GetFeatureCount();
	static const char *GetFeatureName(uint FeatureIndex);
	static uint GetFeatureIndex(const string &FeatureName);
	static void InitScoreTable();
	static uint GetFeatureBin(uint FeatureIndex, double Diff);
	static double GetScore_Bins(const vector<uint> &Bins);
	static double GetScore(const vector<double> &Features1,
	  const vector<double> &Features2);
	static double GetDiff(uint FeatureIndex, double Value1,
	  double Value2);
	static uint Get3di(const vector<double> &FeatureValues);

private:
	static void FeatureScoreBin(const string &FeatureName, uint Bin,
	  double BinLo, double Score);
	};

#endif // xprof_h

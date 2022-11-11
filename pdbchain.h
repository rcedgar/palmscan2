#pragma once

class PDBChain
	{
public:
	string m_ChainLabel;
	char m_Chain = 0;
	string m_Seq;
	vector<double> m_Xs;
	vector<double> m_Ys;
	vector<double> m_Zs;
	vector<uint> m_MotifPosVec;
	string m_SS;

public:
	void Clear()
		{
		m_ChainLabel.clear();
		m_Chain = 0;
		m_Seq.clear();
		m_Xs.clear();
		m_Ys.clear();
		m_Zs.clear();
		m_MotifPosVec.clear();
		}

	void GetLabel(string &Label) const;
	uint GetSeqLength() const;
	void FromLines(const string &Label, char ChainChar,
	  const vector<string> &Lines);
	void FromCalLines(const vector<string> &Lines);
	void ParseCalLabelLine(const string &Line);
	void ToCal(FILE *f) const;
	void ToCalSeg(FILE *f, uint Pos, uint n) const;
	void ToCal(const string &FileName) const;
	void ToPDB(const string &FileName) const;
	void CopyTriForm(const vector<double> &t,
	  const vector<vector<double> > &R,
	  PDBChain &XChain) const;
	void LogMe(bool WithCoords = false) const;
	const char *GetMotifSeq(uint MotifIndex, string &s) const;
	const char *GetMotifSeqNoFail(uint MotifIndex, string &s) const;
	void GetXYZ(uint Pos, double &x, double &y, double &z) const;
	void GetPt(uint Pos, vector<double> &Pt) const;
	void SetPt(uint Pos, const vector<double> &Pt);
	double GetDist(uint Pos1, uint Pos2) const;
	double GetDist2(uint Pos1, uint Pos2) const;
	void GetMotifCoords(vector<vector<double> > &MotifCoords) const;
	void GetMotifDists(double &AB, double &BC, double &AC) const;
	void GetMotifDists2(double &AB, double &BC, double &AC) const;
	void GetSubSeq(uint StartPos, uint n,
	  bool FailOnOverflow, string &MotifSeq) const;
	void GetSS(string &SS) const;
	void SetSS()
		{
		GetSS(m_SS);
		}

public:
	static void ReadChainsFromLines(const string &Label,
	  const vector<string> &Lines, vector<PDBChain *> &Chains);
	static void ReadChainsFromFile(const string &FileName,
	  vector<PDBChain *> &Chains);
	};

#pragma once

class PDBChain
	{
public:
	string m_Label;
	string m_Seq;
	vector<double> m_Xs;
	vector<double> m_Ys;
	vector<double> m_Zs;
	vector<vector<string> > m_ATOMs;
	vector<uint> m_MotifPosVec;
	string m_SS;

public:
	void Clear()
		{
		m_Label.clear();
		m_Seq.clear();
		m_Xs.clear();
		m_Ys.clear();
		m_Zs.clear();
		m_MotifPosVec.clear();
		m_ATOMs.clear();
		}

	uint GetSeqLength() const;
	void FromPDBLines(const string &Label,
	  const vector<string> &Lines, bool SaveAtoms);
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
	void GetMotifSeq(uint MotifIndex, string &s) const;
	void GetSubSeq(uint Pos, uint n, string &s) const;
	void GetXYZ(uint Pos, double &x, double &y, double &z) const;
	void GetPt(uint Pos, vector<double> &Pt) const;
	void SetPt(uint Pos, const vector<double> &Pt);
	double GetDist(uint Pos1, uint Pos2) const;
	double GetDist2(uint Pos1, uint Pos2) const;
	void GetMotifCoords(vector<vector<double> > &MotifCoords) const;
	void GetMotifDists(double &AB, double &BC, double &AC) const;
	void GetMotifDists2(double &AB, double &BC, double &AC) const;
	uint GetMotifPos(uint MotifIndex) const;
	void GetSubSeq(uint StartPos, uint n,
	  bool FailOnOverflow, string &MotifSeq) const;
	void XFormATOM(string &ATOM, const vector<double> &t,
	  const vector<vector<double> > &R) const;
	void GetSS(string &SS) const;
	void SetSS()
		{
		GetSS(m_SS);
		}
	bool CheckMotifCoords(bool FailOnError = true) const;
	bool CheckPPCMotifCoords(bool FailOnError = true) const;
	void GetPPC(uint PosA, uint PosB, uint PosC, PDBChain &PPC) const;

public:
	static uint GetMotifLength(uint MotifIndex);
	static uint GetPalmPrintLength(uint PosA, uint PosC, uint L);
	static void ChainsFromLines(const string &Label,
	  const vector<string> &Lines, vector<PDBChain *> &Chains,
	  bool SaveAstoms);
	static void ReadChainsFromFile(const string &FileName,
	  vector<PDBChain *> &Chains, bool SaveAtoms);
	static void AppendChainToLabel(string &Label, char Chain);
	static char GetChainCharFromPDBAtomLine(const string &Line);
	static bool IsPDBAtomLine(const string &Line);
	};

void ReadChains(const string &FileName,
  vector<PDBChain *> &Structures, bool SaveAtoms = false);
void ReadChains(const vector<string> &FileNames,
  vector<PDBChain *> &Structures, bool SaveAtoms);
void GetLabelFromFileName(const string &FileName, string &Label);

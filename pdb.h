#pragma once

class PDB
	{
public:
	string m_Label;
	string m_Seq;
	vector<double> m_Xs;
	vector<double> m_Ys;
	vector<double> m_Zs;
	vector<uint> m_MotifPosVec;

public:
	void Clear()
		{
		m_Xs.clear();
		m_Ys.clear();
		m_Zs.clear();
		m_MotifPosVec.clear();
		}

	uint GetSeqLength() const;
	void FromFile(const string &FileName);
	void ToCal(FILE *f) const;
	void ToCal(const string &FileName) const;
	void ToPDB(const string &FileName) const;
	void LogMe(bool WithCoords = false) const;
	const char *GetMotifSeq(uint MotifIndex, string &s) const;
	void GetXYZ(uint Pos, double &x, double &y, double &z) const;
	void GetPt(uint Pos, vector<double> &Pt) const;
	void SetPt(uint Pos, const vector<double> &Pt);
	double GetDist(uint Pos1, uint Pos2) const;
	void GetMotifCoords(vector<vector<double> > &MotifCoords) const;
	void GetMotifDists(double &AB, double &BC, double &AC) const;
	void GetMotifSeq(uint MotifStartPos, uint MotifLength,
	  bool FailOnOverflow, string &MotifSeq) const;
	void GetSS(uint StartPos, uint n, string &ss) const;
	};

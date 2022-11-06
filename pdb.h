#pragma once

struct TriParams
	{
	double LAB;
	double LBC;
	double LAC;
	double Radius;
	uint NABmin;
	uint NABmax;
	uint NBCmin;
	uint NBCmax;
	uint NACmin;
	uint NACmax;
	};

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
	void LogMe() const;
	const char *GetMotifSeq(uint MotifIndex, string &s) const;
	void GetXYZ(uint Pos, double &x, double &y, double &z) const;
	double GetDist(uint Pos1, uint Pos2) const;
	void GetMotifTriangle(double &LAB, double &LBC, double &LAC) const;
	void SearchTriangle(const TriParams &TP,
	  vector<uint> &PosAs,
	  vector<uint> &PosBs,
	  vector<uint> &PosCs,
	  vector<double> &RMSDs);
	void GetMotifSeqFromMidPos(uint MidPos, uint MotifLength,
	  bool FailOnOverflow, string &MotifSeq) const;
	};

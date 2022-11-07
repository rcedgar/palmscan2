#pragma once

class PDB;

class TriSearcher
	{
public:
	double MaxRMSD = 0;
	double Radius = 0;
	uint NABmin = 0;
	uint NABmax = 0;
	uint NBCmin = 0;
	uint NBCmax = 0;
	uint NACmin = 0;
	uint NACmax = 0;

	const PDB *m_Query = 0;
	const PDB *m_Ref = 0;
	vector<uint> m_PosAs;
	vector<uint> m_PosBs;
	vector<uint> m_PosCs;
	vector<double> m_RMSDs;

public:
	void Clear()
		{
		m_Query = 0;
		m_Ref = 0;
		m_PosAs.clear();
		m_PosBs.clear();
		m_PosCs.clear();
		m_RMSDs.clear();
		}

	void LogMe(FILE *f = g_fLog) const;
	void Search(const PDB &Query, const PDB &Ref);
	void GetOrder(vector<uint> &Order) const;
	};

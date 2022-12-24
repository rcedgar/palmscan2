#include "myutils.h"
#include "cmsearcher.h"

vector<const CMP *> CMSearcher::m_Profs;

void CMSearcher::Search(PDBChain &Query)
	{
	ClearSearch();
	m_Query = &Query;
	const uint RefCount = SIZE(m_Profs);
	for (uint RefIndex = 0; RefIndex < RefCount; ++RefIndex)
		{
		m_PS.m_Prof = m_Profs[RefIndex];
		m_PS.Search(Query);

		uint PosA, PosB, PosC;
		double Score = m_PS.GetPSSMStarts(PosA, PosB, PosC);
		Log("Score = %8.3f  %4u  %4u  %4u  >%s\n",
		  Score, PosA, PosB, PosC, m_Profs[RefIndex]->m_Label.c_str());//@@
		if (Score < m_Score)
			{
			m_Score = Score;
			m_PosA = PosA;
			m_PosB = PosB;
			m_PosC = PosC;
			}
		}
	}

bool CMSearcher::MeansFromFile(FILE *f, string &Label,
  vector<vector<double> > &Means)
	{
	Means.resize(PPSPL);
	for (uint i = 0; i < PPSPL; ++i)
		Means[i].resize(PPSPL, DBL_MAX);

	string Line;
	vector<string> Fields;
	for (uint i = 1; i < PPSPL; ++i)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			{
			if (i == 1)
				return false;
			Die("Premature EOF in CMP file");
			}
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == i+3);
		if (i == 1)
			Label = Fields[0];
		else
			asserta(Fields[0] == Label);
		asserta(StrToUint(Fields[1]) == i);
		for (uint j = 0; j <= i; ++j)
			{
			double Mean = StrToFloat(Fields[j+2]);
			if (i == j)
				asserta(Mean == 0);
			Means[i][j] = Mean;
			Means[j][i] = Mean;;
			}
		}
	return true;
	}

void CMSearcher::ProfsFromFile(const string &FileName)
	{
	asserta(FileName != "");
	m_Profs.clear();
	FILE *f = OpenStdioFile(FileName);
	for (;;)
		{
		CMP *Prof = new CMP;
		bool Ok = MeansFromFile(f, Prof->m_Label, Prof->m_Means);
		if (!Ok)
			{
			delete Prof;
			break;
			}
		Prof->m_StdDevs.clear();
		m_Profs.push_back(Prof);
		}
	CloseStdioFile(f);
	}

double CMSearcher::GetPSSMStarts(uint &PosA, uint &PosB, uint &PosC) const
	{
	PosA = m_PosA;
	PosB = m_PosB;
	PosC = m_PosC;
	return m_Score;
	}

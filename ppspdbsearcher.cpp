#include "myutils.h"
#include "ppspdbsearcher.h"

void PPSPDBSearcher::Search(PDBChain &Query)
	{
	ClearSearch();
	const uint RefCount = SIZE(m_Profs);
	for (uint RefIndex = 0; RefIndex < RefCount; ++RefIndex)
		{
		m_PS.m_Prof = m_Profs[RefIndex];
		m_PS.Search(Query);

		uint PosA, PosB, PosC;
		double Score = m_PS.GetPSSMStarts(PosA, PosB, PosC);
		if (Score > m_Score)
			{
			m_Score = Score;
			m_PosA = PosA;
			m_PosB = PosB;
			m_PosC = PosC;
			}
		}
	}

void PPSPDBSearcher::ProfsFromFile(const string &FileName)
	{
	asserta(FileName != "");
	m_Profs.clear();
	FILE *f = OpenStdioFile(FileName);
	for (;;)
		{
		PPSP *Prof = new PPSP;
		bool Ok = Prof->FromFile(f);
		if (!Ok)
			{
			delete Prof;
			break;
			}
		m_Profs.push_back(Prof);
		}
	CloseStdioFile(f);
	}

#include "myutils.h"
#include "pdbchain.h"
#include "calreader.h"

void PDBChain::FromCalLines(const vector<string> &Lines)
	{
	Clear();

	if (Lines.empty())
		return;

	const string &FirstLine = Lines[0];
	if (FirstLine.size() < 2 || FirstLine[0] != '>')
		Die("Invalid first line of CAL record '%s'",
		  FirstLine.c_str());

	m_ChainLabel = FirstLine.substr(1, string::npos);
	size_t n = m_ChainLabel.size();
	m_Chain = 0;
	if (n > 2)
		if (m_ChainLabel[n-2] == '_' || m_ChainLabel[n-2] == '.')
			{
			m_Chain = m_ChainLabel[n-1];
			m_ChainLabel = m_ChainLabel.substr(0, n-2);
			}

/***
>102l
M       43.619  -1.924  8.869
N       40.445  -0.876  10.670
I       38.254  2.240   11.220
F       40.340  3.621   14.036
***/
	const uint N = SIZE(Lines);
	vector<string> Fields;
	for (uint LineNr = 1; LineNr < N; ++LineNr)
		{
		const string &Line = Lines[LineNr];
		Split(Line, Fields, '\t');
		if (Fields.size() != 4 || Fields[0].size() != 1)
			Die("Invalid .cal record '%s'", Line.c_str());

		char aa = Fields[0][0];
		double X = StrToFloat(Fields[1]);
		double Y = StrToFloat(Fields[2]);
		double Z = StrToFloat(Fields[3]);

		m_Seq.push_back(aa);
		m_Xs.push_back(X);
		m_Ys.push_back(Y);
		m_Zs.push_back(Z);
		}
	}

void cmd_calinfo()
	{
	const string &FileName = opt_calinfo;

	CalReader CR;
	CR.Open(FileName);
	PDBChain Chain;
	while (CR.GetNext(Chain))
		{
		Log("%s\n", Chain.m_ChainLabel.c_str());
		}
	}

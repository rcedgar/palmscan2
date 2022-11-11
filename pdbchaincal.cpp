#include "myutils.h"
#include "pdbchain.h"
#include "calreader.h"
#include "abcxyz.h"

//	>AIZ00432 A:1:ISGDNTKWGPIH B:119:QGIHHATSSLLTSL C:160:GSSDDYAK
void PDBChain::ParseCalLabelLine(const string &Line)
	{
	if (Line[0] != '>')
		Die("Invalid .cal/.ppc file, expected '>' in line: %s",
		  Line.c_str());

	vector<string> Fields;
	Split(Line, Fields, ' ');
	uint PosA = UINT_MAX;
	uint PosB = UINT_MAX;
	uint PosC = UINT_MAX;
	if (Fields.size() > 1)
		{
		if (Fields.size() == 4 &&
		  StartsWith(Fields[1], "A:") &&
		  StartsWith(Fields[2], "B:") &&
		  StartsWith(Fields[3], "C:"))
			{
			m_ChainLabel = Fields[0].substr(1, string::npos);
			PosA = (uint) atoi(Fields[1].c_str() + 2);
			PosB = (uint) atoi(Fields[2].c_str() + 2);
			PosC = (uint) atoi(Fields[3].c_str() + 2);
			asserta(PosA > 0);
			asserta(PosB > 0);
			asserta(PosC > 0);
			asserta(PosA < PosB);
			asserta(PosB < PosC);

			m_MotifPosVec.push_back(PosA - 1);
			m_MotifPosVec.push_back(PosB - 1);
			m_MotifPosVec.push_back(PosC - 1);
			}
		else
			Die("Invalid .ppc label %s", Line.c_str());
		}
	else
		m_ChainLabel = Line.substr(1, string::npos);

	size_t n = m_ChainLabel.size();
	m_Chain = 0;
	if (n > 2)
		if (m_ChainLabel[n-2] == '_' || m_ChainLabel[n-2] == '.')
			{
			m_Chain = m_ChainLabel[n-1];
			m_ChainLabel = m_ChainLabel.substr(0, n-2);
			}
	}

void PDBChain::FromCalLines(const vector<string> &Lines)
	{
	Clear();

	if (Lines.empty())
		return;

	const string &FirstLine = Lines[0];
	ParseCalLabelLine(FirstLine);

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

	if (m_MotifPosVec.size() == 3)
		{
		uint QL = SIZE(m_Seq);
		uint PosA = m_MotifPosVec[0];
		uint PosB = m_MotifPosVec[1];
		uint PosC = m_MotifPosVec[2];
		asserta(PosA == 0);
		asserta(PosB > PosA + AL && PosB + BL < PosC);
		asserta(PosC + CL == QL);
		}
	else
		asserta(m_MotifPosVec.empty());
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

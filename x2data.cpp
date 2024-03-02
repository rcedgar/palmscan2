#include "myutils.h"
#include "xprof.h"
#include "x2data.h"
#include <set>

/***
Input is created by xfeatures

palmscan2 \
  -xcluster d:/int/scop40/out/xfeatures.tsv
  -log xcluster.log

# head -n5 /d/int/dave_grant/scop40/aligned_xfeatvecs_0.6_0.8.tsv | columns.py
Q  Q_Ang_m2_p2  Q_Ang_m3_p3  Q_ED_p4  Q_ED_m4  Q_NU  Q_ND  R  R_Ang_m2_p2  R_Ang_m3_p3  R_ED_p4  R_ED_m4  R_NU  R_ND    QDom  QPos    RDom  RPos
N         83.1         62.8     9.57     11.7     0    14  L         13.9           22     11.6        0     1    19  d1a04a1    4  d1fsea_    3
L         72.9         48.9      5.6     8.06    19     4  L         18.5         10.6     5.74     12.9    21     5  d1a04a1    6  d1fsea_    4
T         34.4         74.9     6.27     9.21     2    16  T         34.9         83.4     6.27       13     3    13  d1a04a1    7  d1fsea_    5
P          125          107     6.25     9.57     0    15  K          125          109      6.2     13.1     1    13  d1a04a1    8  d1fsea_    6

***/

void X2Data::FromTsv(const string &FileName)
	{
	m_Aminos.clear();
	m_Aminos1.clear();
	m_Aminos2.clear();
	m_FeatureValuesVec.clear();
	m_FeatureValuesVec1.clear();
	m_FeatureValuesVec2.clear();

	FILE *f = OpenStdioFile(FileName);
	string HdrLine;
	bool Ok = ReadLineStdioFile(f, HdrLine);
	vector<string> HdrFields;
	Split(HdrLine, HdrFields, '\t');
	asserta(SIZE(HdrFields) == 2*XFEATS + 6);
	asserta(HdrFields[0] == "Q");
	asserta(HdrFields[XFEATS+1] == "R");
	for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
		{
		asserta(HdrFields[FeatureIndex+1] == 
		  "Q_" +  (string) XProf::GetFeatureName(FeatureIndex));
		asserta(HdrFields[XFEATS+2+FeatureIndex] == 
		  "R_" + (string) XProf::GetFeatureName(FeatureIndex));
		}

	set<pair<string, uint> > DoneSet;

	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2*XFEATS + 6);

		{
		asserta(SIZE(Fields[0]) == 1);
		char Amino = Fields[0][0];
		vector<double> Values;
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			double Value = StrToFloat(Fields[FeatureIndex+1]);
			Values.push_back(Value);
			}

		m_Aminos1.push_back(Amino);
		m_FeatureValuesVec1.push_back(Values);

		string Dom = Fields[2*XFEATS+2];
		uint Coord = StrToUint(Fields[2*XFEATS+3]);
		pair<string, uint> Pair(Dom, Coord);
		if (DoneSet.find(Pair) == DoneSet.end())
			{
			m_Aminos.push_back(Amino);
			m_FeatureValuesVec.push_back(Values);
			m_Doms.push_back(Dom);
			m_Coords.push_back(Coord);
			DoneSet.insert(Pair);
			}
		}

		{
		asserta(SIZE(Fields[XFEATS+1]) == 1);
		char Amino = Fields[XFEATS+1][0];
		vector<double> Values;
		for (uint FeatureIndex = 0; FeatureIndex < XFEATS; ++FeatureIndex)
			{
			double Value = StrToFloat(Fields[XFEATS+2+FeatureIndex]);
			Values.push_back(Value);
			}

		m_Aminos2.push_back(Amino);
		m_FeatureValuesVec2.push_back(Values);

		string Dom = Fields[2*XFEATS+4];
		uint Coord = StrToUint(Fields[2*XFEATS+5]);
		pair<string, uint> Pair(Dom, Coord);
		if (DoneSet.find(Pair) == DoneSet.end())
			{
			m_Aminos.push_back(Amino);
			m_FeatureValuesVec.push_back(Values);
			m_Doms.push_back(Dom);
			m_Coords.push_back(Coord);
			DoneSet.insert(Pair);
			}
		}
		}
	const uint N = SIZE(m_FeatureValuesVec);
	const uint M = SIZE(m_FeatureValuesVec1);
	asserta(SIZE(m_Aminos) == N);
	asserta(SIZE(m_Doms) == N);
	asserta(SIZE(m_Coords) == N);

	asserta(SIZE(m_Aminos1) == M);
	asserta(SIZE(m_Aminos2) == M);
	asserta(SIZE(m_FeatureValuesVec2) == M);

	ProgressLog("%s unique feature vecs\n", IntToStr(N));
	ProgressLog("%s total feature vecs\n", IntToStr(M));
	CloseStdioFile(f);
	}

#include "myutils.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "cddata.h"
#include "cmp.h"
#include <map>

static void ReadMotifDataCoords(
  CDInfo &Info,
  vector<vector<uint> > &MotifCoordsVec,
  vector<vector<string> > &MotifSeqsVec,
  vector<string> &Labels,
  map<string, uint> &LabelToIndex)
	{
	MotifCoordsVec.clear();
	MotifSeqsVec.clear();
	Labels.clear();
	LabelToIndex.clear();

	vector<string> MotifNames;
	vector<uint> MotifLengths;

	MotifNames.push_back("A");
	MotifNames.push_back("B");
	MotifNames.push_back("C");
	MotifNames.push_back("D");
	MotifNames.push_back("E");
	MotifNames.push_back("F2");

	MotifLengths.push_back(12);
	MotifLengths.push_back(14);
	MotifLengths.push_back(10);
	MotifLengths.push_back(7);
	MotifLengths.push_back(7);
	MotifLengths.push_back(7);

	const uint NM = 6;

	Info.Init(MotifNames, MotifLengths);

	if (!optset_motif_coords)
		Die("Must specify -motif_coords");

	FILE *f = OpenStdioFile(opt_motif_coords);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2*NM + 1);

		vector<uint> MotifCoords;
		vector<string> MotifSeqs;

		const string &Label = Fields[0];

		for (uint MotifIndex = 0; MotifIndex < NM; ++MotifIndex)
			{
			const string &s = Fields[1 + 2*MotifIndex];
			if (s != ".")
				{
				uint Pos = StrToUint(s);
				asserta(Pos > 0);
				MotifCoords.push_back(Pos-1);
				}
			else
				MotifCoords.push_back(UINT_MAX);
			MotifSeqs.push_back(Fields[2 + 2*MotifIndex]);
			}

		uint Index = SIZE(Labels);
		asserta(LabelToIndex.find(Label) == LabelToIndex.end());
		LabelToIndex[Label] = Index;
		Labels.push_back(Label);
		MotifCoordsVec.push_back(MotifCoords);
		MotifSeqsVec.push_back(MotifSeqs);
		}
	}

static void GetData(const PDBChain &Q, const CDInfo &Info,
  const vector<uint> &MotifCoords, 
  const vector<string> &MotifSeqs,
  CDData &Data)
	{
	Log("\n_________________________\n");
	Log(">%s\n", Q.m_Label.c_str());

	const uint NM = Info.m_MotifCount;
	for (uint MotifIndex = 0; MotifIndex < NM; ++MotifIndex)
		{
		uint ML = Info.m_MotifLengths[MotifIndex];
		uint Pos = MotifCoords[MotifIndex];
		if (Pos == UINT_MAX)
			continue;
		const string &MotifSeq = MotifSeqs[MotifIndex];
		string Motif;
		Q.GetSubSeq(Pos, ML, Motif);
		Log(" %3.3s", Info.m_MotifNames[MotifIndex].c_str());
		Log(" %4u", Pos + 1);
		Log("  %16.16s", MotifSeq.c_str());
		Log("  %16.16s", Motif.c_str());
		Log("\n");
		asserta(Motif == MotifSeq);
		}

	const uint Size = Info.GetSize();
	for (uint Ix1 = 0; Ix1 < Size; ++Ix1)
		{
		uint MotifIndex1 = Info.m_MotifIndexes[Ix1];
		uint Coord1 = Info.m_MotifCoords[Ix1];
		for (uint Ix2 = 0; Ix2 < Size; ++Ix2)
			{
			uint MotifIndex2 = Info.m_MotifIndexes[Ix2];
			uint Coord2 = Info.m_MotifCoords[Ix2];
			}
		}
	}

void cmd_cdp_train()
	{
	const string &QueryFN = opt_cdp_train;

	CDInfo Info;
	vector<vector<string> > MotifSeqsVec;
	vector<vector<uint> > MotifCoordsVec;
	vector<string> Labels;
	map<string, uint> LabelToIndex;
	ReadMotifDataCoords(Info, MotifCoordsVec, 
	  MotifSeqsVec, Labels, LabelToIndex);

	PDBChain Q;
	ChainReader CR;
	CR.Open(QueryFN);

	CDData Data;
	vector<CDData> DataVec;
	while (CR.GetNext(Q))
		{
		const string &Label = Q.m_Label;
		map<string, uint>::const_iterator p = LabelToIndex.find(Label);
		if (p == LabelToIndex.end())
			{
			Log("Not found >%s\n", Label.c_str());
			continue;
			}

		uint Index = p->second;
		const vector<uint> &MotifCoords = MotifCoordsVec[Index];
		const vector<string> &MotifSeqs = MotifSeqsVec[Index];

		GetData(Q, Info, MotifCoords, MotifSeqs, Data);
		}
	}

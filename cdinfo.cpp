#include "myutils.h"
#include "cdinfo.h"

void CDInfo::LogMe() const
	{
	Log("%u motifs: \n", m_MotifCount);
	for (uint i = 0; i < m_MotifCount; ++i)
		Log(" %s(%s)",
		  m_MotifNames[i].c_str(),
		  m_MotifLengths[i]);
	}

void CDInfo::Init(const vector<string> &MotifNames,
	vector<uint> &MotifLengths)
	{
	Clear();
	m_MotifCount = SIZE(MotifNames);
	asserta(SIZE(MotifLengths) == m_MotifCount);
	m_MotifNames = MotifNames;
	m_MotifLengths = MotifLengths;
	uint Offset = 0;
	for (uint i = 0; i < m_MotifCount; ++i)
		{
		uint ML = m_MotifLengths[i];
		m_Offsets.push_back(ML);
		Offset += ML;
		}

	uint Size = GetSize();
	for (uint Ix = 0; Ix < Size; ++Ix)
		{
		uint MotifIndex, i;
		GetCoord(Ix, MotifIndex, i);
		m_MotifIndexes.push_back(MotifIndex);
		m_MotifCoords.push_back(i);
		}
	}

void CDInfo::GetCoord(uint Ix, uint &MotifIndex, uint &i) const
	{
	for (MotifIndex = 0; ; ++MotifIndex)
		{
		if (MotifIndex + 1 == m_MotifCount ||
		  Ix < m_Offsets[MotifIndex+1])
			{
			i = Ix - m_Offsets[MotifIndex];
			return;
			}
		}
	asserta(false);
	}

uint CDInfo::GetIx(uint MotifIndex, uint i) const
	{
	assert(MotifIndex < SIZE(m_Offsets));
	assert(i < m_MotifLengths[i]);
	uint Ix = m_Offsets[MotifIndex] + i;
	return Ix;
	}

void CDInfo::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, "cdinfo\t%u", m_MotifCount);
	for (uint i = 0; i < m_MotifCount; ++i)
		fprintf(f, "%u\t%u\t%s\n",
		  i,
		  m_MotifLengths[i],
		  m_MotifNames[i].c_str());
	}

void CDInfo::FromTsv(FILE *f)
	{
	Clear();
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 2);
	asserta(Fields[0] == "cdinfo");
	m_MotifCount = StrToUint(Fields[0]);
	asserta(m_MotifCount > 0 && m_MotifCount < 99);
	for (uint i = 0; i < m_MotifCount; ++i)
		{
		Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		asserta(StrToUint(Fields[0]) == i);
		uint ML = StrToUint(Fields[1]);
		asserta(ML > 0 && ML < 99);
		m_MotifLengths.push_back(ML);
		m_MotifNames.push_back(Fields[2]);
		}
	}

uint CDInfo::GetSize() const
	{
	assert(m_MotifCount > 0 && SIZE(m_Offsets) == m_MotifCount);
	uint ML = m_MotifLengths[m_MotifCount-1];
	uint Off = m_Offsets[m_MotifCount-1];
	uint Size = ML + Off;
	return Size;
	}

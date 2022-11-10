#include "myutils.h"
#include "calreader.h"
#include "omplock.h"

void CalReader::Open(const string &FileName)
	{
	Clear();
	m_f = OpenStdioFile(FileName);
	bool Ok = ReadLineStdioFile(m_f, m_PendingLine);
	if (!Ok)
		Close();
	if (m_PendingLine[0] != '>')
		Die("Invalid .cal file, does not start with '>': %s",
		  FileName.c_str());
	}

bool CalReader::GetNext(PDBChain &Chain)
	{
	if (m_f == 0)
		return false;

	Lock("CalReader::GetNext");
	if (m_EOF)
		{
		Unlock("CalReader::GetNext");
		return false;
		}
	asserta(!m_PendingLine.empty() && m_Lines.empty());
	asserta(m_PendingLine[0] == '>');

	m_Lines.push_back(m_PendingLine);
	m_PendingLine.clear();
	string Line;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(m_f, Line);
		if (!Ok)
			{
			m_EOF = true;
			break;
			}
		if (Line[0] == '>')
			{
			m_PendingLine = Line;
			break;
			}
		m_Lines.push_back(Line);
		}

	Chain.FromCalLines(m_Lines);
	m_Lines.clear();
	Unlock("CalReader::GetNext");
	return true;
	}

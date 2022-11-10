#pragma once

#include "pdbchain.h"

class CalReader
	{
public:
	FILE *m_f = 0;
	string m_PendingLine;
	vector<string> m_Lines;
	bool m_EOF = false;

public:
	void Clear()
		{
		m_f = 0;
		m_PendingLine.clear();
		m_Lines.clear();
		m_EOF = false;
		}

	void Close()
		{
		if (m_f != 0)
			CloseStdioFile(m_f);
		Clear();
		}

	void Open(const string &FileName);
	bool GetNext(PDBChain &Chain);
	};

#pragma once

#include "calreader.h"

enum CR_TYPE
	{
	CR_None,
	CR_PDB,
	CR_CAL,
	CR_Files,
	};

class ChainReader
	{
public:
	CR_TYPE m_Type = CR_None;
	string m_FileName;
	vector<string> m_FileNames;
	vector<string> m_Lines;
	vector<PDBChain *> m_FilesChains;
	uint m_EndOfInput = true;
	CalReader m_CR;
	bool m_SaveAtoms = false;
	uint m_FileIndex = UINT_MAX;
	uint m_LineIndex = UINT_MAX;
	uint m_FilesChainIndex = UINT_MAX;

public:
	void Clear()
		{
		m_Type = CR_None;
		m_FileName.clear();
		m_FileNames.clear();
		m_Lines.clear();
		m_FilesChains.clear();
		m_EndOfInput = true;
		m_CR.Clear();
		m_SaveAtoms = false;
		m_FileIndex = UINT_MAX;
		m_LineIndex = UINT_MAX;
		m_FilesChainIndex = UINT_MAX;
		}

	void Open(const string &FileName, bool SaveAtoms);
	bool GetNext(PDBChain &Chain);
	bool GetNext_CAL(PDBChain &Chain);
	bool GetNext_PDB(PDBChain &Chain);
	bool GetNext_Files(PDBChain &Chain);
	double GetPctDone() const;
	double GetPctDone_CAL() const;
	double GetPctDone_PDB() const;
	double GetPctDone_Files() const;
	void GetStrPctDone(string &s) const;

public:
	static bool KeepPDBAtomLine(const string &Line);
	};

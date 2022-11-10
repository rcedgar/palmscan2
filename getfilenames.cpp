#include "myutils.h"
#include "pdbchain.h"

void GetFileNames(const string &SpecFileName, vector<string> &FileNames)
	{
	FileNames.clear();
	if (EndsWith(SpecFileName, ".files"))
		{
		FILE *f = OpenStdioFile(SpecFileName);
		string FileName;
		while (ReadLineStdioFile(f, FileName))
			{
			StripWhiteSpace(FileName);
			FileNames.push_back(FileName);
			}
		}
	else
		FileNames.push_back(SpecFileName);
	}

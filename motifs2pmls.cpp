#include "myutils.h"

const char *GetPymolColor(uint MotifIndex);

void cmd_motifs2pmls()
	{
	const string &TsvFN = opt_motifs2pmls;

	vector<string> Lines;
	ReadLinesFromFile(TsvFN, Lines);
	const uint LineCount = SIZE(Lines);
	asserta(LineCount > 1);
	vector<string> HdrFields;
	Split(Lines[0], HdrFields, '\t');
	asserta(!HdrFields.empty());
	const uint MotifCount = SIZE(HdrFields) - 1;
	asserta(MotifCount > 0);
	for (uint i = 1; i < LineCount; ++i)
		{
		const string &Line = Lines[i];
		vector<string> Fields;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == MotifCount + 1);
		const string &Label = Fields[0];
		string PMLFN = string(opt_output) + Label + ".pml";
		string PDBFN = string(opt_input) + Label + ".pdb";
		FILE *f = CreateStdioFile(PMLFN);
		fprintf(f, "#!/usr/bin/pymol\n");
		fprintf(f, "\n");
		fprintf(f, "cmd.load(\"%s\")\n", PDBFN.c_str());
		fprintf(f, "color gray70\n");
		fprintf(f, "bg white\n");
		for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
			{
			char M = 'A' + MotifIndex;
			const char *Color = GetPymolColor(MotifIndex);
			const char *MotifSeq = Fields[MotifIndex+1].c_str();
			fprintf(f, "select Mtf%c, pepseq %s\n", M, MotifSeq);
			fprintf(f, "color %s, Mtf%c\n", Color, M);
			}
		fprintf(f, "deselect\n");
		fprintf(f, "util.performance(0)\n");
		fprintf(f, "\n");
		CloseStdioFile(f);
		}
	}

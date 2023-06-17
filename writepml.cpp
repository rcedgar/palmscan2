#include "myutils.h"
#include "structprof.h"
#include "chainreader.h"
#include "outputfiles.h"
#include "cmpsearcher.h"
#include "gsprof.h"
#include "abcxyz.h"

extern uint g_APos;
extern uint g_BPos;
extern uint g_CPos;
extern uint g_DPos;
extern uint g_EPos;
extern uint g_F1Pos;
extern uint g_F2Pos;

void WritePML(const PDBChain &Chain, const string &PDBFileName)
	{
	if (g_fpml == 0)
		return;
	string Label;
	GetLabelFromFileName(PDBFileName, Label);
//	fprintf(g_fpml, "window hide\n");
	fprintf(g_fpml, "cmd.load(\"%s\")\n", PDBFileName.c_str());
	fprintf(g_fpml, "select %s\n", Label.c_str());
	fprintf(g_fpml, "color gray40, %s\n", Label.c_str());
	if (optset_png)
		fprintf(g_fpml, "bg_color white\n");

	if (g_APos != UINT_MAX)
		{
		uint Res1 = Chain.GetResidueNr(g_APos);
		uint Res2 = Chain.GetResidueNr(g_APos + AL - 1);
		fprintf(g_fpml, "cmd.select(\"motifA\", \"%s and resi %u-%u\")\n",
			Label.c_str(), Res1, Res2);
		fprintf(g_fpml, "color tv_blue, motifA\n");

		fprintf(g_fpml, "cmd.select(\"motifA_D\", \"%s and resi %d and name CA\")\n",
			Label.c_str(), Res1+3);
		fprintf(g_fpml, "show sphere, motifA_D\n");
		fprintf(g_fpml, "set sphere_transparency, 0.3\n");
		}
	if (g_BPos != UINT_MAX)
		{
		int Res1 = Chain.GetResidueNr(g_BPos);
		int Res2 = Chain.GetResidueNr(g_BPos + BL - 1);
		fprintf(g_fpml, "cmd.select(\"motifB\", \"%s and resi %d-%d\")\n",
			Label.c_str(), Res1, Res2);
		fprintf(g_fpml, "color tv_green, motifB\n");

		fprintf(g_fpml, "cmd.select(\"motifB_G\", \"%s and resi %d and name CA\")\n",
			Label.c_str(), Res1+1);
		fprintf(g_fpml, "show sphere, motifB_G\n");
		fprintf(g_fpml, "set sphere_transparency, 0.3\n");
		}
	if (g_CPos != UINT_MAX)
		{
		int Res1 = Chain.GetResidueNr(g_CPos);
		int Res2 = Chain.GetResidueNr(g_CPos + CL - 1);
		fprintf(g_fpml, "cmd.select(\"motifC\", \"%s and resi %d-%d\")\n",
			Label.c_str(), Res1, Res2);
		fprintf(g_fpml, "color tv_red, motifC\n");

		fprintf(g_fpml, "cmd.select(\"motifC_D\", \"%s and resi %d and name CA\")\n",
			Label.c_str(), Res1+3);
		fprintf(g_fpml, "show sphere, motifC_D\n");
		fprintf(g_fpml, "set sphere_transparency, 0.3\n");
		}
	if (g_DPos != UINT_MAX)
		{
		int Res1 = Chain.GetResidueNr(g_DPos);
		int Res2 = Chain.GetResidueNr(g_DPos + 6);
		fprintf(g_fpml, "cmd.select(\"motifD\", \"%s and resi %u-%u\")\n",
			Label.c_str(), Res1, Res2);
		fprintf(g_fpml, "color yellow, motifD\n");

		fprintf(g_fpml, "cmd.select(\"motifD_X\", \"%s and resi %d and name CA\")\n",
			Label.c_str(), Res1+3);
		fprintf(g_fpml, "show sphere, motifD_X\n");
		fprintf(g_fpml, "set sphere_transparency, 0.3\n");
		}
	if (g_EPos != UINT_MAX)
		{
		int Res1 = Chain.GetResidueNr(g_EPos);
		int Res2 = Chain.GetResidueNr(g_EPos + 6);
		fprintf(g_fpml, "cmd.select(\"motifE\", \"%s and resi %u-%u\")\n",
			Label.c_str(), Res1, Res2);
		fprintf(g_fpml, "color orange, motifE\n");

		fprintf(g_fpml, "cmd.select(\"motifE_X\", \"%s and resi %d and name CA\")\n",
			Label.c_str(), Res1+3);
		fprintf(g_fpml, "show sphere, motifE_X\n");
		fprintf(g_fpml, "set sphere_transparency, 0.3\n");
		}
	if (g_F2Pos != UINT_MAX)
		{
		int Res1 = Chain.GetResidueNr(g_F2Pos);
		int Res2 = Chain.GetResidueNr(g_F2Pos + 6);
		fprintf(g_fpml, "cmd.select(\"motifF\", \"%s and resi %u-%u\")\n",
			Label.c_str(), Res1, Res2);
		fprintf(g_fpml, "color cyan, motifF\n");

		fprintf(g_fpml, "cmd.select(\"motifF_R\", \"%s and resi %d and name CA\")\n",
			Label.c_str(), Res1+2);
		fprintf(g_fpml, "show sphere, motifF_R\n");
		fprintf(g_fpml, "set sphere_transparency, 0.3\n");
		}
	fprintf(g_fpml, "deselect\n");
	if (optset_png)
		{
		fprintf(g_fpml, "ray\n");
		fprintf(g_fpml, "png %s\n", opt_png);
		fprintf(g_fpml, "quit\n");
		}
	}

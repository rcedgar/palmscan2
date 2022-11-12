#include "myutils.h"
#include "timing.h"
#include "outputfiles.h"

int g_Frame = 0;

int main(int argc, char **argv)
	{
	MyCmdLine(argc, argv);
	LogProgramInfoAndCmdLine();
	if (!opt_quiet)
		{
		PrintProgramInfo(stdout);
		PrintCopyright(stdout);
		}

	if (optset_frame)
		{
		g_Frame = atoi(opt_frame);
		asserta(g_Frame >= -3 && g_Frame <= 3 && g_Frame != 0);
		}

	ProgressPrefix(false);
	Progress("%s\n", g_ShortCmdLine.c_str());
	ProgressPrefix(true);

	OpenOutputFiles();
	if (0)
		;
#define x(opt, f)	else if (optset_##opt) { void f(); f(); }
	x(build_pssm, BuildPSM)
	x(build_rdrp_model, BuildRdRpModel)
	x(search_pssms, cmd_search_pssms)
	x(cluster_cl, cmd_cluster_cl)
	x(alignabc, cmd_alignabc)
	x(pdbss, cmd_pdbss)
	x(pdb2cal, cmd_pdb2cal)
	x(search3d, cmd_search3d)
	x(search3d_cal, cmd_search3d_cal)
	x(search_pssms, cmd_search_pssms)
	x(search3d_pssms, cmd_search3d_pssms)
	x(search3d_ppc, cmd_search3d_ppc)
	x(remove_dupes, cmd_remove_dupes)
	x(cal2fa, cmd_cal2fa)
	x(getchains, cmd_getchains)
	x(cluster_ppc, cmd_cluster_ppc)
#undef x
	else
		Die("No command specified");

	CloseOutputFiles();

	LogElapsedTimeAndRAM();
	return 0;
	}

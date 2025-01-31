#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "fakechain.h"

void cmd_fake()
	{
	const string &InputFN = g_Arg1;

	FakeChain FC;
	ReadChains(InputFN, FC.m_Library);

	FC.m_Chain.m_Label = "FC";
	bool Ok = FC.MakeFake(150);
	asserta(Ok);

	FC.LogMe();
	FC.m_Chain.ToPDB("fake.pdb");
	}

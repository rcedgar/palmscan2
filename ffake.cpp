#include "myutils.h"
#include "chainfaker.h"

void cmd_ffake()
	{
	ChainFaker CF;
	CF.m_Trace = true;
	CF.ReadSCOP40(g_Arg1);

	PDBChain FakeChain;
	CF.MakeFake(0, FakeChain);
	}

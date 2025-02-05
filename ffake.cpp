#include "myutils.h"
#include "chainfaker.h"

void cmd_ffake()
	{
	ChainFaker CF;
	CF.m_Trace = true;
	CF.ReadSCOP40(g_Arg1);

	const uint ChainCount = SIZE(CF.m_SCOP40);
	uint ChainIdx = randu32()%ChainCount;

	PDBChain FakeChain;
	CF.MakeFake(ChainIdx, FakeChain);
	}

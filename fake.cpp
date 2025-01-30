#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "fakechain.h"

bool FakeChain::MakeFake(const vector<PDBChain *> &Frags,
						 uint L, PDBChain &Fake)
	{
	FakeChain FC;
	const uint FragCount = SIZE(Frags);
	for (uint Try = 0; Try < 100; ++Try)
		{
		uint FragIdx = randu32()%FragCount;
		const PDBChain &Frag = *Frags[FragIdx];
		FC.TryAppendFrag(Frag);
		}

	return false;
	}

void cmd_fake()
	{
	const string &InputFN = g_Arg1;

	vector<PDBChain *> Frags;
	ReadChains(InputFN, Frags);
	const uint N = SIZE(Frags);

	FakeChain FC;
	bool Ok = FC.AppendBest(Frags, 100);
	asserta(Ok);
	FC.Validate();

	Ok = FC.AppendBest(Frags, 100);
	asserta(Ok);
	FC.Validate();

	Ok = FC.AppendBest(Frags, 100);
	asserta(Ok);
	FC.Validate();

	Ok = FC.AppendBest(Frags, 100);
	asserta(Ok);
	FC.Validate();

	FC.m_Chain.ToPDB("fake.pdb");
	}

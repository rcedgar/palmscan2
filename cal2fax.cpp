#include "myutils.h"
#include "pdbchain.h"
#include "xbinner.h"
#include "xprof.h"
#include "seqdb.h"
#include "alpha.h"
#include "outputfiles.h"

void cmd_cal2fax()
	{
	const string &CalFN = opt_cal2fax;
	vector<PDBChain *> Chains;
	ReadChains(CalFN, Chains);
	const uint ChainCount = SIZE(Chains);

	XProf::InitScoreTable();
	XBinner::InitCentroids();
	XBinner XB;

	XProf XP;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Converting cal to X");
		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;

		XP.Init(Chain);
		const uint L = Chain.GetSeqLength();
		string XSeq;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			vector<double> Values;
			XP.GetFeatures(Pos, Values);
			uint XLetter = XB.GetLetter(Values);
			char XChar = g_LetterToCharAmino[XLetter];
			XSeq += XChar;
			}
		SeqToFasta(g_ffasta, Label, XSeq);
		}
	}

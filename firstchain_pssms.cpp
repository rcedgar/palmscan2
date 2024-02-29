#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "rdrpsearcher.h"
#include "outputfiles.h"

void cmd_firstchain_pssms()
	{
	const string &InputFN = opt_firstchain_pssms;

	RdRpModel Model;
	GetRdrpModel(Model);

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);

	const uint ChainCount = SIZE(Chains);
	uint ConvertedCount = 0;
	for (uint i = 0; i < ChainCount; ++i)
		{
		PDBChain &Chain = *Chains[i];
		const string &QLabel = Chain.m_Label;
		const string &QSeq = Chain.m_Seq;
		const uint QL = SIZE(QSeq);

		RdRpSearcher RS;
		RS.Init(Model);
		RS.Search(QLabel, QSeq);
		RS.WriteOutput();
		if (!RS.IsHit())
			continue;

		uint APos = RS.GetMotifPos(A);
		uint BPos = RS.GetMotifPos(B);
		uint CPos = RS.GetMotifPos(C);

		Chain.m_MotifPosVec.resize(3);
		Chain.m_MotifPosVec[0] = APos;
		Chain.m_MotifPosVec[1] = BPos;
		Chain.m_MotifPosVec[2] = CPos;

		Chain.ToPDB(opt_output);

		Chain.ToPML_Seqs(g_fpml, opt_output);
		++ConvertedCount;
		break;
		}

	if (ConvertedCount == 1)
		ProgressLog("Converted ok\n");
	else if (ConvertedCount == 0)
		ProgressLog("No PSSM hits");
	else
		asserta(false);
	}

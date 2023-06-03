#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "rdrpsearcher.h"
#include "outputfiles.h"

void GetRdrpModel(RdRpModel &Model)
	{
	if (optset_model)
		{
		const string &ModelFileName = opt_model;
		Model.FromModelFile(ModelFileName);
		}
	else
		{
		vector<string> Lines;
		Model.GetDefaultModelLines(Lines);
		Model.FromLines(Lines);
		}
	}

void cmd_palmcore_pssms()
	{
	const string &InputFN = opt_palmcore_pssms;

	RdRpModel Model;
	GetRdrpModel(Model);

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);

	const uint ChainCount = SIZE(Chains);
	uint PermutedCount = 0;
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
		if (CPos < APos)
			{
			++PermutedCount;
			continue;
			}

		++ConvertedCount;
		Chain.SetMotifPosVec(APos, BPos, CPos);

		PDBChain PC;
		Chain.GetPalmCore(PC);

		PC.ToPDB(opt_output);
		PC.ToPML_Seqs(g_fpml, opt_output);
		break;
		}

	ProgressLog("%u chains converted\n", ConvertedCount);
	if (PermutedCount > 0)
		Warning("%u permuted domains skipped", PermutedCount);
	}

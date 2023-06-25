#include "myutils.h"
#include "shapes.h"

/***
Create a new Shapes object by combining palmprint shapes
object (Shapes_ABC) with profile for one new motif.
***/
void Shapes::Train_AddMotif(const Shapes &Shapes_ABC,
  const vector<PDBChain *> &Chains,
  const vector<vector<string> > &ABCSeqsVec,
  const vector<string> &NewMotifSeqs)
	{
	const uint N = SIZE(Chains);
	asserta(SIZE(ABCSeqsVec) == N);
	asserta(SIZE(NewMotifSeqs) == N);

	bool First = true;
	bool BeforeA = false;
	uint LX = UINT_MAX;
	vector<PDBChain *> Chains2;
	vector<uint> PosXs;
	vector<uint> PosAs;
	vector<uint> PosBs;
	vector<uint> PosCs;
	for (uint i = 0; i < N; ++i)
		{
		PDBChain *Chain = Chains[i];
		const string &Seq = Chain->m_Seq;
		const string &SeqX = NewMotifSeqs[i];

		const vector<string> &ABCSeqs = ABCSeqsVec[i];
		asserta(SIZE(ABCSeqs) == 3);

		const string &SeqA = ABCSeqs[0];
		const string &SeqB = ABCSeqs[1];
		const string &SeqC = ABCSeqs[2];

		if (SeqX == "" || SeqX == ".")
			continue;
		if (SeqA == "" || SeqA == ".")
			continue;
		if (SeqB == "" || SeqB == ".")
			continue;
		if (SeqC == "" || SeqC == ".")
			continue;

		size_t stPosX = Seq.find(SeqX);
		size_t stPosA = Seq.find(SeqA);
		size_t stPosB = Seq.find(SeqB);
		size_t stPosC = Seq.find(SeqC);
		asserta(stPosX != string::npos);
		asserta(stPosA != string::npos);
		asserta(stPosB != string::npos);
		asserta(stPosC != string::npos);

		uint PosX = uint(stPosX);
		uint PosA = uint(stPosA);
		uint PosB = uint(stPosB);
		uint PosC = uint(stPosC);

		if (First)
			{
			BeforeA = (PosX < PosA);
			First = false;
			LX = SIZE(SeqX);
			}
		else
			{
			asserta(SIZE(SeqX) == LX);
			if (BeforeA)
				assert(PosX < PosA);
			else
				assert(PosX > PosA);
			}

		Chains2.push_back(Chain);
		PosXs.push_back(PosX);
		PosAs.push_back(PosA);
		PosBs.push_back(PosB);
		PosCs.push_back(PosC);
		}

	const uint N2 = SIZE(Chains2);
	if (N2 < 10)
		Die("N2=%u", N2);
	asserta(SIZE(PosXs) == N2);
	asserta(SIZE(PosAs) == N2);
	asserta(SIZE(PosBs) == N2);
	asserta(SIZE(PosCs) == N2);

	uint MinDist = UINT_MAX;
	uint MaxDist = UINT_MAX;
	for (uint i = 0; i < N2; ++i)
		{
		const PDBChain &Chain = *Chains2[i];
		uint PosX = PosXs[i];
		uint PosA = PosAs[i];
		uint PosB = PosBs[i];
		uint PosC = PosCs[i];
		if (BeforeA)
			{
			asserta(PosA > PosX);
			uint Dist = PosA - PosX;
			if (i == 0)
				{
				MinDist = Dist;
				MaxDist = Dist;
				}
			else
				{
				MinDist = min(MinDist, Dist);
				MaxDist = min(MaxDist, Dist);
				}
			}
		}
	}

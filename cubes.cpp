#include "myutils.h"
#include "pdbchain.h"
#include "fakechain.h"

//static void GetCube(double x, double y, double z,
//					double Size,
//					vector<int> &Cube)
//	{
//	Cube.resize(3);
//	Cube[0] = int(round(x/Size));
//	Cube[1] = int(round(y/Size));
//	Cube[2] = int(round(z/Size));
//	}
//
///***
//3 axes with possible delta -1, 0, +1
//total possible deltas = 3^3 = 27
//***/
//bool GetCubeDeltas(const PDBChain &Chain, double Size, 
//			  vector<uint> &Deltas)
//	{
//	Deltas.clear();
//	const uint L = Chain.GetSeqLength();
//	vector<int> PrevCube;
//	vector<int> Cube;
//	for (uint i = 0; i < L; ++i)
//		{
//		double x = Chain.m_Xs[i];
//		double y = Chain.m_Ys[i];
//		double z = Chain.m_Zs[i];
//		if (i == 0)
//			{
//			GetCube(x, y, x, Size, PrevCube);
//			continue;
//			}
//		GetCube(x, y, x, Size, Cube);
//
//	// Represent -1,0,+1 as 0,1,2
//		int Diff0 = Cube[0] - PrevCube[0] + 1;
//		int Diff1 = Cube[1] - PrevCube[1] + 1;
//		int Diff2 = Cube[2] - PrevCube[2] + 1;
//		PrevCube = Cube;
//		if (Diff0 == 1 && Diff1 == 1 && Diff2 == 1)
//			continue;
//		if (Diff0 < 0 || Diff0 > 2) return false;
//		if (Diff1 < 0 || Diff1 > 2) return false;
//		if (Diff2 < 0 || Diff2 > 2) return false;
//		uint Delta = uint(Diff0) + 3*uint(Diff1) + 9*uint(Diff2);
//		assert(Delta < 28);
//		if (!Deltas.empty() && Delta == Deltas.back())
//			continue;
//		Log("%2d %2d %2d\n", int(Diff0)-1, int(Diff1)-1, int(Diff2)-1);
//		Deltas.push_back(Delta);
//		}
//	return true;
//	}
//
//static void TestDeltas()
//	{
//	const string &InputFN = g_Arg1;
//
//	vector<PDBChain *> Chains;
//	ReadChains(InputFN, Chains);
//	const double Size = 5;
//
//	const uint N = SIZE(Chains);
//	vector<uint> Deltas;
//	uint NotOk = 0;
//	for (uint i = 0; i < N; ++i)
//		{
//		ProgressStep(i, N, "Processing");
//		const PDBChain &Chain = *Chains[i];
//		const uint L = Chain.GetSeqLength();
//		const string &Label = Chain.m_Label;
//		bool Ok = GetCubeDeltas(Chain, Size, Deltas);
//		if (!Ok)
//			{
//			++NotOk;
//			continue;
//			}
//		const uint C = SIZE(Deltas);
//		Log("\n");
//		Log(">%s(%u, %u)\n", Label.c_str(), L, C);
//		for (uint c = 0; c < C; ++c)
//			{
//			uint Delta = Deltas[c];
//			uint Diff2 = Delta/9;
//			uint Diff1 = (Delta - Diff2*9)/3;
//			uint Diff0 = Delta%3;
//			asserta(Diff0 <= 2);
//			asserta(Diff1 <= 2);
//			asserta(Diff2 <= 2);
//			Log("%08x  %2d  %2d  %2d\n", Delta,
//				int(Diff0)-1,
//				int(Diff1)-1,
//				int(Diff2)-1);
//			}
//		}
//	ProgressLog("%u / %u not ok\n", NotOk, N);
//	}

void cmd_cubes()
	{
	const string &InputFN = g_Arg1;

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);
	const uint N = SIZE(Chains);

	FakeChain FC;

	vector<uint64> hs;
	hs.push_back(FC.GetHash(0, 0, 0));
	hs.push_back(FC.GetHash(1, 0, 0));
	hs.push_back(FC.GetHash(-1, 0, 0));
	hs.push_back(FC.GetHash(0, 0, 1));
	hs.push_back(FC.GetHash(0, 0, -1));
	hs.push_back(FC.GetHash(1, 2, 3));
	hs.push_back(FC.GetHash(2, 2, 3));

	bool Ok = FC.TryAppendFrag(*Chains[0]);
	asserta(Ok);

	for (uint k = 0; k < 3; ++k)
		{
		for (uint i = 1; i < N; ++i)
			{
			Ok = FC.TryAppendFrag(*Chains[i]);
			if (Ok)
				{
				Progress("Success %u\n", k);
				break;
				}
			}
		}
	if (!Ok)
		Progress("Fail\n");

	FC.m_Chain.ToPDB("cubes_out.pdb");
	}
#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"

static uint MinLen = 16;
static uint MaxLen = 32;
static uint RandTick = 64;

static bool IsCut(const string &SS, uint Pos)
	{
	if (randu32()%RandTick == 0)
		return true;
	const uint L = SIZE(SS);
	if (Pos < 1 || Pos + MinLen > L)
		return false;
	char c_1 = SS[Pos-1];
	char c0 = SS[Pos];
	bool PossibleCut = (c_1 != '~' && c0 == '~');
	if (!PossibleCut)
		return false;
	uint n = 0;
	for (uint i = Pos+1; i + 1 < L; ++i)
		{
		if (SS[i] == '~')
			++n;
		else
			break;
		}
	return n >= 5;
	}

static void GetCuts(const string &SS, vector<uint> &Cuts)
	{
	string scut;
	for (uint Pos = 0; Pos < SIZE(SS); ++Pos)
		{
		if (IsCut(SS, Pos))
			{
			scut += '|';
			Cuts.push_back(Pos);
			}
		else
			scut += '_';
		}
	}

static bool MakeFrag(PDBChain &Q, uint CutPos)
	{
	uint L = Q.GetSeqLength();
	bool Fwd = true;
	if (opt_randrev)
		Fwd = (randu32()%2 == 0);
	uint n = MinLen + randu32()%(MaxLen - MinLen);
	int Lo = -1;
	int Hi = -1;
	if (Fwd)
		{
		Lo = int(CutPos);
		Hi = Lo + int(n);
		}
	else
		{
		Lo = int(CutPos) - n;
		Hi = Lo + n;
		}
	if (Lo < 0 || Hi >= int(L))
		return false;

	PDBChain Frag;
	Q.GetRange(uint(Lo), uint(Hi), Frag);

	if (opt_flip)
		{
		PDBChain FragR;
		Frag.GetReverse(FragR);
		Ps(FragR.m_Label, "%s|%d-%d(%d%c)", Q.m_Label.c_str(), Lo, Hi, Hi-Lo+1, pom(Fwd));
		if (opt_shuffle)
			random_shuffle(FragR.m_Seq.begin(), FragR.m_Seq.end());
		FragR.ZeroOrigin();
		FragR.ToCal(g_fcal);
		}
	else
		{
		if (opt_shuffle)
			random_shuffle(Frag.m_Seq.begin(), Frag.m_Seq.end());
		Ps(Frag.m_Label, "%s[%d-%d]%d%c", Q.m_Label.c_str(), Lo, Hi, Hi-Lo+1, pom(Fwd));
		Frag.ZeroOrigin();
		Frag.ToCal(g_fcal);
		}
	return true;
	}

static void Frag(PDBChain &Q)
	{
	Q.SetSS();
	string SSS;
	vector<uint> Cuts;
	GetCuts(Q.m_SS, Cuts);
	const uint N = SIZE(Cuts);
	for (uint i = 0; i < N; ++i)
		{
		uint CutPos = Cuts[i];
		MakeFrag(Q, CutPos);
		}
	}

void cmd_fragment_library()
	{
	const string &QFN = opt_fragment_library;

	ChainReader CR;
	CR.Open(QFN);

	PDBChain Q;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;
		Frag(Q);
		}
	}

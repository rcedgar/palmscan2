#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"

static uint s_FragIdx;
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

static void GetSSum(const string &SS, uint Lo, uint Hi, string &SSum)
	{
	const uint L = SIZE(SS);
	asserta(Lo < Hi);
	asserta(Hi < L);
	uint ns = 0;
	uint nh = 0;
	uint nl = 0;
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		char c = SS[Pos];
		if (c == 'h')
			++nh;
		else if (c == 's')
			++ns;
		else if (c == '~' || c == 't')
			++nl;
		else
			Die("Bad ss '%c'", c);
		}
	Ps(SSum, "H%uS%uL%u", nh, ns, nl);
	}

static bool MakeFrag(const string &FN, PDBChain &Q, uint CutPos)
	{
	uint L = Q.GetSeqLength();
	uint n = MinLen + randu32()%(MaxLen - MinLen);
	int Lo = int(CutPos);
	int Hi = Lo + int(n);
	if (Lo < 0 || Hi >= int(L))
		return false;

	PDBChain Frag;
	Q.GetRange(uint(Lo), uint(Hi), Frag);

	if (opt_shuffle)
		random_shuffle(Frag.m_Seq.begin(), Frag.m_Seq.end());

	string SSum;
	GetSSum(Q.m_SS, uint(Lo), uint(Hi), SSum);
	Frag.ZeroOrigin();

	string &Label = Frag.m_Label;
	Ps(Label, "frag%u", s_FragIdx);
	Psa(Label, "|%s", FN.c_str());
	Psa(Label, "|%d", Lo);
	Psa(Label, "|%d", Hi);
	Psa(Label, "|%s", SSum.c_str());
	Frag.ToCal(g_fcal);

	s_FragIdx++;
	return true;
	}

static void Chop(const string &FN, PDBChain &Q)
	{
	Q.SetSS();
	vector<uint> Cuts;
	GetCuts(Q.m_SS, Cuts);
	const uint N = SIZE(Cuts);
	for (uint i = 0; i < N; ++i)
		{
		uint CutPos = Cuts[i];
		MakeFrag(FN, Q, CutPos);
		}
	}

void cmd_fragment_library()
	{
	string PDBDir = string(opt_fragment_library);
	if (!EndsWith(PDBDir, "/"))
		PDBDir += "/";

	vector<string> FNs;
	mylistdir(PDBDir, FNs);

	const uint N = SIZE(FNs);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Fragging");
		const string &FN = FNs[i];
		if (FN == "." || FN == "..")
			continue;
		string PathName = PDBDir + FN;

		vector<string> Lines;
		ReadLinesFromFile(PathName, Lines);

		string Label;
		Label = string(BaseName(FN.c_str()));

		PDBChain Q;
		Q.FromPDBLines(Label, Lines);

		Chop(FN, Q);
		}
	}

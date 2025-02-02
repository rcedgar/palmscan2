#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include "chainreader.h"

static uint s_FragIdx;
static uint MinLen = 6;
static uint MaxLen = 32;
static uint RandTick = 64;
static uint Frag_NH = 0;
static uint Frag_NS = 0;
static uint Frag_NL = 0;
static uint LoopPct = 5;

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

static void GetSSns(const string &SS, uint Lo, uint Hi,
					uint &nh, uint &ns, uint &nl)
	{
	const uint L = SIZE(SS);
	asserta(Lo < Hi);
	asserta(Hi < L);
	ns = 0;
	nh = 0;
	nl = 0;
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
	}

static bool MakeFrag(const string &FN, const PDBChain &Q, uint CutPos)
	{
	uint L = Q.GetSeqLength();
	uint n = MinLen + randu32()%(MaxLen - MinLen);
	int Lo = int(CutPos);
	int Hi = Lo + int(n) - 1;
	if (Lo < 0 || Hi >= int(L))
		return false;

	PDBChain Frag;
	Q.GetRange(uint(Lo), uint(Hi), Frag);

	string SSum;
	uint nh, ns, nl;
	GetSSns(Q.m_SS, Lo, Hi, nh, ns, nl);
	Ps(SSum, "H%uS%uL%u", nh, ns, nl);

	if (opt_strands)
		{
		if (nh > 0)
			return false;
		if (nl > ns)
		if (randu32()%100 > LoopPct)
			return false;
		}

	if (nl > nh && nl > ns)
		{
		if (randu32()%100 > LoopPct)
			return false;
		}

	if (opt_shuffle)
		random_shuffle(Frag.m_Seq.begin(), Frag.m_Seq.end());

	Frag_NH += nh;
	Frag_NS += ns;
	Frag_NL += nl;

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

static void Chop(const string &FN, const PDBChain &Q)
	{
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

	uint NH = 0;
	uint NS = 0;
	uint NL = 0;
	vector<string> FNs;
	mylistdir(PDBDir, FNs);

	const uint FileCount = SIZE(FNs);
	for (uint i = 0; i < FileCount; ++i)
		{
		ProgressStep(i, FileCount, "Fragging");
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
		Q.SetSS();
		uint L = Q.GetSeqLength();
		if (L < 16)
			continue;

		uint nh, ns, nl;
		GetSSns(Q.m_SS, 0, L-1, nh, ns, nl);
		NH += nh;
		NS += ns;
		NL += nl;

		Chop(FN, Q);
		}

	uint N = NS + NH + NL;
	ProgressLog("NS %u (%.1f%%), NH %u (%.1f%%), NL %u (%.1f%%)\n",
				NS, GetPct(NS, N),
				NH, GetPct(NH, N),
				NL, GetPct(NL, N));

	uint Frag_N = Frag_NS + Frag_NH + Frag_NL;
	ProgressLog("Frag_NS %u (%.1f%%), Frag_NH %u (%.1f%%), Frag_NL %u (%.1f%%)\n",
				Frag_NS, GetPct(Frag_NS, Frag_N),
				Frag_NH, GetPct(Frag_NH, Frag_N),
				Frag_NL, GetPct(Frag_NL, Frag_N));
	}

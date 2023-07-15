#include "myutils.h"
#include "dssp.h"

/***
    .-- sequential resnumber, including chain breaks as extra residues
    |    .-- original PDB resname, not nec. sequential, may contain letters
    |    | .-- one-letter chain ID, if any
    |    | | .-- amino acid sequence in one letter code
    |    | | |  .-- secondary structure summary based on columns 19-38
    |    | | |  | xxxxxxxxxxxxxxxxxxxx recommend columns for secstruc details
    |    | | |  | .-- 3-turns/helix
    |    | | |  | |.-- 4-turns/helix
    |    | | |  | ||.-- 5-turns/helix
    |    | | |  | |||.-- geometrical bend
    |    | | |  | ||||.-- chirality
    |    | | |  | |||||.-- beta bridge label
    |    | | |  | ||||||.-- beta bridge label
    |    | | |  | |||||||   .-- beta bridge partner resnum
    |    | | |  | |||||||   |   .-- beta bridge partner resnum
    |    | | |  | |||||||   |   |.-- beta sheet label
    |    | | |  | |||||||   |   ||   .-- solvent accessibility
    |    | | |  | |||||||   |   ||   |
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC
    |    | | |  | |||||||   |   ||   |
   35   47 A I  E     +     0   0    2

 DSSP Code       mmCIF Code      Description
     H         HELX\_RH\_AL\_P   Alphahelix
     B              STRN         Betabridge
     E              STRN           Strand
     G         HELX\_RH\_3T\_P    Helix\_3
     I         HELX\_RH\_PI\_P    Helix\_5
     P         HELX\_LH\_PP\_P   Helix\_PPII
     T          TURN\_TY1\_P        Turn
     S              BEND            Bend
 ' ' (space)        OTHER           Loop
***/

static int GetInt(const string &Line, uint ColFrom, uint ColTo)
	{
	string s;
	for (uint Col = ColFrom; Col <= ColTo; ++Col)
		{
		char c = Line[Col];
		if (!isspace(c))
			s +=c ;
		}
	int i = StrToInt(s);
	return i;
	}

static void GetStr(const string &Line, uint ColFrom, uint ColTo, string &s)
	{
	s.clear();
	for (uint Col = ColFrom; Col <= ColTo; ++Col)
		{
		char c = Line[Col];
		if (!isspace(c))
			s +=c ;
		}
	}

static uint GetUint(const string &Line, uint ColFrom, uint ColTo)
	{
	string s;
	for (uint Col = ColFrom; Col <= ColTo; ++Col)
		{
		char c = Line[Col];
		if (!isspace(c))
			s +=c ;
		}
	if (s.empty())
		return UINT_MAX;
	uint i = StrToUint(s);
	return i;
	}

void DSSP::FromLines(const string &Label, const vector<string> &Lines)
	{
	m_Label = Label;

/***
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA

          1         2         3         4         5         6
0123456789012345678901234567890123456789012345678901234567890
  662  662 A M  S  < S-     0   0   35     -4,-2.6  -526,-0.1  -313,-0.2  -525,-0.1  -0.155  92.2 -84.0 -69.0 154.4    4.7  -27.2   24.3
    6    6 A A  E     -A  379   0A  38    373,-0.2   373,-0.2    -2,-0.1   375,-0.1  -0.994   5.8-164.2-131.7 140.5  -10.4   -3.3   47.9
    1    1 A P              0   0  123      0, 0.0   237,-0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0  29.7  -20.8   -1.5   36.7
    2    2 A R        -     0   0  233      2,-0.0     2,-0.3     1,-0.0     0, 0.0   0.874 360.0-139.4  85.1  99.5  -18.9   -3.3   39.4
***/
	Clear();
	const uint N = SIZE(Lines);
	for (uint i = 0; i < N; ++i)
		{
		const string &Line = Lines[i];
		char ChainId = GetChainIdFromLine(Line);
		if (i == 0)
			m_ChainId = ChainId;
		else
			asserta(ChainId == m_ChainId);

		uint SequentialResNr = GetUint(Line, 0, 4);
		string OriginalResNr;
		GetStr(Line, 0, 4, OriginalResNr);
		char SS = Line[16];
		if (SS == ' ')
			SS = 'L';
		uint SolventAcc = GetUint(Line, 34, 37);
		uint Beta1 = GetUint(Line, 25, 28);
		uint Beta2 = GetUint(Line, 29, 32);
		if (Beta1 == 0)
			Beta1 = UINT_MAX;
		if (Beta2 == 0)
			Beta2 = UINT_MAX;

		m_SequentialResNrs.push_back(SequentialResNr);
		m_OriginalResNrs.push_back(OriginalResNr);
		m_Seq.push_back(Line[13]);
		m_SS.push_back(SS);
		m_SolventAccs.push_back(SolventAcc);
		m_BetaPartnerResNr1s.push_back(Beta1);
		m_BetaPartnerResNr2s.push_back(Beta2);
 		}
	}

void DSSP::PrintSeq(FILE *f) const
	{
	if (m_Label != "")
		fprintf(f, ">%s\n", m_Label.c_str());

	const uint L = GetSeqLength();
	asserta(SIZE(m_SS) == L);
	const unsigned ROWLEN = 80;
	unsigned BlockCount = (L + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned From = BlockIndex*ROWLEN;
		unsigned To = From + ROWLEN;
		if (To >= L)
			To = L;

		fputc('\n', f);

		for (unsigned Pos = From; Pos < To; ++Pos)
			{
			if (Pos%10 == 0)
				{
				string s;
				Ps(s, "%u", Pos);
				uint n = SIZE(s);
				fputs(s.c_str(), f);
				for (uint i = n; i < 10; ++i)
					fputc(' ', f);
				}
			}
		fputc('\n', f);

		for (unsigned Pos = From; Pos < To; ++Pos)
			{
			if (Pos%10 == 0)
				fputc(' ', f);
			else
				fprintf(f, "%u", Pos%10);
			}
		fputc('\n', f);

		for (unsigned Pos = From; Pos < To; ++Pos)
			fputc(tolower(m_SS[Pos]), f);
		fputc('\n', f);

		for (unsigned Pos = From; Pos < To; ++Pos)
			fputc(m_Seq[Pos], f);
		fputc('\n', f);
		}
	}

void DSSP::LogMe() const
	{
	const uint L = GetSeqLength();
	asserta(SIZE(m_SequentialResNrs) == L);
	asserta(SIZE(m_OriginalResNrs) == L);
	asserta(SIZE(m_SS) == L);
	asserta(SIZE(m_SolventAccs) == L);
	asserta(SIZE(m_BetaPartnerResNr1s) == L);
	asserta(SIZE(m_BetaPartnerResNr2s) == L);

	for (uint i = 0; i < L; ++i)
		{
		uint Beta1 = m_BetaPartnerResNr1s[i];
		uint Beta2 = m_BetaPartnerResNr2s[i];

		Log("%c", m_Seq[i]);
		Log("  %c", m_SS[i]);
		Log("  %4u ", m_SequentialResNrs[i]);
		Log("  %4.4s", m_OriginalResNrs[i]);
		if (Beta1 == UINT_MAX)
			Log("      ");
		else
			Log("  %4u", Beta1);
		if (Beta2 == UINT_MAX)
			Log("      ");
		else
			Log("  %4u", Beta2);
		Log("\n");
		}
	}

char DSSP::GetChainIdFromLine(const string &Line)
	{
	asserta(SIZE(Line) > 11);
	return Line[11];
	}

void DSSP::GetBetaPairs(vector<uint> &Start1s, vector<uint> &End1s,
  vector<uint> &Start2s, vector<uint> &End2s) const
	{
	const uint L = GetSeqLength();
	uint Start1 = UINT_MAX;
	uint Start2 = UINT_MAX;
	uint End1 = UINT_MAX;
	uint End2 = UINT_MAX;
	uint RunLength = 0;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		uint Beta1 = m_BetaPartnerResNr1s[Pos];
		if (Beta1 == UINT_MAX)
			{
			if (RunLength > 1)
				{
				Start1s.push_back(Start1);
				Start2s.push_back(Start2);
				End1s.push_back(End1);
				End2s.push_back(End2);
				}
			RunLength = 0;
			continue;
			}
		asserta(Beta1 >= 1 && Beta1 <= L);
		if (RunLength == 0)
			{
			Start1 = Pos;
			Start2 = Beta1;
			}
		else
			{
			End1 = Pos;
			End2 = Beta1;
			}
		++RunLength;
		}
	}

void DSSP::GetSSEls(vector<SSEl *> &Els) const
	{
	Els.clear();
	}

static uint GetSpan(uint i, uint j)
	{
	uint L = max(i, j) - min(i, j) + 1;
	return L;
	}

void DSSP::LogBetaPairs(const vector<uint> &Start1s, const vector<uint> &End1s,
  const vector<uint> &Start2s, const vector<uint> &End2s) const
	{
	const uint N = SIZE(Start1s);
	asserta(SIZE(End1s) == N);
	asserta(SIZE(Start2s) == N);
	asserta(SIZE(End2s) == N);

	Log("\n");
	Log("%u beta pairs\n", N);
	for (uint i = 0; i < N; ++i)
		{
		uint Start1 = Start1s[i];
		uint End1 = End1s[i];
		uint Start2 = Start2s[i];
		uint End2 = End2s[i];
		uint L1 = GetSpan(Start1, End1);
		uint L2 = GetSpan(Start2, End2);
		Log("(%4u  - %4u %u)  (%4u   - %4u %u)\n",
		  Start1, End1, L1, Start2, End2, L2);
		}
	}

void cmd_dssp()
	{
	string FN = opt_dssp;
	vector<string> Lines;
	ReadLinesFromFile(FN, Lines);
	DSSP D;
	D.FromLines(FN, Lines);
	//D.LogMe();
	D.PrintSeq(g_fLog);

	vector<uint> Start1s;
	vector<uint> Start2s;
	vector<uint> End1s;
	vector<uint> End2s;
	D.GetBetaPairs(Start1s, End1s, Start2s, End2s);
	D.LogBetaPairs(Start1s, End1s, Start2s, End2s);
	}

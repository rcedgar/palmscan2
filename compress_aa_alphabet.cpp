#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "sort.h"
#include <set>

/***
Input
=====
Pair-wise alignments of scop40 domains, created by
  palmscan2 -tm_scop d:/int/scop40/out/domains_scop.cal -output scop40.fa2 -tsv scop40.tsv -log tm.log -threads 10

K=8 0.122262        AST,C,DN,EHKQR,FWY,G,ILMV,P
***/

static void GetAlphaStr(vector<vector<uint> > &GroupToLetters, string &AlphaStr)
	{
	AlphaStr.clear();
	vector<string> Strs;
	for (uint i = 0; i < SIZE(GroupToLetters); ++i)
		{
		const vector<uint> &Letters = GroupToLetters[i];
		vector<char> Chars;
		for (uint j = 0; j < SIZE(Letters); ++j)
			Chars.push_back(g_LetterToCharAmino[Letters[j]]);
		sort(Chars.begin(), Chars.end());
		string Str;
		for (uint j = 0; j < SIZE(Letters); ++j)
			Str += Chars[j];
		Strs.push_back(Str);
		}
	sort(Strs.begin(), Strs.end());
	for (uint i = 0; i < SIZE(Strs); ++i)
		{
		if (i > 0)
			AlphaStr += ",";
		AlphaStr += Strs[i];
		}
	}

static void SelectRandomPair(const vector<uint> &Letter1s,
  const vector<uint> &Letter2s, uint &Letter1, uint &Letter2)
	{
	const uint PairCount = SIZE(Letter1s);
	uint r = randu32()%PairCount;
	Letter1 = Letter1s[r];
	Letter2 = Letter2s[r];
	}

void BuildAlphabet(
  const vector<uint> &Letter1s,
  const vector<uint> &Letter2s,
  uint N,
  vector<vector<uint> > &GroupToLetters,
  vector<uint> &LetterToGroup)
	{
	GroupToLetters.clear();
	GroupToLetters.resize(N);
	set<uint> Assigned;
	LetterToGroup.clear();
	LetterToGroup.resize(20, UINT_MAX);
	for (uint i = 0; i < N; ++i)
		{
		uint Letter = UINT_MAX;
		for (uint j = 0; j < 100; ++j)
			{
			uint r = randu32()%20;
			if (Assigned.find(r) == Assigned.end())
				{
				Letter = r;
				break;
				}
			}
		asserta(Letter != UINT_MAX);
		Assigned.insert(Letter);
		GroupToLetters[i].push_back(Letter);
		LetterToGroup[Letter] = i;
		}
	const uint MAXITERS = 1000;
	for (uint Iter = 0; Iter < 1000; ++Iter)
		{
		if (SIZE(Assigned) == 20)
			return;
		uint Letter1, Letter2;
		SelectRandomPair(Letter1s, Letter2s, Letter1, Letter2);
		uint Group1 = LetterToGroup[Letter1];
		uint Group2 = LetterToGroup[Letter2];

		if (Group1 != UINT_MAX && Group2 != UINT_MAX)
			{
		// Both letters already in groups, merge these groups
			continue;
			}
		else if (Group1 != UINT_MAX && Group2 == UINT_MAX)
			{
		// Letter1 in group, Letter2 not in group, assign Letter2 to Group1
			LetterToGroup[Letter2] = Group1;
			GroupToLetters[Group1].push_back(Letter2);
			Assigned.insert(Letter2);
			}
		else if (Group1 == UINT_MAX && Group2 != UINT_MAX)
			{
		// Letter2 in group, Letter1 not in group, assign Letter1 to Group2
			LetterToGroup[Letter1] = Group2;
			GroupToLetters[Group2].push_back(Letter1);
			Assigned.insert(Letter1);
			}
		else if (Group1 == UINT_MAX && Group2 == UINT_MAX)
			{
		// Neither in group
			continue;
			}
		else
			asserta(false);
		}
	for (uint i = 0; i < 20; ++i)
		asserta(LetterToGroup[i] != UINT_MAX);
	asserta(SIZE(GroupToLetters) == N);
	}

static double GetExpectedScore(const vector<uint> &Counts, 
  const vector<vector<uint> > &CountMx)
	{
	const uint K = SIZE(Counts);
	asserta(SIZE(CountMx) == K);
	for (uint i = 0; i < K; ++i)
		asserta(SIZE(CountMx) == K);
	uint Total = 0;
	for (uint i = 0; i < K; ++i)
		Total += Counts[i];

	uint Total2 = 0;
	for (uint i = 0; i < K; ++i)
		for (uint j = 0; j < K; ++j)
			Total2 += CountMx[i][j];

	asserta(Total2 == Total);
	vector<double> Freqs;
	double SumFreq = 0;
	for (uint i = 0; i < K; ++i)
		{
		double Freq = double(Counts[i])/Total;
		Freqs.push_back(Freq);
		SumFreq += Freq;
		}
	asserta(feq(SumFreq, 1.0));

	double SumFreq2 = 0;
	double ES = 0;
	for (uint i = 0; i < K; ++i)
		{
		for (uint j = 0; j < K; ++j)
			{
			uint n = CountMx[i][j];
			double ObsFreq = double(n)/double(Total);
			SumFreq2 += ObsFreq;
			double ExpFreq = double(Freqs[i]*Freqs[j]);
			double Ratio = ObsFreq/ExpFreq;
			double Score = log(Ratio);
			ES += ObsFreq*Score;
			}
		}
	asserta(feq(SumFreq, 1.0));
	return ES;
	}

static void PeturbAlphabet(uint K, const vector<uint> &InLetterToGroup,
  vector<vector<uint> > &GroupToLetters, vector<uint> &LetterToGroup)
	{
	asserta(SIZE(InLetterToGroup) == 20);
	GroupToLetters.clear();
	GroupToLetters.resize(K);
	LetterToGroup = InLetterToGroup;

	uint r1 = UINT_MAX;
	uint r2 = UINT_MAX;
	do
		{
		r1 = randu32()%20;
		r2 = randu32()%20;
		if (r2 == r1)
			r2 = (r2 + 1)%20;
		}
	while(LetterToGroup[r1] == LetterToGroup[r2]);
	if (randu32()%2 == 0)
		swap(LetterToGroup[r1], LetterToGroup[r2]);
	else
		LetterToGroup[r1] = LetterToGroup[r2];
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		uint k = LetterToGroup[Letter];
		asserta(k < K);
		GroupToLetters[k].push_back(Letter);
		}
	}
static uint AlphabetFromAlphaStr(const string &AlphaStr,
  vector<vector<uint> > &GroupToLetters, vector<uint> &LetterToGroup)
	{
	GroupToLetters.clear();
	LetterToGroup.clear();
	LetterToGroup.resize(20, UINT_MAX);
	vector<string> Strs;
	Split(AlphaStr, Strs, ',');
	const uint K = SIZE(Strs);
	GroupToLetters.resize(K);
	for (uint k = 0; k < K; ++k)
		{
		const string &Str = Strs[k];
		for (uint i = 0; i < SIZE(Str); ++i)
			{
			char c = Str[i];
			uint Letter = g_CharToLetterAmino[c];
			asserta(Letter < 20);
			LetterToGroup[Letter] = k;
			GroupToLetters[k].push_back(Letter);
			}
		}
	for (uint i = 0; i < 20; ++i)
		asserta(LetterToGroup[i] != UINT_MAX);
	return K;
	}

void cmd_compress_aa_alphabet()
	{
	const uint K = opt_k;

	SeqDB Input;
	Input.FromFasta(opt_compress_aa_alphabet, true);

	asserta(optset_mintm && optset_maxtm);
	double MinTM = opt_mintm;
	double MaxTM = opt_maxtm;

	uint TriangleSize = (20*19)/2;

	vector<uint> OffDiagCountVec(TriangleSize, 0);
	vector<vector<uint> > TriangleIdx(20);
	vector<uint> IdxToi;
	vector<uint> IdxToj;
	for (uint i = 0; i < 20; ++i)
		TriangleIdx[i].resize(20);
	uint Idx = 0;
	for (uint i = 0; i < 20; ++i)
		{
		TriangleIdx[i][i] = UINT_MAX;
		for (uint j = 0; j < i; ++j)
			{
			TriangleIdx[i][j] = Idx;
			TriangleIdx[j][i] = Idx;
			IdxToi.push_back(i);
			IdxToj.push_back(j);
			++Idx;
			}
		}

	const uint SeqCount = Input.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;
	uint OffDiagCount = 0;
	vector<uint> CountVec(20);
	vector<uint> Letter1s;
	vector<uint> Letter2s;
	vector<vector<uint> > CountMx(20);
	for (uint i = 0; i < 20; ++i)
		CountMx[i].resize(20);
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Processing");
		const string &QRow = Input.GetSeq(2*PairIndex);
		const string &RRow = Input.GetSeq(2*PairIndex+1);
		const string &QLabel = Input.GetLabel(2*PairIndex);
		const string &RLabel = Input.GetLabel(2*PairIndex+1);
		const uint ColCount = SIZE(QRow);
		asserta(SIZE(QRow) == SIZE(RRow));
		vector<string> Fields;
		Split(QLabel, Fields, '/');
		asserta(SIZE(Fields) == 4);
		const string &Dom = Fields[0];
		const string &Fam = Fields[1];
		const string sTM = Fields[2];
		const string sPctId = Fields[3];
		double TM = StrToFloat(sTM);
		if (TM < MinTM || TM > MaxTM)
			continue;

		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char q = QRow[Col];
			char r = RRow[Col];
			if (isgap(q) || isgap(r))
				continue;
			uint iq = g_CharToLetterAmino[q];
			uint ir = g_CharToLetterAmino[r];
			if (iq >= 20 || ir >= 20)
				continue;
			LetterPairCount += 2;
			CountVec[iq] += 1;
			CountVec[ir] += 1;
			CountMx[iq][ir] += 1;
			CountMx[ir][iq] += 1;
			if (iq != ir)
				{
				++OffDiagCount;
				Letter1s.push_back(iq);
				Letter2s.push_back(ir);
				uint Idx = TriangleIdx[iq][ir];
				asserta(Idx < SIZE(OffDiagCountVec));
				OffDiagCountVec[Idx] += 1;
				}
			}
		}

	double SumFreq = 0;
	vector<double> Freqs;
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		double Freq = double(CountVec[i])/(LetterPairCount);
		SumFreq += Freq;
		Log("%c %.4f\n", c, Freq);
		Freqs.push_back(Freq);
		}
	asserta(feq(SumFreq, 1.0));

	double SumFreq2 = 0;
	vector<double> Freqs2;
	for (uint i = 0; i < 20; ++i)
		{
		for (uint j = 0; j < i; ++j)
			{
			uint Idx = TriangleIdx[i][j];
			uint n = OffDiagCountVec[Idx];
			double Freq2 = double(n)/OffDiagCount;
			SumFreq2 += Freq2;
			Freqs2.push_back(Freq2);
			}
		}
	asserta(feq(SumFreq2, 1.0));
	vector<uint> Order(TriangleSize);
	QuickSortOrderDesc(Freqs2.data(), TriangleSize, Order.data());
	for (uint k = 0; k < TriangleSize; ++k)
		{
		Idx = Order[k];
		uint i = IdxToi[Idx];
		uint j = IdxToj[Idx];
		uint n = OffDiagCountVec[Idx];
		double Freq2 = double(n)/OffDiagCount;
		Log("%c %c %8.6f\n",
			g_LetterToCharAmino[i],
			g_LetterToCharAmino[j],
			Freq2);
		}

	string BestAlphaStr;
	double BestES = 0;
	uint BestA = 0;
	vector<vector<uint> > BestGroupToLetters;
	vector<uint> BestLetterToGroup;
	const int MAXALPHAS = 1000*1000*5;
	const uint MAXSTALLED = 250000;
	const char *Action = "Rand";
	uint PeturbTriesRemaining = 0;
	uint Stalled = 0;
	string InitAlphaStr = opt_alphastr;
	for (int AlphaIdx = 0; AlphaIdx < MAXALPHAS; ++AlphaIdx)
		{
		ProgressStep(AlphaIdx, MAXALPHAS, "%s K=%u %.4f %s (%.1f%%)",
		  Action, K, BestES, BestAlphaStr.c_str(), GetPct(BestA, MAXALPHAS));

		vector<vector<uint> > GroupToLetters;
		vector<uint> LetterToGroup;
		if (InitAlphaStr != "")
			{
			uint K2 = AlphabetFromAlphaStr(InitAlphaStr,
			  GroupToLetters, LetterToGroup);
			InitAlphaStr.clear();
			asserta(K2 == K);
			}
		else if (PeturbTriesRemaining == 0)
			{
			Action = "Rand";
			BuildAlphabet(Letter1s, Letter2s, K, GroupToLetters, LetterToGroup);
			}
		else
			{
			Action = "Peturb";
			PeturbAlphabet(K, BestLetterToGroup, GroupToLetters, LetterToGroup);
			asserta(LetterToGroup != BestLetterToGroup);
			--PeturbTriesRemaining;
			}
		vector<vector<uint > > GroupCountMx(K);
		for (uint i = 0; i < K; ++i)
			GroupCountMx[i].resize(K);
		vector<uint> GroupCounts(K);
		for (uint i = 0; i < 20; ++i)
			{
			uint Groupi = LetterToGroup[i];
			for (uint j = 0; j < 20; ++j)
				{
				uint Groupj = LetterToGroup[j];
				uint n = CountMx[i][j];
				GroupCounts[Groupi] += n;
				GroupCounts[Groupj] += n;
				GroupCountMx[Groupi][Groupj] += 2*n;
				}
			}

		string AlphaStr;
		GetAlphaStr(GroupToLetters, AlphaStr);
		double ES = GetExpectedScore(GroupCounts, GroupCountMx);
		if (ES > BestES)
			{
			Stalled = 0;
			PeturbTriesRemaining = 256;
			BestA = AlphaIdx;
			BestES = ES;
			BestAlphaStr = AlphaStr;
			BestLetterToGroup = LetterToGroup;
			BestGroupToLetters = GroupToLetters;
			Log("%s K=%u %.8f %.8f %s (%u, %.3f%%)\n",
			  Action, K, ES, BestES, BestAlphaStr.c_str(), BestA, GetPct(BestA, MAXALPHAS));
			}
		else
			++Stalled;
		if (0) //(Stalled > MAXSTALLED)
			{
			ProgressStep(MAXALPHAS-1, MAXALPHAS, "STALLED K=%u %.4f %s (%.1f%%)",
			  K, BestES, BestAlphaStr.c_str(), GetPct(BestA, MAXALPHAS));
			break;
			}
		}

	FILE *fOut = CreateStdioFile(opt_output);
	if (fOut)
		fprintf(fOut, "%.6f\t%s\n", BestES, BestAlphaStr.c_str());
	CloseStdioFile(fOut);
	}

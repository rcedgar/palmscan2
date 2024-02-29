#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "sort.h"
#include <set>

static void GetKmers(const string &Seq, uint L, uint k, vector<uint> &Kmers)
	{
	Kmers.clear();
	for (uint Pos = 0; Pos + k < L; ++Pos)
		{
		const char *s = Seq.c_str() + Pos;
		uint Kmer = StrToWordAmino(s, k);
		if (Kmer != UINT_MAX)
			{
			Kmers.push_back(Kmer);
			//Log("[%5u]", Pos);
			//Log("  %*.*s", k, k, s);
			//Log("  %*.*s", k, k, WordToStrAmino(Kmer, k));
			//Log("\n");
			}
		}
	}

// K-mer stats on foldseek 3di sequences, correlating 
//  with SCOPclassifications
void cmd_threedix()
	{
	SeqDB Seqs;
	Seqs.FromFasta(opt_threedix);

	uint k = 5;
	if (optset_k)
		k = opt_k;

	const string &TsvFN= "d:/a/res/scop/out/legacy_domains.tsv";
	FILE *f = OpenStdioFile(TsvFN);
	string Line;
	vector<string> Fields;
/**
Format is class.fold.superfamily.family
a.1.1.1 d1dlwa_
a.1.1.1 d1uvya_
a.1.1.1 d1dlya_
...
**/
	set<string> FamSet;
	vector<string> Fams;
	map<string, uint> FamToIdx;
	vector<vector<string> > FamToDoms;
	map<string, string> DomToFam;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		const string Fam = Fields[0];
		const string Dom = Fields[1];
		if (FamSet.find(Fam) == FamSet.end())
			{
			uint Idx = SIZE(Fams);
			FamSet.insert(Fam);
			Fams.push_back(Fam);
			FamToIdx[Fam] = Idx;
			vector<string> Doms;
			Doms.push_back(Dom);
			FamToDoms.push_back(Doms);
			}
		DomToFam[Dom] = Fam;
		}

	const uint SeqCount = Seqs.GetSeqCount();
	uint NotFound = 0;
	uint Found = 0;
	vector<uint> SeqFamIdxs;
	uint DictSize = myipow(20, k);
	Progress("Dict size %s\n", IntToStr(DictSize));
	vector<vector<uint> > KmerToFamIdxs(DictSize);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Pass 1, %u not found", NotFound);
		const string &Seq = Seqs.GetSeq(SeqIndex);
		const string &Dom = Seqs.GetLabel(SeqIndex);
		uint L = Seqs.GetSeqLength(SeqIndex);
		map<string, string>::const_iterator iter = DomToFam.find(Dom);
		if (iter == DomToFam.end())
			{
			//Log("Not in scop >%s\n", Dom.c_str());
			SeqFamIdxs.push_back(UINT_MAX);
			++NotFound;
			continue;
			}
		const string &Fam = iter->second;
		uint FamIdx = FamToIdx[Fam];
		SeqFamIdxs.push_back(FamIdx);
		++Found;
		vector<uint> Kmers;
		GetKmers(Seq, L, k, Kmers);
		for (uint i = 0; i < SIZE(Kmers); ++i)
			{
			uint Kmer = Kmers[i];
			KmerToFamIdxs[Kmer].push_back(FamIdx);
			}
		}
	const uint FamCount = SIZE(FamSet);
	ProgressLog("%u fams, %u doms found, %u doms missing\n", FamCount, Found, NotFound);

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Pass 2");
		uint FamIdx = SeqFamIdxs[SeqIndex];
		const string &Seq = Seqs.GetSeq(SeqIndex);
		const string &Dom = Seqs.GetLabel(SeqIndex);
		uint L = Seqs.GetSeqLength(SeqIndex);

		vector<uint> Kmers;
		GetKmers(Seq, L, k, Kmers);

		map<uint, uint> FamIdxToCount;
		vector<uint> FamIdxs;

		for (uint i = 0; i < SIZE(Kmers); ++i)
			{
			uint Kmer = Kmers[i];
			const vector<uint> &KmerFamIdxs = KmerToFamIdxs[Kmer];
			const uint N = SIZE(KmerFamIdxs);
			asserta(N > 0);
			for (uint j = 0; j < N; ++j)
				{
				uint KmerFamIdx = KmerFamIdxs[j];
				asserta(KmerFamIdx < FamCount);
				if (FamIdxToCount.find(KmerFamIdx) == FamIdxToCount.end())
					{
					FamIdxs.push_back(KmerFamIdx);
					FamIdxToCount[KmerFamIdx] = 1;
					}
				else
					FamIdxToCount[KmerFamIdx] += 1;
				}
			}
		vector<uint> FamCounts;
		const uint F = SIZE(FamIdxs);
		for (uint j = 0; j < F; ++j)
			{
			uint KmerFamIdx  = FamIdxs[j];
			uint n = FamIdxToCount[KmerFamIdx];
			FamCounts.push_back(n);
			}
		vector<uint> Order(F);
		QuickSortOrderDesc(FamCounts.data(), F, Order.data());
		Log(">%s", Dom.c_str());
		for (uint j = 0; j < F; ++j)
			{
			uint jj = Order[j];
			const string &Fam = Fams[FamIdxs[jj]];
			Log(" %s(%u)", Fam.c_str(), FamCounts[jj]);
			}
		Log("\n");
		}
	}

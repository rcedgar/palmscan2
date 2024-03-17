#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "quarts.h"
#include "outputfiles.h"

/***
With CMAA:
Min=-190, LoQ=147, Med=204, HiQ=279, Max=1.29e+03, Avg=227, StdDev=126
Min=-1.61e+03, LoQ=-321, Med=-188, HiQ=-109, Max=166, Avg=-225, StdDev=155

No CMAA, 8 NUDX bins
Min=-201, LoQ=137, Med=193, HiQ=264, Max=1.15e+03, Avg=213, StdDev=119
Min=-1.36e+03, LoQ=-308, Med=-179, HiQ=-100, Max=164, Avg=-213, StdDev=149

Without CMAA size = 8x4x4x4 = 512
	case 0: return m_NUDX_Bins=8;
	case 1: return 4;
	case 2: return 4;
	case 3: return 4;

No CMAA, 4 NUDX bins, size = 256
Min=5.1, LoQ=178, Med=236, HiQ=316, Max=1.53e+03, Avg=263, StdDev=130
Min=-696, LoQ=-148, Med=-73.7, HiQ=-10.7, Max=363, Avg=-80.6, StdDev=95
***/

void GetScopDomFromLabel(const string &Label, string &Dom);

static void LogHist(const char *Name,
  const vector<double> &Scores, vector<uint> &BinToCount)
	{
	BinToCount.clear();
	BinToCount.resize(21);
	const uint N = SIZE(Scores);
	uint MaxCount = 0;
	for (uint i = 0; i < N; ++i)
		{
		double Score = Scores[i];
		int Bin = int(20.0*(Score + 3.0)/6.0);
		if (Bin < 0)
			Bin = 0;
		if (Bin > 20)
			Bin = 20;
		uint Count = BinToCount[Bin] + 1;
		BinToCount[Bin] = Count;
		MaxCount = max(Count, MaxCount);
		}
	const uint W = 40;
	Log("\n");
	Log("HISTO(%s)\n", Name);
	for (int Bin = 0; Bin <= 20; ++Bin)
		{
		uint n = BinToCount[Bin];
		double w = (n*W)/MaxCount;
		double x = 3.0*(Bin - 10)/10.0;
		Log("%+7.2f   %10u  ", x, n);
		for (uint i = 0; i < w; ++i)
			Log("*");
		Log("\n");
		}
	}

void cmd_testdss()
	{
	SeqDB Input;
	Input.FromFasta(opt_testdss, true);

	double MinTM = 0.6;
	double MaxTM = 0.8;
	if (optset_mintm)
		MinTM = opt_mintm;
	if (optset_maxtm)
		MaxTM = opt_maxtm;

	vector<PDBChain *> Chains;
	ReadChains(opt_train_cal, Chains);
	const uint ChainCount = SIZE(Chains);
	map<string, uint> DomToChainIndex;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const string &Label = Chains[ChainIndex]->m_Label;
		vector<string> Fields;
		string Dom;
		GetScopDomFromLabel(Label, Dom);
		DomToChainIndex[Dom] = ChainIndex;
		}

	const uint SeqCount = Input.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;

	DSS DSS_Q;
	DSS DSS_R;
	vector<double> AlnScores;
	vector<double> ShuffledAlnScores;
	double SumAlnScores = 0;
	double SumShuffledAlnScores = 0;
	double MeanAlnScore = 0;
	double MeanShuffledAlnScore = 0;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount,
		  "Mean aln score %.3g / %.3g",
		  MeanAlnScore, MeanShuffledAlnScore);
		const string &QLabel = Input.GetLabel(2*PairIndex);
		const string &RLabel = Input.GetLabel(2*PairIndex+1);
		vector<string> Fields;
		Split(QLabel, Fields, '/');
		asserta(SIZE(Fields) == 4);
		const string &QDom = Fields[0];
		const string &Fam = Fields[1];
		const string sTM = Fields[2];
		const string sPctId = Fields[3];
		double TM = StrToFloat(sTM);
		if (TM < MinTM || TM > MaxTM)
			continue;

		string RDom;
		GetScopDomFromLabel(RLabel, RDom);
		uint QChainIndex = DomToChainIndex[QDom];
		uint RChainIndex = DomToChainIndex[RDom];
		const PDBChain &QChain = *Chains[QChainIndex];
		const PDBChain &RChain = *Chains[RChainIndex];
		uint QL = QChain.GetSeqLength();
		uint RL = RChain.GetSeqLength();
		DSS_Q.Init(QChain);
		DSS_R.Init(RChain);
		const string &QRow = Input.GetSeq(2*PairIndex);
		const string &RRow = Input.GetSeq(2*PairIndex+1);
		const uint ColCount = SIZE(QRow);
		asserta(SIZE(QRow) == SIZE(RRow));

		uint QPos = 0;
		uint RPos = 0;
		double AlnScore = 0;
		double ShuffledAlnScore = 0;
		uint MatchCount = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char q = QRow[Col];
			char r = RRow[Col];
			if (!isgap(q) && !isgap(r))
				{
				uint ShuffledQPos = randu32()%QL;

				uint i0 = DSS_Q.GetFeature(0, QPos);
				uint i1 = DSS_Q.GetFeature(1, QPos);
				uint i2 = DSS_Q.GetFeature(2, QPos);
				uint i3 = DSS_Q.GetFeature(3, QPos);
				uint i4 = DSS_Q.GetFeature(4, QPos);

				uint is0 = DSS_Q.GetFeature(0, ShuffledQPos);
				uint is1 = DSS_Q.GetFeature(1, ShuffledQPos);
				uint is2 = DSS_Q.GetFeature(2, ShuffledQPos);
				uint is3 = DSS_Q.GetFeature(3, ShuffledQPos);
				uint is4 = DSS_Q.GetFeature(4, ShuffledQPos);

				uint j0 = DSS_R.GetFeature(0, RPos);
				uint j1 = DSS_R.GetFeature(1, RPos);
				uint j2 = DSS_R.GetFeature(2, RPos);
				uint j3 = DSS_R.GetFeature(3, RPos);
				uint j4 = DSS_R.GetFeature(4, RPos);

				//double Score = DSS::GetScore(
				//  i0, i1, i2, i3, i4,
				//  j0, j1, j2, j3, j4);

				//double ShuffledScore = DSS::GetScore(
				//  is0, is1, is2, is3, is4,
				//  j0, j1, j2, j3, j4);

				double Score = DSS::GetScore_NoCMAA(
				  i0, i1, i2, i3,
				  j0, j1, j2, j3);

				double ShuffledScore = DSS::GetScore_NoCMAA(
				  is0, is1, is2, is3,
				  j0, j1, j2, j3);

				AlnScore += Score;
				ShuffledAlnScore += ShuffledScore;
				++MatchCount;
				}
			if (!isgap(q))
				++QPos;
			if (!isgap(r))
				++RPos;
			}
		if (MatchCount == 0)
			MatchCount = 1;
		AlnScores.push_back(AlnScore/MatchCount);
		ShuffledAlnScores.push_back(ShuffledAlnScore/MatchCount);
		SumAlnScores += AlnScore;
		SumShuffledAlnScores += ShuffledAlnScore;
		MeanAlnScore = SumAlnScores/(PairIndex+1);
		MeanShuffledAlnScore = SumShuffledAlnScores/(PairIndex+1);
		}

	QuartsDouble QD;
	GetQuartsDouble(AlnScores, QD);
	QD.LogMe();
	GetQuartsDouble(ShuffledAlnScores, QD);
	QD.LogMe();

	vector<uint> BinToCount;
	vector<uint> ShuffledBinToCount;
	LogHist("Aligned", AlnScores, BinToCount);
	LogHist("Shuffled", ShuffledAlnScores, ShuffledBinToCount);
	const uint BinCount = SIZE(BinToCount);
	asserta(SIZE(ShuffledBinToCount) == BinCount);
	for (uint i = 0; i < BinCount; ++i)
		Log("%u\t%u\t%u\n", i, BinToCount[i], ShuffledBinToCount[i]);
	}

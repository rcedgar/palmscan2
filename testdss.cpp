#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "quarts.h"
#include "outputfiles.h"

void GetScopDomFromLabel(const string &Label, string &Dom);

void cmd_testdss()
	{
	SeqDB Input;
	Input.FromFasta(opt_testdss, true);

	const uint DSSFeatureIndex = 99;

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

	DSS QX;
	DSS RX;
	vector<double> AlnScores;
	vector<double> ShuffledAlnScores;
	double SumAlnScores = 0;
	double SumShuffledAlnScores = 0;
	double MeanAlnScore = 0;
	double MeanShuffledAlnScore = 0;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Mean aln score %.3g / %.3g",
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
		QX.Init(QChain);
		RX.Init(RChain);
		const string &QRow = Input.GetSeq(2*PairIndex);
		const string &RRow = Input.GetSeq(2*PairIndex+1);
		const uint ColCount = SIZE(QRow);
		asserta(SIZE(QRow) == SIZE(RRow));

		uint QPos = 0;
		uint RPos = 0;
		double AlnScore = 0;
		double ShuffledAlnScore = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char q = QRow[Col];
			char r = RRow[Col];
			if (!isgap(q) && !isgap(r))
				{
				uint i0 = QX.GetFeature(0, QPos);
				uint ShuffledQPos = randu32()%QL;

				uint i1 = QX.GetFeature(1, QPos);
				uint i2 = QX.GetFeature(2, QPos);
				uint i3 = QX.GetFeature(3, QPos);
				uint i4 = QX.GetFeature(4, QPos);

				uint is0 = QX.GetFeature(0, ShuffledQPos);
				uint is1 = QX.GetFeature(1, ShuffledQPos);
				uint is2 = QX.GetFeature(2, ShuffledQPos);
				uint is3 = QX.GetFeature(3, ShuffledQPos);
				uint is4 = QX.GetFeature(4, ShuffledQPos);

				uint j0 = RX.GetFeature(0, RPos);
				uint j1 = RX.GetFeature(1, RPos);
				uint j2 = RX.GetFeature(2, RPos);
				uint j3 = RX.GetFeature(3, RPos);
				uint j4 = RX.GetFeature(4, RPos);

				double Score = DSS::GetScore(
				  i0, i1, i2, i3, i4,
				  j0, j1, j2, j3, j4);

				double ShuffledScore = DSS::GetScore(
				  is0, is1, is2, is3, is4,
				  j0, j1, j2, j3, j4);

				AlnScore += Score;
				ShuffledAlnScore += ShuffledScore;
				}
			if (!isgap(q))
				++QPos;
			if (!isgap(r))
				++RPos;
			}
		AlnScores.push_back(AlnScore);
		ShuffledAlnScores.push_back(ShuffledAlnScore);
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
	}

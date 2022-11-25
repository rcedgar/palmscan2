#include "myutils.h"
#include "seqdb.h"
#include "abcxyz.h"
#include "mx.h"
#include "alpha.h"
#include "mpcluster.h"
#include "outputfiles.h"
#include "sort.h"

extern float **g_SubstMx;

void SetBLOSUM62();

char GetAnnotChar(uint Letter1, uint Letter2)
	{
	if (Letter1 == Letter2)
		return '|';
	float Score = g_SubstMx[Letter1][Letter2];
	if (Score >= 1)
		return '+';
	if (Score > 0)
		return '.';
	return ' ';
	}

void MPCluster::LogLogos(const vector<MotifProfile *> &MPs)
	{
	for (uint i = 0; i < SIZE(MPs); ++i)
		{
		string Logo;
		MPs[i]->GetLogo(Logo);
		Log("[%5u]  %s\n", i, Logo.c_str());
		}
	}

float MPCluster::GetScore(const MotifProfile &MP1,
  const MotifProfile &MP2) const
	{
#if DEBUG
	MP1.ValidateFreqs();
	MP2.ValidateFreqs();
#endif
	float Sum = 0;
	for (uint i = 0; i < MPL; ++i)
		{
		for (uint j = 0; j < 20; ++j)
			{
			float f1 = MP1.m_FreqVec[i][j];
			if (f1 < 1e-6)
				continue;
			byte c1 = g_LetterToCharAmino[j];
			for (uint k = 0; k < 20; ++k)
				{
				byte c2 = g_LetterToCharAmino[k];
				float f2 = MP2.m_FreqVec[i][k];
				float Score = f1*f2*g_SubstMx[c1][c2];
				Sum += Score;
				}
			}
		}
	float Score = Sum/MPL;
	return Score;
	}

void MPCluster::LogPair(const MotifProfile &MP1,
  const MotifProfile &MP2) const
	{
	string Logo1;
	string Logo2;
	MP1.GetLogo(Logo1);
	MP2.GetLogo(Logo2);

	float Score = GetScore(MP1, MP2);

	asserta(Logo1.size() == MPL);
	asserta(Logo2.size() == MPL);
	string Annot;
	for (uint i = 0; i < MPL; ++i)
		{
		char a = GetAnnotChar(Logo1[i], Logo2[i]);
		Annot.push_back(a);
		}
	Log("\n");
	Log("%s\n", Logo1.c_str());
	Log("%s\n", Annot.c_str());
	Log("%s\n", Logo2.c_str());
	Log("Score = %.3f\n", Score);
	}

void MotifProfile::GetLettersFromXxxSeq(const string &Seq,
  vector<uint> &Letters)
	{
	Letters.clear();
	asserta(SIZE(Seq) == AL+3+BL+3+CL);

	asserta(Seq[12] == 'x');
	asserta(Seq[13] == 'x');
	asserta(Seq[14] == 'x');

	asserta(Seq[14+14+1] == 'x');
	asserta(Seq[14+14+1] == 'x');
	asserta(Seq[14+14+2] == 'x');

	const uint PosVec[MPL] =
		{
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
		// 12, 13, 14,
		15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
		// 29, 30, 31
		32, 33, 34, 35, 36, 37, 38, 39
		};

	for (uint k = 0; k < MPL; ++k)
		{
		uint i = PosVec[k];
		byte c = (byte) Seq[i];
		assert(c != 0);
		uint Letter = g_CharToLetterAmino[c];
		Letters.push_back(Letter);
		}
	}

void MotifProfile::GetLettersFromSeq(const string &Seq,
  vector<uint> &Letters)
	{
	Letters.clear();
	asserta(SIZE(Seq) == AL+BL+CL);

	for (uint i = 0; i < MPL; ++i)
		{
		byte c = (byte) Seq[i];
		assert(c != 0);
		uint Letter = g_CharToLetterAmino[c];
		Letters.push_back(Letter);
		}
	}

void MotifProfile::FromXxxSeq(const string &Seq)
	{
	Clear();

	vector<uint> Letters;
	GetLettersFromXxxSeq(Seq, Letters);

	asserta(SIZE(Letters) == MPL);
	for (uint i = 0; i < MPL; ++i)
		{
		byte c = (byte) Seq[i];
		assert(c != 0);
		uint Letter = g_CharToLetterAmino[c];
		if (Letter >= 20)
			{
			for (uint j = 0; j < 20; ++j)
				m_FreqVec[i][j] = 1.0f/20.0f;
			continue;
			}
		m_FreqVec[i][Letter] = 1;
		}
	}

void MotifProfile::FromSeqs(const vector<string> &Seqs)
	{
	Clear();

	const uint SeqCount = SIZE(Seqs);
	asserta(SeqCount > 0);
	vector<uint> Letters;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Seq = Seqs[SeqIndex];
		GetLettersFromSeq(Seq, Letters);
		asserta(SIZE(Letters) == MPL);
		for (uint i = 0; i < MPL; ++i)
			{
			byte c = (byte) Seq[i];
			assert(c != 0);
			uint Letter = g_CharToLetterAmino[c];
			if (Letter >= 20)
				{
				for (uint j = 0; j < 20; ++j)
					m_FreqVec[i][j] += 1.0f/(20.0f*SeqCount);
				continue;
				}
			m_FreqVec[i][Letter] += 1.0f/SeqCount;
			}
		}
	ValidateFreqs();
	}

void MotifProfile::GetMaxLetter(uint i, uint &MaxLetter, float &MaxFreq) const
	{
	asserta(i < SIZE(m_FreqVec));
	const vector<float> &fs = m_FreqVec[i];
	MaxFreq = 0;
	MaxLetter = 0;
	for (uint j = 0; j < 20; ++j)
		{
		if (fs[j] > MaxFreq)
			{
			MaxFreq = fs[j];
			MaxLetter = j;
			}
		}
	}

void MotifProfile::GetLogo(string &Logo) const
	{
	Logo.clear();
	double Sum = 0;
	for (uint i = 0; i < MPL; ++i)
		{
		if (i == AL || i == AL+BL)
			Logo += "  ";
		uint Letter;
		float Freq;
		GetMaxLetter(i, Letter, Freq);
		char c = (char) g_LetterToCharAmino[Letter];
		if (Freq < 0.1)
			c = 'X';
		else if (Freq < 0.25)
			c = '.';
		else if (Freq < 0.5)
			c = tolower(c);
		Logo.push_back(c);
		}
	}

void MotifProfile::LogMe() const
	{
#if 0
	vector<uint> Order(20);
	for (uint i = 0; i < MPL; ++i)
		{
		const vector<float> &fs = m_FreqVec[i];
		QuickSortOrderDesc(fs.data(), 20, Order.data());
		Log("%c[%02u]  ", PosToMotif(i), PosToOffset(i));
		for (uint k = 0; k < 3; ++k)
			{
			uint Letter = Order[k];
			float f = fs[Letter];
			if (f < 1e-6)
				continue;
			if (k > 0)
				Log("  ");
			Log("%c/%.4f", g_LetterToCharAmino[Letter], fs[Letter]);
			}
		Log("\n");
		}
#endif
	string Logo;
	GetLogo(Logo);
	Log("%s\n", Logo.c_str());
	}

void MPCluster::GreedyCluster(const vector<MotifProfile *> &Input,
  float MinScore)
	{
	Clear();
	m_Input = &Input;
	m_MinScore = MinScore;
	const uint InputCount = SIZE(Input);
	for (uint i = 0; i < InputCount; ++i)
		m_PendingIndexes.insert(i);

	uint DoneCount = 0;
	vector<uint> ClusterSizes;
	while (m_PendingIndexes.size() > 0)
		{
		ProgressStep(DoneCount, InputCount + 1, "Clustering");
		uint SizeStart = SIZE(m_PendingIndexes);
		uint CentroidIndex = *m_PendingIndexes.begin();
		++DoneCount;
		m_PendingIndexes.erase(CentroidIndex);
		m_CentroidIndexes.push_back(CentroidIndex);

		MotifProfile &CP = GetProfile(CentroidIndex);
		vector<uint> MemberIndexes;
		MemberIndexes.push_back(CentroidIndex);
		for (set<uint>::const_iterator p = m_PendingIndexes.begin();
		  p != m_PendingIndexes.end(); ++p)
			{
			uint Index = *p;
			MotifProfile &P = GetProfile(Index);
			float Score = GetScore(CP, P);
			if (Score >= MinScore)
				{
				++DoneCount;
				MemberIndexes.push_back(Index);
				}
			}

		uint Size = SIZE(MemberIndexes);
		ClusterSizes.push_back(Size);
		m_CentroidIndexToMemberIndexes.push_back(MemberIndexes);
		for (uint i = 0; i < SIZE(MemberIndexes); ++i)
			{
			uint Index = MemberIndexes[i];
			m_PendingIndexes.erase(Index);
			}
		uint SizeEnd = SIZE(m_PendingIndexes);
		asserta(SizeEnd < SizeStart);
		}
	ProgressStep(InputCount, InputCount + 1, "Clustering");

	uint ClusterCount = SIZE(ClusterSizes);
	m_ClusterSizeOrder.resize(ClusterCount);
	QuickSortOrderDesc(ClusterSizes.data(), ClusterCount,
	  m_ClusterSizeOrder.data());
	}

void MPCluster::LogCluster(uint ClusterIndex) const
	{
	asserta(ClusterIndex < SIZE(m_CentroidIndexes));
	asserta(ClusterIndex < SIZE(m_CentroidIndexToMemberIndexes));

	uint CentroidIndex = m_CentroidIndexes[ClusterIndex];
	const vector<uint> &MemberIndexes = m_CentroidIndexToMemberIndexes[ClusterIndex];
	const uint Size = SIZE(MemberIndexes);

	Log("\n");
	Log("Cluster %u, size %u\n", ClusterIndex, Size);

	const MotifProfile &CP = GetProfile(CentroidIndex);
	string Logo;
	CP.GetLogo(Logo);
	Log("[Centroid]  %s\n", Logo.c_str());
	vector<float> Scores;
	for (uint i = 0; i < Size; ++i)
		{
		uint Index = MemberIndexes[i];
		const MotifProfile &P = GetProfile(Index);
		float Score = GetScore(CP, P);
		Scores.push_back(Score);
		}
	vector<uint> Order(Size);
	QuickSortOrderDesc(Scores.data(), Size, Order.data());

	for (uint k = 0; k < Size; ++k)
		{
		uint i = Order[k];
		uint Index = MemberIndexes[i];
		if (Index == CentroidIndex)
			continue;
		const MotifProfile &P = GetProfile(Index);
		float Score = GetScore(CP, P);
		P.GetLogo(Logo);
		Log("[%8.4f]  %s\n", Score, Logo.c_str());
		}
	}

void MPCluster::LogClusters() const
	{
	const uint ClusterCount = SIZE(m_CentroidIndexes);
	asserta(SIZE(m_CentroidIndexToMemberIndexes) == ClusterCount);
	Log("\n");
	Log("%u clusters\n", ClusterCount);
	for (uint i = 0; i < ClusterCount; ++i)
		LogCluster(i);
	}

void MPCluster::ClusterToTsv(FILE *f, uint OrderIndex) const
	{
	if (f == 0)
		return;

	asserta(OrderIndex < SIZE(m_ClusterSizeOrder));
	uint ClusterIndex = m_ClusterSizeOrder[OrderIndex];
	asserta(ClusterIndex < SIZE(m_CentroidIndexes));
	asserta(ClusterIndex < SIZE(m_CentroidIndexToMemberIndexes));

	uint CentroidIndex = m_CentroidIndexes[ClusterIndex];
	const vector<uint> &MemberIndexes = m_CentroidIndexToMemberIndexes[ClusterIndex];
	const uint Size = SIZE(MemberIndexes);

	const MotifProfile &CP = GetProfile(CentroidIndex);
	string CentroidLogo;
	CP.GetLogo(CentroidLogo);
	fprintf(f, "C\t%u\t%u\t%s\n", 
	  OrderIndex, Size, CentroidLogo.c_str());

	vector<float> Scores;
	for (uint i = 0; i < Size; ++i)
		{
		uint Index = MemberIndexes[i];
		const MotifProfile &P = GetProfile(Index);
		float Score = GetScore(CP, P);
		Scores.push_back(Score);
		}
	vector<uint> Order(Size);
	QuickSortOrderDesc(Scores.data(), Size, Order.data());

	for (uint k = 0; k < Size; ++k)
		{
		uint i = Order[k];
		uint Index = MemberIndexes[i];
		if (Index == CentroidIndex)
			continue;
		const MotifProfile &P = GetProfile(Index);
		float Score = GetScore(CP, P);
		string Logo;
		P.GetLogo(Logo);
		fprintf(f, "M\t%u\t%.2f\t%s\n",
		  OrderIndex, Score, Logo.c_str());
		}
	}

void MPCluster::ClustersToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	const uint ClusterCount = SIZE(m_CentroidIndexes);
	asserta(SIZE(m_CentroidIndexToMemberIndexes) == ClusterCount);
	for (uint i = 0; i < ClusterCount; ++i)
		ClusterToTsv(f, i);
	}

void MPCluster::FindNN(uint &Index1, uint &Index2) const
	{
	asserta(SIZE(m_PendingIndexes) >= 2);
	Index1 = UINT_MAX;
	Index2 = UINT_MAX;
	float MaxScore = -1;
	for (set<uint>::const_iterator p = m_PendingIndexes.begin();
	  p != m_PendingIndexes.end(); ++p)
		{
		uint i1 = *p;
		set<uint>::const_iterator q = p;
		for (;;)
			{
			++q;
			if (q == m_PendingIndexes.end())
				break;
			uint i2 = *q;
			float Score = GetScoreNNPair(i1, i2);
			if (Score > MaxScore)
				{
				MaxScore = Score;
				Index1 = i1;
				Index2 = i2;
				}
			}
		}
	}

void MPCluster::Join(uint Index1, uint Index2)
	{
	MotifProfile &P = CreateProfileNN(Index1, Index2);
	uint Parent = SIZE(m_MPs);
	uint LeftSize = m_Sizes[Index1];
	uint RightSize = m_Sizes[Index2];
	uint ParentSize = LeftSize + RightSize;
	m_MPs.push_back(&P);
	m_Parents.push_back(Parent);
	m_Lefts.push_back(Index1);
	m_Rights.push_back(Index2);
	m_PendingIndexes.erase(Index1);
	m_PendingIndexes.erase(Index2);
	m_PendingIndexes.insert(Parent);
	m_Sizes.push_back(ParentSize);
	}

float MPCluster::GetScoreNNPair(uint i1, uint i2) const
	{
	asserta(i1 < SIZE(m_MPs));
	asserta(i2 < SIZE(m_MPs));

	const MotifProfile &MP1 = *m_MPs[i1];
	const MotifProfile &MP2 = *m_MPs[i2];
	float Score = GetScore(MP1, MP2);
	return Score;
	}

MotifProfile &MPCluster::CreateProfileNN(uint i1, uint i2) const
	{
	asserta(i1 < SIZE(m_MPs));
	asserta(i2 < SIZE(m_MPs));

	const MotifProfile &MP1 = *m_MPs[i1];
	const MotifProfile &MP2 = *m_MPs[i2];

	uint Size1 = m_Sizes[i1];
	uint Size2 = m_Sizes[i2];
	float w1 = (float) sqrt(Size1);
	float w2 = (float) sqrt(Size2);
	float w12 = w1 + w2;
	w1 /= w12;
	w2 /= w12;
	asserta(feq(w1 + w2, 1.0));

	MotifProfile &P = *new MotifProfile;
	for (uint i = 0; i < MPL; ++i)
		{
		const vector<float> &v1 = MP1.m_FreqVec[i];
		const vector<float> &v2 = MP2.m_FreqVec[i];
		vector<float> &v = P.m_FreqVec[i];

		float Sum = 0;
		for (uint j = 0; j < 20; ++j)
			{
			float f1 = v1[j];
			float f2 = v2[j];
			float f = w1*f1 + w2*f2;
			v[j] = f;
			Sum += f;
			}
		asserta(feq(Sum, 1.0f));
		}

	string Logo1;
	string Logo2;
	string LogoP;
	MP1.GetLogo(Logo1);
	MP2.GetLogo(Logo2);
	P.GetLogo(LogoP);
	
	Log("\n");
	Log("Join\n");
	Log("  L[%5u, %6.4f]  %s\n", Size1, w1, Logo1.c_str());
	Log("  R[%5u, %6.4f]  %s\n", Size2, w2, Logo2.c_str());

	return P;
	}

void MPCluster::NNCluster(const vector<MotifProfile *> &Input,
  float MinScore)
	{
	Clear();
	Clear();
	m_Input = &Input;
	const uint N = SIZE(Input);
	for (uint i = 0; i < N; ++i)
		{
		m_MPs.push_back(Input[i]);
		m_PendingIndexes.insert(i);
		m_Sizes.push_back(1);
		}

	uint JoinCount = N - 1;
	for (uint JoinIndex = 0; JoinIndex < JoinCount; ++JoinIndex)
		{
		ProgressStep(JoinIndex, JoinCount, "NN cluster");
		uint Index1, Index2;
		FindNN(Index1, Index2);
		Join(Index1, Index2);
		}
	}

static void ClusterMotifs(const string &InputFileName,
  const string &Strategy)
	{
	SetBLOSUM62();

	asserta(optset_minscore);
	const float MinScore = (float) opt_minscore;

	SeqDB Input;
	Input.FromFasta(InputFileName);
	const uint InputCount = Input.GetSeqCount();
	vector<MotifProfile *> MPs;
	for (uint i = 0; i < InputCount; ++i)
		{
		const string &Seq = Input.GetSeq(i);
		MotifProfile *MP = new MotifProfile;
		MP->FromXxxSeq(Seq);
		MPs.push_back(MP);
		}

	MPCluster MC;
	if (Strategy == "greedy")
		{
		MC.GreedyCluster(MPs, MinScore);
		MC.LogClusters();
		ProgressLog("%u motifs, %u clusters\n",
		  InputCount, SIZE(MC.m_CentroidIndexes));
		MC.ClustersToTsv(g_ftsv);
		}
	else if (Strategy == "nn")
		MC.NNCluster(MPs, MinScore);
	else
		Die("Invalid strategy '%s'", Strategy.c_str());
	}

void cmd_cluster_motifs_greedy()
	{
	const string &InputFileName = opt_cluster_motifs_greedy;
	ClusterMotifs(InputFileName, "greedy");
	}

void cmd_cluster_motifs_nn()
	{
	const string &InputFileName = opt_cluster_motifs_nn;
	ClusterMotifs(InputFileName, "nn");
	}

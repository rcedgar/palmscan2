#include "myutils.h"
#include "abcxyz.h"
#include "filler.h"
#include "motifsettings.h"
#include "quarts.h"

void Filler::FromLines(const vector<string> &Lines)
	{
	InitMx(m_ScoreMx);
	asserta(Lines[0] == "scoremx	80	69	57");
	uint LineNr = 1;
	vector<string> Fields;
	for (int ix = 0; ix < m_RangeX; ++ix)
		{
		for (int iy = 0; iy < m_RangeY; ++iy)
			{
			const string &Line = Lines[LineNr++];
			Split(Line, Fields, '\t');
			asserta(SIZE(Fields) == (uint) m_RangeZ);
			for (int iz = 0; iz < m_RangeZ; ++iz)
				{
				double Score = StrToFloat(Fields[iz]);
				m_ScoreMx[ix][iy][iz] = Score;
				}
			}
		}
	}

void Filler::ScoreMxToFile(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, "scoremx\t%d\t%d\t%d\n", m_RangeX, m_RangeY, m_RangeZ);
	for (int ix = 0; ix < m_RangeX; ++ix)
		{
		for (int iy = 0; iy < m_RangeY; ++iy)
			{
			for (int iz = 0; iz < m_RangeZ; ++iz)
				{
				double Score = m_ScoreMx[ix][iy][iz];
				if (iz > 0)
					fprintf(f, "\t");
				fprintf(f, "%.3g", Score);
				}
			fprintf(f, "\n");
			}
		}
	}

double Filler::CalcZeroScoreFromCounts(double Base,
  int ix, int iy, int iz) const
	{
	const int d = m_Radius;
	uint ZeroCount = 0;
	uint NonZeroCount = 0;
	for (int jx = ix - d; jx <= ix + d; ++jx)
		{
		if (jx < 0 || jx >= m_RangeX)
			continue;
		for (int jy = iy - d; jy <= iy + d; ++jy)
			{
			if (jy < 0 || jy >= m_RangeY)
				continue;
			for (int jz = iz - d; jz <= iz + d; ++jz)
				{
				if (jz < 0 || jz >= m_RangeZ)
					continue;
				int totd = abs(jx - ix) + abs(jy - iy) + abs(jz - iz);
				if (totd > 3)
					continue;
				uint n = m_CountMx[jx][jy][jz];
				if (n == 0)
					++ZeroCount;
				else
					++NonZeroCount;
				}
			}
		}
	if (ZeroCount > 1 && NonZeroCount == 0)
		return Base;
	return DBL_MAX;
	}

double Filler::CalcScoreFromCounts(int ix, int iy, int iz) const
	{
	const int d = m_Radius;
	double Sum = 0;
	uint ZeroCount = 0;
	uint NonZeroCount = 0;
	for (int jx = ix - d; jx <= ix + d; ++jx)
		{
		if (jx < 0 || jx >= m_RangeX)
			continue;
		for (int jy = iy - d; jy <= iy + d; ++jy)
			{
			if (jy < 0 || jy >= m_RangeY)
				continue;
			for (int jz = iz - d; jz <= iz + d; ++jz)
				{
				if (jz < 0 || jz >= m_RangeZ)
					continue;
				int totd = abs(jx - ix) + abs(jy - iy) + abs(jz - iz);
				if (totd > 3)
					continue;
				uint n = m_CountMx[jx][jy][jz];
				if (n == 0)
					++ZeroCount;
				else
					++NonZeroCount;
				double f = double(n)/(1 + totd);
				Sum += f;
				}
			}
		}
	double Score = 100*Sum/(ZeroCount + NonZeroCount);
	return Score;
	}

double Filler::CalcFilledLocMx1(int ix, int iy, int iz) const
	{
	const int d = m_Radius;
	double Sum = 0;
	uint N = 0;
	for (int jx = ix - d; jx <= ix + d; ++jx)
		{
		if (jx < 0 || jx >= m_RangeX)
			continue;
		for (int jy = iy - d; jy <= iy + d; ++jy)
			{
			if (jy < 0 || jy >= m_RangeY)
				continue;
			for (int jz = iz - d; jz <= iz + d; ++jz)
				{
				if (jz < 0 || jz >= m_RangeZ)
					continue;
				int totd = abs(jx - ix) + abs(jy - iy) + abs(jz - iz);
				if (totd > 3)
					continue;
				uint n = m_CountMx[jx][jy][jz];
				double Loc = m_LocMx[jx][jy][jz];
				if (n == 0)
					{
					asserta(Loc == 0);
					continue;
					}
				assert(Loc > 0);
				N += n;
				Sum += Loc;
				}
			}
		}
	if (N == 0)
		return 0;
	double FilledLoc = Sum/N;
	return FilledLoc;
	}

void Filler::CalcFilledLocMx()
	{
	InitMx(m_FilledLocMx);
	for (int ix = 0; ix < m_RangeX; ++ix)
		for (int iy = 0; iy < m_RangeY; ++iy)
			for (int iz = 0; iz < m_RangeZ; ++iz)
				{
				double Loc = CalcFilledLocMx1(ix, iy, iz);
				m_FilledLocMx[ix][iy][iz] = Loc;
				}
	}

void Filler::CalcScoreMx()
	{
	InitMx(m_ScoreMx);
	vector<double> PositiveScores;
	for (int ix = 0; ix < m_RangeX; ++ix)
		for (int iy = 0; iy < m_RangeY; ++iy)
			for (int iz = 0; iz < m_RangeZ; ++iz)
				{
				double Score = CalcScoreFromCounts(ix, iy, iz);
				m_ScoreMx[ix][iy][iz] = Score;
				if (Score > 0)
					PositiveScores.push_back(Score);
				}

	QuartsDouble Q;
	GetQuartsDouble(PositiveScores, Q);
	Log("Positive scores ");
	Q.LogMe();

	double Med = Q.Med;
	for (int ix = 0; ix < m_RangeX; ++ix)
		for (int iy = 0; iy < m_RangeY; ++iy)
			for (int iz = 0; iz < m_RangeZ; ++iz)
				m_ScoreMx[ix][iy][iz] /= Med;

	uint ZN = 0;
	uint N = 0;
	double ZeroScore = -2*Med;
	for (int ix = 0; ix < m_RangeX; ++ix)
		for (int iy = 0; iy < m_RangeY; ++iy)
			for (int iz = 0; iz < m_RangeZ; ++iz)
				{
				++N;
				double Score = CalcZeroScoreFromCounts(ZeroScore, ix, iy, iz);
				if (Score != DBL_MAX)
					{
					m_ScoreMx[ix][iy][iz] = Score;
					++ZN;
					}
				}
	Log("ZN=%u/%u (%.1f%%)\n", ZN, N, GetPct(ZN, N));
	}

void Filler::UpdateCountMx1(int xn, int yn, int zn)
	{
	if (xn < 0 || xn >= m_RangeX)
		return;
	if (yn < 0 || yn >= m_RangeY)
		return;
	if (zn < 0 || zn >= m_RangeZ)
		return;
	m_CountMx[xn][yn][zn] += 1;
	}

void Filler::UpdateLocMx1(int xn, int yn, int zn, double Loc)
	{
	if (xn < 0 || xn >= m_RangeX)
		return;
	if (yn < 0 || yn >= m_RangeY)
		return;
	if (zn < 0 || zn >= m_RangeZ)
		return;
	asserta(Loc > 0);
	m_LocMx[xn][yn][zn] += Loc;
	}

void Filler::UpdateCountMx(const PDBChain &PPX)
	{
	++m_TrainCount;
	const uint L = PPX.GetSeqLength();
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		double x, y, z;
		PPX.GetXYZ(Pos, x, y, z);

		int ix = int(x+0.5);
		int iy = int(y+0.5);
		int iz = int(z+0.5);
		int xn = m_MaxX - ix;
		int yn = m_MaxY - iy;
		int zn = m_MaxZ - iz;
		UpdateCountMx1(xn, yn, zn);
		}
	}

void Filler::UpdateLocMx(const PDBChain &PPX)
	{
	const uint L = PPX.GetSeqLength();
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		double x, y, z;
		PPX.GetXYZ(Pos, x, y, z);
		double Loc = (Pos+1)*100.0/L;

		int ix = int(x+0.5);
		int iy = int(y+0.5);
		int iz = int(z+0.5);
		int xn = m_MaxX - ix;
		int yn = m_MaxY - iy;
		int zn = m_MaxZ - iz;
		UpdateLocMx1(xn, yn, zn, Loc);
		}
	}

double Filler::GetPPScore(const PDBChain &PPX) const
	{
	const uint L = PPX.GetSeqLength();
	double Sum = 0;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		double x, y, z;
		PPX.GetXYZ(Pos, x, y, z);

		int ix = int(x+0.5);
		int iy = int(y+0.5);
		int iz = int(z+0.5);
		int xn = m_MaxX - ix;
		int yn = m_MaxY - iy;
		int zn = m_MaxZ - iz;
		if (xn < 0 || xn >= m_RangeX)
			{
			Sum += m_OverflowScore;
			continue;
			}
		if (yn < 0 || yn >= m_RangeY)
			{
			Sum += m_OverflowScore;
			continue;
			}
		if (zn < 0 || zn >= m_RangeZ)
			{
			Sum += m_OverflowScore;
			continue;
			}
		double Score = m_ScoreMx[xn][yn][zn];
		Sum += Score;
		}
	double PPScore = Sum/L;
	return PPScore;
	}

double Filler::GetPPLocScore(const PDBChain &PPX) const
	{
	const uint L = PPX.GetSeqLength();
	double Sum = 0;
	double SumDiff = 0;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		double Loc = (Pos+1)*100.0/L;

		double x, y, z;
		PPX.GetXYZ(Pos, x, y, z);

		int ix = int(x+0.5);
		int iy = int(y+0.5);
		int iz = int(z+0.5);
		int xn = m_MaxX - ix;
		int yn = m_MaxY - iy;
		int zn = m_MaxZ - iz;
		if (xn < 0 || xn >= m_RangeX)
			continue;
		if (yn < 0 || yn >= m_RangeY)
			continue;
		if (zn < 0 || zn >= m_RangeZ)
			continue;
		double FilledLoc = m_FilledLocMx[xn][yn][zn];
		if (FilledLoc == 0)
			SumDiff += 50;
		else
			SumDiff += fabs(Loc - FilledLoc);
		}
	double PPLocScore = SumDiff/L;
	return PPLocScore;
	}

void Filler::InitMx(vector<vector<vector<uint> > > &Mx)
	{
	Mx.clear();
	Mx.resize(m_RangeX);
	for (int x = 0; x < m_RangeX; ++x)
		{
		Mx[x].resize(m_RangeY);
		for (int y = 0; y < m_RangeY; ++y)
			Mx[x][y].resize(m_RangeZ, 0);
		}
	}

void Filler::InitMx(vector<vector<vector<double> > > &Mx)
	{
	Mx.clear();
	Mx.resize(m_RangeX);
	for (int x = 0; x < m_RangeX; ++x)
		{
		Mx[x].resize(m_RangeY);
		for (int y = 0; y < m_RangeY; ++y)
			Mx[x][y].resize(m_RangeZ, 0);
		}
	}

void Filler::Train(const vector<PDBChain *> &Chains)
	{
	InitMx(m_CountMx);
	InitMx(m_LocMx);

	GetRdrpModel(m_RdRpModel);
	m_RS.Init(m_RdRpModel);

	m_TrainCount = 0;
	m_PermutedCount = 0;
	m_NoPSSMHitCount = 0;

	const uint ChainCount = SIZE(Chains);
	asserta(ChainCount > 0);
	for (uint i = 0; i < ChainCount; ++i)
		{
		ProgressStep(i, ChainCount, "Training");
		const PDBChain &Chain = *Chains[i];
		Train1(Chain);
		}
	CalcScoreMx();
	CalcFilledLocMx();

	vector<double> Scores;
	vector<double> RevScores;
	for (uint i = 0; i < ChainCount; ++i)
		{
		ProgressStep(i, ChainCount, "Testing");
		const PDBChain &Chain = *Chains[i];

		double Score, RevScore;
		Test1(Chain, Score, RevScore);
		Log("%8.3g  %8.3g  >%s\n",
		  Score, RevScore, Chain.m_Label.c_str());
		Scores.push_back(Score);
		RevScores.push_back(RevScore);
		}
	QuartsDouble Q;
	GetQuartsDouble(Scores, Q);
	Log("Pos: ");
	Q.LogMe();
	GetQuartsDouble(RevScores, Q);
	Log("Rev: ");
	Q.LogMe();
	}

void Filler::GetPPX(const PDBChain &Chain,
  uint PosA, uint PosB, uint PosC,
  PDBChain &PPX) const
	{
	const uint L = Chain.GetSeqLength();
	asserta(PosA < PosC && PosC + g_LC <= L);
	PDBChain PP;
	Chain.GetRange(PosA, PosC + g_LC - 1, PP);

	uint PPPosB = PosB - PosA;
	uint PPPosC = PosC - PosA;
	uint PPPosA = 0;

	vector<double> PtC1;
	vector<double> PtC2;
	vector<double> PtA1;
	PP.GetPt(PPPosC + g_OffCd - 5, PtC1);
	PP.GetPt(PPPosC + g_OffCd - 2, PtC2);
	PP.GetPt(PPPosA + g_OffAd, PtA1);

	vector<double> Origin = PtC2;

	vector<double> Axis0;
	Sub_Vecs(PtC1, PtC2, Axis0);
	NormalizeVec(Axis0);

	vector<double> CtoA;
	Sub_Vecs(PtC2, PtA1, CtoA);
	NormalizeVec(CtoA);

	vector<double> Axis2;
	CrossProduct(Axis0, CtoA, Axis2);
	NormalizeVec(Axis2);

	vector<double> Axis1;
	CrossProduct(Axis0, Axis2, Axis1);
	NormalizeVec(Axis1);

	vector<vector<double> > Basis;
	Basis.push_back(Axis0);
	Basis.push_back(Axis1);
	Basis.push_back(Axis2);

	AssertUnitBasisA(Basis);

	vector<vector<double> > R;
	GetBasisR(Basis, R);

	vector<double> t;
	t.push_back(-Origin[0]);
	t.push_back(-Origin[1]);
	t.push_back(-Origin[2]);

	PPX = *new PDBChain;
	PP.GetXFormChain_tR(t, R, PPX);
	}

void Filler::Train1(const PDBChain &Chain)
	{
	const string &Label = Chain.m_Label;
	const string &Seq = Chain.m_Seq;
	m_RS.Search(Label, Seq);

	string PSSM_A, PSSM_B, PSSM_C;
	m_RS.GetShapesTrainABC(PSSM_A, PSSM_B, PSSM_C);
		
	uint L = Chain.GetSeqLength();

	size_t stPosA = Seq.find(PSSM_A);
	size_t stPosB = Seq.find(PSSM_B);
	size_t stPosC = Seq.find(PSSM_C);
	if (stPosA == string::npos || stPosB == string::npos ||
		stPosC == string::npos)
		{
		++m_NoPSSMHitCount;
		return;
		}

	asserta(stPosA < L && stPosB < L && stPosC < L);

	uint PosA = uint(stPosA);
	uint PosB = uint(stPosB);
	uint PosC = uint(stPosC);

	if (PosC < PosA)
		{
		++m_PermutedCount;
		return;
		}

	PDBChain &PPX = *new PDBChain;
	GetPPX(Chain, PosA, PosB, PosC, PPX);

	UpdateCountMx(PPX);
	UpdateLocMx(PPX);
	}

void Filler::Test1(const PDBChain &Chain,
  double &Score, double &RevScore)
	{
	const string &Label = Chain.m_Label;
	const string &Seq = Chain.m_Seq;
	m_RS.Search(Label, Seq);

	string PSSM_A, PSSM_B, PSSM_C;
	m_RS.GetShapesTrainABC(PSSM_A, PSSM_B, PSSM_C);
		
	uint L = Chain.GetSeqLength();

	size_t stPosA = Seq.find(PSSM_A);
	size_t stPosB = Seq.find(PSSM_B);
	size_t stPosC = Seq.find(PSSM_C);
	if (stPosA == string::npos || stPosB == string::npos ||
		stPosC == string::npos)
		return;

	asserta(stPosA < L && stPosB < L && stPosC < L);

	uint PosA = uint(stPosA);
	uint PosB = uint(stPosB);
	uint PosC = uint(stPosC);

	if (PosC < PosA)
		return;

	Score = GetScore(Chain, PosA, PosB, PosC);

	PDBChain Rev;
	Chain.GetReverse(Rev);
	RevScore = GetScore(Rev, PosA, PosB, PosC);

	PDBChain PPX;
	GetPPX(Chain, PosA, PosB, PosC, PPX);
	double LocScore = GetPPLocScore(PPX);
	double RevLocScore = GetPPLocScore(Rev);
	Log("LocScore\t%.3g\t%.3g\t%.3g\t%.3g\n",
	  Score, RevScore, LocScore, RevLocScore);
	}

double Filler::GetScore(const PDBChain &Chain,
  uint PosA, uint PosB, uint PosC) const
	{
	if (PosA == UINT_MAX || PosB == UINT_MAX || PosC == UINT_MAX)
		return 0;
	if (PosC < PosA)
		Die("Filler::GetScorePermuted");

	PDBChain PPX;
	GetPPX(Chain, PosA, PosB, PosC, PPX);

	double Score = GetPPScore(PPX);
	return Score;
	}

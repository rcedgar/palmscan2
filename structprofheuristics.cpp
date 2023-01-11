#include "myutils.h"
#include "structprof.h"
#include "abcxyz.h"

/***
D1 = max dist to G in 50aa window after motif C
  (start of pre-D helix).
D2 = min dist to G in 50aa window after D1.
D3 = min dist to A in +/- 10aa window around D2.
***/
uint StructProf::FindMofifD_Hueuristics() const
	{
	const uint L = m_Chain->GetSeqLength();

	uint Pos_aD = m_Chain->GetMotifPos(A) + 3;
	uint Pos_bG = m_Chain->GetMotifPos(B) + 1;
	uint Pos_cD = m_Chain->GetMotifPos(C) + 3;

	double Dist1 = DBL_MAX;
	uint Lo1 = Pos_cD;
	uint Hi1 = Lo1 + 50;
	if (Hi1 >= L)
		Hi1 = L - 1;
	uint D1 = SearchDist(Pos_bG, Lo1, Hi1, true, 5.0, Dist1);
	Log("Lo1 = %u, Hi1 = %u, Dist1=%.1f, D1=%u\n",
	  Lo1, Hi1, Dist1, D1);//@@

	double Dist2 = DBL_MAX;
	uint Lo2 = D1;
	uint Hi2 = D1 + 50;
	if (Hi2 >= L)
		Hi2 = L - 1;
	uint D2 = SearchDist(Pos_bG, Lo2, Hi2, false, 5.0, Dist2, true);
	Log("Lo2 = %u, Hi2 = %u, Dist2=%.1f, D2=%u\n",
	  Lo2, Hi2, Dist2, D2);//@@

	uint CN = GetCavityNumber(D2);
	if (CN > 6)
		return UINT_MAX;

	Log("D[%u] = ", D2+1);
	for (int i = -3; i <= 3; ++i)
		{
		char c = m_Chain->m_Seq[int(D2)+i];
		Log("%c", c);
		}
	Log("\n");

	return D2;
	}

uint StructProf::FindMofifE_Hueuristics(uint Pos_MotifD) const
	{
	const uint L = m_Chain->GetSeqLength();
	uint Pos_bG = m_Chain->GetMotifPos(B) + 1;
	uint Lo1 = Pos_MotifD + 8;
	uint Hi1 = Lo1 + 25;
	if (Hi1 >= L)
		Hi1 = L - 1;
	double Dist1 = DBL_MAX;
	uint Pos1 = SearchDist(Pos_bG, Lo1, Hi1, true, 5.0, Dist1);
	if (Pos1 >= Hi1 - 5)
		return UINT_MAX;

	double DistE;
	uint PosE = SearchDist(Pos_bG, Pos1+1, Hi1, false, 5.0, DistE);
	Log("Lo1 = %u, Hi1 = %u, Dist1=%.1f, PosE=%u\n",
	  Lo1, Hi1, Dist1, PosE);//@@
	Log("E[%u] = ", PosE+1);
	for (int i = -3; i <= 3; ++i)
		{
		char c = m_Chain->m_Seq[int(PosE)+i];
		Log("%c", c);
		}
	Log("\n");
	return PosE;
	}

uint StructProf::FindMofifF_Hueuristics(uint Pos_MotifA) const
	{
	if (Pos_MotifA < 25)
		return UINT_MAX;
	uint Pos_bG = m_Chain->GetMotifPos(B) + 1;
	double Dist1 = DBL_MAX;
	uint Lo1 = (Pos_MotifA < 100 ? 0 : Pos_MotifA - 100);
	uint Hi1 = Pos_MotifA - 10;
	uint PosF = SearchDist(Pos_bG, Lo1, Hi1, false, 99.0, Dist1);
	Log("Lo1 = %u, Hi1 = %u, Dist1=%.1f, PosF=%u\n",
	  Lo1, Hi1, Dist1, PosF);//@@
	Log("F[%u] = ", PosF+1);
	for (int i = -3; i <= 3; ++i)
		{
		char c = m_Chain->m_Seq[int(PosF)+i];
		Log("%c", c);
		}
	Log("\n");
	return PosF;
	}

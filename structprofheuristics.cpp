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
	uint Pos_aD = m_Chain->GetMotifPos(A) + 3;
	uint Pos_bG = m_Chain->GetMotifPos(B) + 1;
	uint Pos_cD = m_Chain->GetMotifPos(C) + 3;

	double Dist1 = DBL_MAX;
	uint Lo1 = Pos_cD;
	uint Hi1 = Lo1 + 50;
	uint D1 = SearchDist(Pos_bG, Lo1, Hi1, true, 5.0, Dist1);
	Log("Lo1 = %u, Hi1 = %u, Dist1=%.1f, D1=%u\n",
	  Lo1, Hi1, Dist1, D1);//@@

	double Dist2 = DBL_MAX;
	uint Lo2 = D1;
	uint Hi2 = D1 + 50;
	uint D2 = SearchDist(Pos_bG, Lo2, Hi2, false, 5.0, Dist2);
	Log("Lo2 = %u, Hi2 = %u, Dist2=%.1f, D2=%u\n",
	  Lo2, Hi2, Dist2, D2);//@@

	double Dist3 = DBL_MAX;
	uint Lo3 = D2 - 10;
	uint Hi3 = D2 + 10;
	uint D3 = SearchDist(Pos_aD, Lo3, Hi3, false, 20.0, Dist3);
	Log("Lo3 = %u, Hi3 = %u, Dist3=%.1f, D3=%u\n",
	  Lo3, Hi3, Dist3, D3);//@@

	Log("D[%u] = ", D3+1);
	for (int i = -3; i <= 3; ++i)
		{
		char c = m_Chain->m_Seq[int(D3)+i];
		Log("%c", c);
		}
	Log("\n");

	return D3;
	}

uint StructProf::FindMofifE_Hueuristics(uint Pos_MotifD) const
	{
	uint Pos_bG = m_Chain->GetMotifPos(B) + 1;
	double Dist1 = DBL_MAX;
	uint Lo1 = Pos_MotifD + 8;
	uint Hi1 = Lo1 + 25;
	uint PosE = SearchDist(Pos_bG, Lo1, Hi1, false, 5.0, Dist1);
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

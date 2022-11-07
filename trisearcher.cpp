#include "myutils.h"
#include "trisearcher.h"
#include "pdb.h"
#include "sort.h"
#include "abcxyz.h"

void TriSearcher::LogMe(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, "TriSearcher() Radius %.1f\n", Radius);
	fprintf(f, "NAB %u..%u", NABmin, NABmax);
	fprintf(f, ", NBC %u..%u", NBCmin, NBCmax);
	fprintf(f, ", NAC %u..%u\n", NACmin, NACmax);

	const uint N = SIZE(m_RMSDs);
	if (N > 0)
		{
		vector<uint> Order;
		GetOrder(Order);
		Log("%u hits\n", N);
		Log(" PosA   PosB   PosC     RMSD\n");
		//   12345  12345  12345  1234567  
		for (uint k = 0; k < N; ++k)
			{
			uint i = Order[k];
			uint PosA = m_PosAs[i];
			uint PosB = m_PosBs[i];
			uint PosC = m_PosCs[i];
			double RMSD = m_RMSDs[i];

			string A, B, C;
			m_Query->GetMotifSeq(PosA, AL, false, A);
			m_Query->GetMotifSeq(PosB, BL, false, B);
			m_Query->GetMotifSeq(PosC, CL, false, C);

			Log("%5u", PosA);
			Log("  %5u", PosB);
			Log("  %5u", PosC);
			Log("  %7.2f", RMSD);
			Log("  %s", A.c_str());
			Log("  %s", B.c_str());
			Log("  %s", C.c_str());
			Log("\n");
			}
		}
	}

void TriSearcher::GetOrder(vector<uint> &Order) const
	{
	const uint N = SIZE(m_PosAs);
	if (N == 0)
		{
		Order.clear();
		return;
		}
	Order.resize(N);
	QuickSortOrder(m_RMSDs.data(), N, Order.data());
	}

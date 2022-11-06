#include "myutils.h"
#include "pdb.h"

#define TRACE	0

void PDB::SearchTriangle(const TriParams &TP,
  vector<uint> &PosAs,
  vector<uint> &PosBs,
  vector<uint> &PosCs,
  vector<double> &RMSDs)
	{
	asserta(TP.NABmin <= TP.NABmax);
	asserta(TP.NBCmin <= TP.NBCmax);
	asserta(TP.NACmin <= TP.NACmax);

	PosAs.clear();
	PosBs.clear();
	PosCs.clear();
	RMSDs.clear();

	const uint L = GetSeqLength();
#if TRACE
	Log("\n");
	Log("SearchTriangle()\n");
	Log("L=%u\n", L);
#endif
	for (uint PosA = 0; PosA < L; ++PosA)
		{
		if (PosA + TP.NACmin > L)
			{
#if TRACE
			Log("PosA=%u + NACmin > L\n", PosA);
#endif
			break;
			}

		for (uint PosC = PosA + TP.NACmin; PosC <= PosA + TP.NACmax; ++PosC)
			{
			if (PosC >= L)
				{
#if TRACE
				Log("PosC=%u > L\n", PosC);
#endif
				break;
				}

			double AC = GetDist(PosA, PosC);
#if TRACE
			Log("PosA=%u, PosC=%u, AC=%.1f LAC=%.1f\n", PosA, PosC, AC, TP.LAC);
#endif
			if (AC < TP.LAC - TP.Radius || AC > TP.LAC + TP.Radius)
				continue;

			for (uint PosB = PosA + TP.NABmin; PosB <= PosA + TP.NABmax; ++PosB)
				{
				if (PosB >= PosC)
					continue;
				uint NBC = PosC - PosB;
				if (NBC < TP.NBCmin || NBC > TP.NBCmax)
					continue;

				double AB = GetDist(PosA, PosB);
#if TRACE
				Log("PosA=%u, PosB=%u, AB=%.1f LAB=%.1f\n", PosA, PosB, AB, TP.LAB);
#endif
				if (AB < TP.LAB - TP.Radius || AB > TP.LAB + TP.Radius)
					continue;

				double BC = GetDist(PosB, PosC);
#if TRACE
				Log("PosB=%u, PosC=%u, BC=%.1f LBC=%.1f\n", PosB, PosC, BC, TP.LBC);
#endif
				if (BC < TP.LBC - TP.Radius || BC > TP.LBC + TP.Radius)
					continue;

				double dAB = TP.LAB - AB;
				double dBC = TP.LBC - BC;
				double dAC = TP.LAC - AC;
				double d2 = dAB*dAB + dBC*dBC + dAC*dAC;
				double RMSD = sqrt(d2);

				PosAs.push_back(PosA);
				PosBs.push_back(PosB);
				PosCs.push_back(PosC);
				RMSDs.push_back(RMSD);
				}
			}
		}
	}

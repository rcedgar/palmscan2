#include "myutils.h"
#include "pdbchain.h"
#include "ppcaligner.h"
#include "calreader.h"
#include "outputfiles.h"

void cmd_classify()
	{
	const string &QFN = opt_classify;
	const string &RefFN = opt_ref;

	vector<PDBChain *> Rs;
	ReadChains(RefFN, Rs);
	const uint NR = SIZE(Rs);
	asserta(NR > 0);

	vector<bool> RefIsRdRp;
	for (uint i = 0; i < NR; ++i)
		{
		const string &Label = Rs[i]->m_Label;
		bool IsRdRp = StartsWith(Label, "rdrp.");
		RefIsRdRp.push_back(IsRdRp);
		Log("%c  %s\n", pom(IsRdRp), Label.c_str());
		}

	PPCAligner PA;

	CalReader CR;
	CR.Open(QFN);

	uint RecCount = 0;
	uint HitCount = 0;
	PDBChain Q;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;
		if (RecCount%10000 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u RdRp hits\r",
			  sPct.c_str(), HitCount, RecCount);
			}
		++RecCount;
		if (Q.m_MotifPosVec.empty())
			{
			Warning("Missing motif coords >%s", Q.m_Label.c_str());
			continue;
			}
		PA.SetQuery(Q);

		double RMSD_TopRdRp = DBL_MAX;
		double RMSD_TopOther = DBL_MAX;
		uint RefIndex_TopRdRp = UINT_MAX;
		uint RefIndex_TopOther = UINT_MAX;
		for (uint RefIndex = 0; RefIndex < NR; ++RefIndex)
			{
			const PDBChain &R = *Rs[RefIndex];
			PA.SetRef(R);
			double RMSD = PA.GetMotifRMSD();
			bool IsRdRp = RefIsRdRp[RefIndex];
			if (IsRdRp)
				{
				if (RMSD < RMSD_TopRdRp)
					{
					RMSD_TopRdRp = RMSD;
					RefIndex_TopRdRp = RefIndex;
					}
				}
			else
				{
				if (RMSD < RMSD_TopOther)
					{
					RMSD_TopOther = RMSD;
					RefIndex_TopOther = RefIndex;
					}
				}
			}

		if (g_ftsv == 0)
			continue;

		if (RefIndex_TopRdRp == UINT_MAX && RefIndex_TopOther == UINT_MAX)
			{
			fprintf(g_ftsv, "%s", Q.m_Label.c_str());
			fprintf(g_ftsv, "\tnot_rdrp.nohit");
			fprintf(g_ftsv, "\t.");
			fprintf(g_ftsv, "\t.");
			fprintf(g_ftsv, "\t.");
			fprintf(g_ftsv, "\t.");
			fprintf(g_ftsv, "\n");
			}
		else if (RefIndex_TopRdRp == UINT_MAX && RefIndex_TopOther != UINT_MAX)
			{
			fprintf(g_ftsv, "%s", Q.m_Label.c_str());
			fprintf(g_ftsv, "\tnot_rdrp.decoyhit");
			fprintf(g_ftsv, "\t.");
			fprintf(g_ftsv, "\t.");
			fprintf(g_ftsv, "\t%s", Rs[RefIndex_TopOther]->m_Label.c_str());
			fprintf(g_ftsv, "\t%.2f", RMSD_TopOther);
			fprintf(g_ftsv, "\n");
			}
		else if (RefIndex_TopRdRp != UINT_MAX && RefIndex_TopOther == UINT_MAX)
			{
			++HitCount;
			fprintf(g_ftsv, "%s", Q.m_Label.c_str());
			fprintf(g_ftsv, "\trdrp.nodecoyhit");
			fprintf(g_ftsv, "\t%s", Rs[RefIndex_TopRdRp]->m_Label.c_str());
			fprintf(g_ftsv, "\t%.2f", RMSD_TopRdRp);
			fprintf(g_ftsv, "\t.");
			fprintf(g_ftsv, "\t.");
			fprintf(g_ftsv, "\n");
			}
		else if (RefIndex_TopRdRp != UINT_MAX && RefIndex_TopOther != UINT_MAX)
			{
			string Type;
			if (RMSD_TopRdRp <= RMSD_TopOther)
				{
				Type = "rdrp.tophit";
				++HitCount;
				}
			else
				{
				if (RMSD_TopRdRp - RMSD_TopOther < 0.5)
					Type = "rdrp-like";
				else
					Type = "not_rdrp.tophit";
				}

			fprintf(g_ftsv, "%s", Q.m_Label.c_str());
			fprintf(g_ftsv, "\t%s", Type.c_str());
			fprintf(g_ftsv, "\t%s", Rs[RefIndex_TopRdRp]->m_Label.c_str());
			fprintf(g_ftsv, "\t%.2f", RMSD_TopRdRp);
			fprintf(g_ftsv, "\t%s", Rs[RefIndex_TopOther]->m_Label.c_str());
			fprintf(g_ftsv, "\t%.2f", RMSD_TopOther);
			fprintf(g_ftsv, "\n");
			}
		else
			asserta(false);
		}
	Progress("100%% done, %u / %u RdRp hits\n",
		HitCount, RecCount);
	}

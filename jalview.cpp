#include "myutils.h"
#include "outputfiles.h"
#include "shapesearcher.h"

void ShapeSearcher::ToJalview(FILE *f) const
	{
	if (f == 0)
		return;
	static bool HdrDone = false;
#pragma omp critical
	{
	if (!HdrDone)
		{
		HdrDone = true;
		fprintf(f, "A	0000ff\n");
		fprintf(f, "B	00ff00\n");
		fprintf(f, "C	ff0000\n");
		fprintf(f, "D	ffff00\n");
		fprintf(f, "E	ffa000\n");
		fprintf(f, "F1	1e90ff\n");
		fprintf(f, "F2	00bfff\n");
		fprintf(f, "H	dcdcdc\n");
		fprintf(f, "J	ff00ff\n");
		fprintf(f, "STARTGROUP	Motifs\n");
		}

	string FoundMotifsStr;
	GetFoundMotifsStr(FoundMotifsStr);

	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		uint Pos = m_ShapePosVec[i];
		if (Pos == UINT_MAX)
			continue;

		string Seq;
		GetShapeSeq(i, Seq);

		fprintf(f, "-");
		fprintf(f, "\t%s", m_Query->m_Label.c_str());
		fprintf(f, "\t-1");
		fprintf(f, "\t%u", Pos+1);
		fprintf(f, "\t%u", Pos+SIZE(Seq));
		fprintf(f, "\t%s\n", GetShapeName(i));
		}
	fprintf(f, "\n");
	}
	}

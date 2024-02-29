#include "myutils.h"
#include "chainreader.h"
#include "shapesearcher.h"

void cmd_shapes_score1()
	{
	const string &QueryFN = opt_shapes_score1;

	if (!optset_shapes)
		Die("Must specify -shapes");
	const string &ShapesFileName = opt_shapes;

	asserta(optset_shapename);

	uint Lo = UINT_MAX;
	uint Hi = UINT_MAX;
	if (optset_pos)
		{
		Lo = opt_pos;
		Hi = opt_pos;
		}
	else if (optset_lo && optset_hi)
		{
		Lo = opt_lo;
		Hi = opt_hi;
		}

	Shapes S;
	S.FromFile(ShapesFileName);

	uint ShapeIndex = S.GetShapeIndex(opt_shapename);
	double MinScore = 0.5;
	if (optset_minscore)
		MinScore = opt_minscore;

	ShapeSearcher SS;
	SS.Init(S);

	vector<PDBChain *> Chains;
	ReadChains(QueryFN, Chains);

	const uint ChainCount = SIZE(Chains);
	asserta(ChainCount > 0);
	if (ChainCount > 1)
		Warning(">1 chains, first only");
	const PDBChain &Query = *Chains[0];

	SS.SetQuery(Query);

	vector<uint> HitPosVec;
	vector<double> HitScores;
	SS.SearchShapeSelf(ShapeIndex, 0.1, Lo, Hi, 0, UINT_MAX,
	  HitPosVec, HitScores);

	uint HitCount = SIZE(HitPosVec);

	Log("%u hits\n", HitCount);
	}

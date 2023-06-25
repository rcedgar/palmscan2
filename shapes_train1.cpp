#include "myutils.h"
#include "shapesearcher.h"
#include "rdrpsearcher.h"
#include <map>

/***
shapes_train1

Three-pass training on new motif (X).
ABC annotated by PSSMs.
Motif X is annotated (A0) on a few training examples (F).
F is a subset of a larger training set (T).
	Pass 1. Train Shapes S1 on F.
	Pass 2. Search T using S1, gives new training set A2.
	Pass 3. Train S2 on A2.
	Pass 3. Search T using S2 giving A3, compare A3 vs. A2 for consistency.

Input:
	Shapes profile with ABC
	Training set T.
	tsv file with a motif annotations A0 for F.

Output:
	tsv file with ABCX for T.
***/

void cmd_shapes_train1()
	{
	vector<string> MotifSeqsVec_A0;
	}

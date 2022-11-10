#include "myutils.h"
#include "trisearcher.h"

void InitTS(TriSearcher& TS)
	{
	TS.Radius = 1.5;
	TS.MaxTriRMSD2 = 1.5;
	TS.MaxMotifRMSD2 = 5.5;
	TS.NABmin = 10;
	TS.NABmax = 80;
	TS.NBCmin = 10;
	TS.NBCmax = 80;
	TS.NACmin = 80;
	TS.NACmax = 200;
	}

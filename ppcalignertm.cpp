#include "myutils.h"
#include "ppcaligner.h"
#include "abcxyz.h"

double AlignTmPpc(const PDBChain &Query, 
  const PDBChain &Ref, string &Path);

double PPCAligner::GetTMScore(const TSHit &Hit, string &Path) const
	{
	const PDBChain &Q = *Hit.m_Query;
	const PDBChain &R = *Hit.m_Ref;
	double TM = AlignTmPpc(Q, R, Path);
	return TM;
	}

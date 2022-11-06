#include "myutils.h"
#include "rdrpmodel.h"
#include "xdpmem.h"
#include "objmgr.h"
#include "alnparams.h"
#include "pathinfo.h"
#include "sort.h"
#include <set>

void SeqToFasta(FILE *f, const char *Label, const char *Seq, unsigned L);

// PSSMs must be in current directory, filenames GroupName.X where X=A,B,C.
void RdRpModel::FromPSSMs(const string &PSSMDir,
  const vector<string> &GroupNames)
	{
	Clear();

	m_GroupNames = GroupNames;
	const uint GroupCount = SIZE(m_GroupNames);
	m_PSSMs.resize(GroupCount*3);
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		string FileNameA = PSSMDir + m_GroupNames[GroupIndex] + ".A";
		string FileNameB = PSSMDir + m_GroupNames[GroupIndex] + ".B";
		string FileNameC = PSSMDir + m_GroupNames[GroupIndex] + ".C";

		m_PSSMs[3*GroupIndex + 0].FromFile(FileNameA);
		m_PSSMs[3*GroupIndex + 1].FromFile(FileNameB);
		m_PSSMs[3*GroupIndex + 2].FromFile(FileNameC);
		}
	}


uint RdRpModel::GetPSSMLength(uint GroupIndex, uint MotifIndex) const
	{
	const PSSM &P = GetPSSM(GroupIndex, MotifIndex);
	uint L = P.GetColCount();
	return L;
	}

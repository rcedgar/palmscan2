#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include <set>

void ReadLinesFromFile(const string &FileName, vector<string> &Lines);

static bool MatchLabelSubstr(const string &QueryLabel,
  const string &RefLabel)
	{
	return (QueryLabel.find(RefLabel) != string::npos);
	}

static bool MatchLabelSet(const string &QueryLabel,
  const set<string> &RefLabels, bool MatchSubstr)
	{
	if (!MatchSubstr)
		{
		bool Found = (RefLabels.find(QueryLabel) != RefLabels.end());
		return Found;
		}
	for (set<string>::const_iterator p = RefLabels.begin();
	  p != RefLabels.end(); ++p)
		{
		bool Found = MatchLabelSubstr(QueryLabel, *p);
		if (Found)
			return true;
		}
	return false;
	}

void cmd_getchains()
	{
	const string &InputFileName = opt_getchains;

	asserta(opt_labels != "");

	bool MatchSubstr = false;
	if (optset_label_substr_match)
		MatchSubstr = true;

	vector<string> Labels;
	ReadLinesFromFile(opt_labels, Labels);
	set<string> LabelSet;
	const uint LabelCount = SIZE(Labels);
	vector<string> Fields;
	for (uint i = 0; i < LabelCount; ++i)
		{
		const string &Label = Labels[i];
		LabelSet.insert(Label);
		}

	vector<PDBChain *> Chains;
	ReadChains(InputFileName, Chains);
	const uint N = SIZE(Chains);

	uint n = 0;
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N+1, "Searching %u / %u found", n, LabelCount);

		PDBChain &Q = *Chains[i];
		const string &Label = Q.m_Label;
		if (MatchLabelSet(Label, LabelSet, MatchSubstr))
			{
			Q.ToCal(g_fcal);
			++n;
			}
		}
	ProgressStep(N, N+1, "Searching %u / %u found", n, LabelCount);
	}

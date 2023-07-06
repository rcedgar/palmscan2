#include "myutils.h"
#include "shapesearcher.h"

void ShapeSearcher::GetFirstIndexes(vector<uint> &Indexes) const
	{
	Indexes.clear();
	asserta(SIZE(m_IncludeShapes) == m_ShapeCount);
	asserta(SIZE(m_SelfTopHits) == m_ShapeCount);
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		if (m_SelfTopHits[i].size() == 0)
			Indexes.push_back(UINT_MAX);
		else
			Indexes.push_back(0);
		}
	}

bool ShapeSearcher::GetNextIndexes(vector<uint> &Indexes) const
	{
	asserta(SIZE(Indexes) == m_ShapeCount);
	asserta(SIZE(m_IncludeShapes) == m_ShapeCount);
	asserta(SIZE(m_SelfTopHits) == m_ShapeCount);
	bool Ok = false;
	for (uint i = 0; i < m_ShapeCount; ++i)
		{
		uint n = m_SelfTopHits[i].size();
		if (n == 0)
			continue;
		uint CurrentIndex = Indexes[i];
		if (CurrentIndex + 1 < n)
			{
			Indexes[i] = CurrentIndex + 1;
			Ok = true;
			continue;
			}
		}
	return Ok;
	}

void ShapeSearcher::Search(const vector<bool> &IncludeShapes)
	{
	m_IncludeShapes = IncludeShapes;
	asserta(m_ShapeCount > 0);
	asserta(SIZE(m_IncludeShapes) == m_ShapeCount);
	m_SelfTopHits.clear();
	m_SelfTopHits.resize(m_ShapeCount);

	uint QL = GetQL();
	if (QL < 32)
		return;
	vector<double> HitScores;
	uint Total = 0;
	for (uint ShapeIndex = 0; ShapeIndex < m_ShapeCount; ++ShapeIndex)
		{
		SearchShapeSelfTop(ShapeIndex, m_MinSelfScore, m_MaxTopHitCount,
		  m_SelfTopHits[ShapeIndex]);
		uint n = SIZE(m_SelfTopHits[ShapeIndex]);
		Total += n;
		}
	if (Total == 0)
		return;

	vector<uint> Indexes;
	GetFirstIndexes(Indexes);
	for (;;)
		{
		bool Ok = GetNextIndexes(Indexes);
		if (!Ok)
			break;
		}
	}

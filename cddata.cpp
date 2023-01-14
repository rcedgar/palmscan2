#include "myutils.h"
#include "cddata.h"

void CDData::Init(const CDInfo &Info)
	{
	m_Info = &Info;
	uint Size = m_Info->GetSize();
	m_Data.clear();
	m_Data.resize(Size);
	for (uint i = 0; i < Size; ++i)
		m_Data[i].resize(Size, DBL_MAX);
	}

double CDData::Get(uint MotifIndex, uint i) const
	{
	uint Ix = m_Info->GetIx(MotifIndex, i);
	assert(Ix < SIZE(m_Data));
	return m_Data[Ix];
	}

void CDData::Set(uint MotifIndex, uint i, double Value)
	{
	uint Ix = m_Info->GetIx(MotifIndex, i);
	assert(Ix < SIZE(m_Data));
	m_Data[Ix] = Value;
	}

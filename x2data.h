#pragma once

class X2Data
	{
public:
// Uniques
	vector<char> m_Aminos;
	vector<vector<double> > m_FeatureValuesVec;
	vector<string> m_Doms;
	vector<uint> m_Coords;

// Aligned pairs
	vector<char> m_Aminos1;
	vector<char> m_Aminos2;
	vector<vector<double> > m_FeatureValuesVec1;
	vector<vector<double> > m_FeatureValuesVec2;

public:
	void FromTsv(const string &FileName);
	};

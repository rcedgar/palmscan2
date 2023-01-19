#include "myutils.h"
#include "xprof.h"
#include "xprofdata.h"

bool XProfData::FromCfv(FILE *f)
	{
	Clear();

	asserta(f != 0);
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	if (!Ok)
		return false;
	Split(Line, Fields, '\t');
	uint NF = SIZE(Fields);
	asserta(NF > 5);
	m_Label = Fields[1];
	uint L = StrToUint(Fields[2]);
	asserta(L > 0);
	uint FeatureCount = StrToUint(Fields[3]);
	asserta(FeatureCount == XProf::GetFeatureCount());
	asserta(FeatureCount + 4 == NF);
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount;
	  ++FeatureIndex)
		{
		string Name = string(XProf::GetFeatureName(FeatureIndex));
		asserta(Fields[4+FeatureIndex] == Name);
		}

	m_PosToFeatureVec.clear();
	m_PosToIntFeatureVec.clear();

	vector<double> FeatureVec(FeatureCount);
	vector<uint> IntFeatureVec(FeatureCount);
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2*FeatureCount + 5);
		asserta(StrToUint(Fields[0]) == Pos+1);
		asserta(SIZE(Fields[1]) == 1);
		m_Seq += Fields[1][0];
		double x = StrToFloat(Fields[2]);
		double y = StrToFloat(Fields[3]);
		double z = StrToFloat(Fields[4]);

		m_Xs.push_back(x);
		m_Ys.push_back(y);
		m_Zs.push_back(z);

		for (uint FeatureIndex = 0; FeatureIndex < FeatureCount;
		  ++FeatureIndex)
			{
			const string &Floatf = Fields[2*FeatureIndex + 5];
			const string &Intf = Fields[2*FeatureIndex + 6];
			double Value = (Floatf == "." ? DBL_MAX : StrToFloat(Floatf));
			uint iValue = (Intf == "." ? UINT_MAX : StrToUint(Intf));
			FeatureVec[FeatureIndex] = Value;
			IntFeatureVec[FeatureIndex] = iValue;
			}
		m_PosToFeatureVec.push_back(FeatureVec);
		m_PosToIntFeatureVec.push_back(IntFeatureVec);
		}
	return true;
	}

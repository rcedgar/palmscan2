#include "myutils.h"
#include "pdb.h"
#include "abcxyz.h"

void GetThreeFromOne(char aa, string &AAA);

/***
PDB ATOM record format
http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.

MODEL     1                                                                     

          1         2         3         4         5         6         7         
01234567890123456789012345678901234567890123456789012345678901234567890123456789
                              xxxxxxxxyyyyyyyyzzzzzzzz
ATOM      1  N   PHE A   1      34.582  19.022  -8.646  1.00 24.35           N  
ATOM      2  CA  PHE A   1      33.319  19.558  -8.153  1.00 24.35           C  
ATOM      3  C   PHE A   1      32.243  19.483  -9.229  1.00 24.35           C  
ATOM      4  CB  PHE A   1      33.494  21.006  -7.685  1.00 24.35           C  
***/

void PDB::LogMe(bool WithCoords) const
	{
	const size_t L = m_Xs.size();
	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	Log("\n");
	Log(">%s\n", m_Label.c_str());
	Log("%s\n", m_Seq.c_str());
	if (!m_MotifPosVec.empty())
		{
		asserta(m_MotifPosVec.size() == 3);
		string A, B, C;
		Log("A:%u:%s  B:%u:%s  C:%u:%s\n",
		  m_MotifPosVec[0]+1, GetMotifSeq(0, A),
		  m_MotifPosVec[1]+1, GetMotifSeq(1, B),
		  m_MotifPosVec[2]+1, GetMotifSeq(2, C));
		}

	if (WithCoords)
		{
		Log("\n");
		Log(" Pos  S         X         Y         Z\n");
		for (size_t i = 0; i < L; ++i)
			{
			Log("%4u", i);
			Log("  %c", m_Seq[i]);
			Log("  %8.3f", m_Xs[i]);
			Log("  %8.3f", m_Ys[i]);
			Log("  %8.3f", m_Zs[i]);
			Log("\n");
			}
		}
	}

void PDB::ToCal(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	ToCal(f);
	CloseStdioFile(f);
	}

void PDB::ToCal(FILE *f) const
	{
	const size_t L = m_Xs.size();
	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	fprintf(f, ">%s", m_Label.c_str());
	if (!m_MotifPosVec.empty())
		{
		asserta(m_MotifPosVec.size() == 3);
		string A, B, C;
		fprintf(f, " A:%u:%s B:%u:%s C:%u:%s",
		  m_MotifPosVec[0] + 1, GetMotifSeq(0, A), 
		  m_MotifPosVec[1] + 1, GetMotifSeq(1, B), 
		  m_MotifPosVec[2] + 1, GetMotifSeq(2, C));
		}
	fputc('\n', f);

	for (size_t i = 0; i < L; ++i)
		{
		fputc(m_Seq[i], f);
		fputc('\t', f);
		fprintf(f, "\t%.3f", m_Xs[i]);
		fprintf(f, "\t%.3f", m_Ys[i]);
		fprintf(f, "\t%.3f", m_Zs[i]);
		fputc('\n', f);
		}
	}

void PDB::ToPDB(const string &FileName) const
	{
	if (FileName == "")
		return;

	FILE *f = CreateStdioFile(FileName);

	const size_t L = m_Xs.size();
	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	fprintf(f, "TITLE %s\n", m_Label.c_str());
	if (!m_MotifPosVec.empty())
		{
		asserta(m_MotifPosVec.size() == 3);
		string A, B, C;
		fprintf(f, "REMARK PALMPRINT A:%u:%s B:%u:%s C:%u:%s\n",
		  m_MotifPosVec[0] + 1, GetMotifSeq(0, A), 
		  m_MotifPosVec[1] + 1, GetMotifSeq(1, B), 
		  m_MotifPosVec[2] + 1, GetMotifSeq(2, C));

		string SSA, SSB, SSC;
		GetSS(m_MotifPosVec[0], AL, SSA);
		GetSS(m_MotifPosVec[1], BL, SSB);
		GetSS(m_MotifPosVec[2], CL, SSC);

		fprintf(f, "REMARK PPMOTIFSS A:%u:%s B:%u:%s C:%u:%s\n",
		  m_MotifPosVec[0] + 1, SSA.c_str(),
		  m_MotifPosVec[1] + 1, SSB.c_str(),
		  m_MotifPosVec[2] + 1, SSC.c_str());
		}

	for (uint i = 0; i < L; ++i)
		{
		char aa = m_Seq[i];
		string sAAA;
		GetThreeFromOne(aa, sAAA);
		const char *AAA = sAAA.c_str();

		fprintf(f, "ATOM  ");			//  1 -  6        Record name   "ATOM  "
		fprintf(f, "%5u", i+1);			//  7 - 11        Integer       serial       Atom  serial number.
		fprintf(f, " ");				// 12
		fprintf(f, " CA ");				// 13 - 16        Atom          name         Atom name.
		fprintf(f, " ");				// 17             Character     altLoc       Alternate location indicator.
		fprintf(f, "%3.3s", AAA);		// 18 - 20        Residue name  resName      Residue name.
		fprintf(f, " ");				// 21
		fprintf(f, "A");				// 22             Character     chainID      Chain identifier.
		fprintf(f, "%4u", i+1);			// 23 - 26        Integer       resSeq       Residue sequence number.
		fprintf(f, " ");				// 27             AChar         iCode        Code for insertion of residues.
		fprintf(f, "   ");				// 28 - 30
		fprintf(f, "%8.3f", m_Xs[i]);	// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		fprintf(f, "%8.3f", m_Ys[i]);	// 39 - 46        Real(8.3)     y            Orthogonal coordinates for X in Angstroms.
		fprintf(f, "%8.3f", m_Zs[i]);	// 47 - 57        Real(8.3)     z            Orthogonal coordinates for X in Angstroms.
		fprintf(f, "%6.2f", 1.0);		// 55 - 60        Real(6.2)     occupancy    Occupancy.
		fprintf(f, "%6.2f", 0.0);		// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		fprintf(f, "          ");		// 67 - 76
		fprintf(f, " C");				// 77 - 78        LString(2)    element      Element symbol, right-justified.
		fprintf(f, "  ");				// 79 - 80        LString(2)    charge       Charge  on the atom.

		fprintf(f, "\n");
		}
	}

void PDB::FromFile(const string &FileName)
	{
	Clear();
	m_Label = FileName;

	FILE *f = OpenStdioFile(FileName);

	string Line;
	while(ReadLineStdioFile(f, Line))
		{
		const size_t L = Line.size();
		if (L < 50)
			continue;

		if (Line.substr(0, 16) == "REMARK PALMPRINT")
			{
		//                  0 1              2  3                4   5        6
		// REMARK PALMPRINT A:1:VATDVSDHDTFW B:73:SGQGATDLMGTLLM C:130:SKSDDAML
			vector<string> Fields;
			Split(Line, Fields, ':');
			if (SIZE(Fields) != 7)
				continue;
			uint PosA = StrToUint(Fields[1]);
			uint PosB = StrToUint(Fields[3]);
			uint PosC = StrToUint(Fields[5]);

			asserta(PosA > 0);
			asserta(PosB > 0);
			asserta(PosC > 0);

			m_MotifPosVec.clear();
			m_MotifPosVec.push_back(PosA-1);
			m_MotifPosVec.push_back(PosB-1);
			m_MotifPosVec.push_back(PosC-1);
			}

		if (Line.substr(0, 6) != "ATOM  ")
			continue;

		string AtomName = Line.substr(12, 4);
		StripWhiteSpace(AtomName);
		if (AtomName != "CA")
			continue;

		string AAA = Line.substr(17, 3);
		char aa = GetOneFromThree(AAA);
		m_Seq.push_back(aa);

		string sX, sY, sZ;
		sX = Line.substr(30, 8);
		sY = Line.substr(38, 8);
		sZ = Line.substr(46, 8);

		StripWhiteSpace(sX);
		StripWhiteSpace(sY);
		StripWhiteSpace(sZ);

		double X = StrToFloat(sX);
		double Y = StrToFloat(sY);
		double Z = StrToFloat(sZ);

		m_Xs.push_back(X);
		m_Ys.push_back(Y);
		m_Zs.push_back(Z);
		}

	CloseStdioFile(f);
	}

const char *PDB::GetMotifSeq(uint i, string &s) const
	{
	asserta(m_MotifPosVec.size() == 3);

	s.clear();
	uint Pos = m_MotifPosVec[i];
	size_t L = m_Seq.size();
	uint ML = 0;
	switch (i)
		{
	case 0:		ML = AL; break;
	case 1:		ML = BL; break;
	case 2:		ML = CL; break;
	default:	asserta(false);
		}
	asserta(Pos + ML <= L);
	for (uint k = 0; k < ML; ++k)
		{
		char c = m_Seq[Pos+k];
		s += c;
		}
	return s.c_str();
	}

void PDB::GetPt(uint Pos, vector<double> &Pt) const
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));

	Resize3(Pt);
	Pt[X] = m_Xs[Pos];
	Pt[Y] = m_Ys[Pos];
	Pt[Z] = m_Zs[Pos];
	}

void PDB::SetPt(uint Pos, const vector<double> &Pt)
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));

	m_Xs[Pos] = Pt[X];
	m_Ys[Pos] = Pt[Y];
	m_Zs[Pos] = Pt[Z];
	}

void PDB::GetXYZ(uint Pos, double &x, double &y, double &z) const
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));
	x = m_Xs[Pos];
	y = m_Ys[Pos];
	z = m_Zs[Pos];
	}

double PDB::GetDist(uint Pos1, uint Pos2) const
	{
	double x1, y1, z1;
	double x2, y2, z2;
	GetXYZ(Pos1, x1, y1, z1);
	GetXYZ(Pos2, x2, y2, z2);
	double d = GetDist3D(x1, y1, z1, x2, y2, z2);
	return d;
	}

double PDB::GetDist2(uint Pos1, uint Pos2) const
	{
	double x1 = m_Xs[Pos1];
	double y1 = m_Ys[Pos1];
	double z1 = m_Zs[Pos1];

	double x2 = m_Xs[Pos2];
	double y2 = m_Ys[Pos2];
	double z2 = m_Zs[Pos2];

	double dx = x1 - x2;
	double dy = y1 - y2;
	double dz = z1 - z2;

	double d = GetDist(Pos1, Pos2);

	double d2 = dx*dx + dy*dy + dz*dz;
	asserta(feq(d*d, d2));
	return d2;
	}

void PDB::GetMotifDists(double &AB, double &BC, double &AC) const
	{
	vector<vector<double> > MotifCoords;
	GetMotifCoords(MotifCoords);
	AB = GetDist_Mxij(MotifCoords, A, B);
	BC = GetDist_Mxij(MotifCoords, B, C);
	AC = GetDist_Mxij(MotifCoords, A, C);
	}

void PDB::GetMotifDists2(double &AB, double &BC, double &AC) const
	{
	vector<vector<double> > MotifCoords;
	GetMotifCoords(MotifCoords);
	AB = GetDist2_Mxij(MotifCoords, A, B);
	BC = GetDist2_Mxij(MotifCoords, B, C);
	AC = GetDist2_Mxij(MotifCoords, A, C);
	}

void PDB::GetMotifCoords(vector<vector<double> > &MotifCoords) const
	{
	asserta(m_MotifPosVec.size() == 3);

	Resize3x3(MotifCoords);

	uint PosA = m_MotifPosVec[0];
	uint PosB = m_MotifPosVec[1];
	uint PosC = m_MotifPosVec[2];

	MotifCoords[A][X] = m_Xs[PosA];
	MotifCoords[A][Y] = m_Ys[PosA];
	MotifCoords[A][Z] = m_Zs[PosA];

	MotifCoords[B][X] = m_Xs[PosB];
	MotifCoords[B][Y] = m_Ys[PosB];
	MotifCoords[B][Z] = m_Zs[PosB];

	MotifCoords[C][X] = m_Xs[PosC];
	MotifCoords[C][Y] = m_Ys[PosC];
	MotifCoords[C][Z] = m_Zs[PosC];
	}

uint PDB::GetSeqLength() const
	{
	return SIZE(m_Seq);
	}

void PDB::GetMotifSeq(uint MotifStartPos, uint MotifLength,
  bool FailOnOverflow, string &MotifSeq) const
	{
	MotifSeq.clear();
	const uint L = GetSeqLength();
	for (uint i = 0; i < MotifLength; ++i)
		{
		uint Pos = MotifStartPos + i;
		if (Pos < 0 || Pos >= L)
			{
			if (FailOnOverflow)
				Die("GetMotifSeqFromMidPos overflow");
			MotifSeq += '!';
			}
		else
			MotifSeq += m_Seq[Pos];
		}
	}

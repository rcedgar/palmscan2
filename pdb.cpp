#include "myutils.h"
#include "pdb.h"
#include "abcxyz.h"

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

void PDB::LogMe() const
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
		Log("A %u=%s, B %u=%s, C %u=%s\n",
		  m_MotifPosVec[0], GetMotifSeq(0, A),
		  m_MotifPosVec[1], GetMotifSeq(1, B),
		  m_MotifPosVec[2], GetMotifSeq(2, C));
		}
	Log("\n");
	Log(" Pos  S  xxxxxxxx  yyyyyyyy  zzzzzzzz\n");
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
	uint ML = (i == 2 ? 8 : 12);
	asserta(Pos >= ML/2);
	uint Lo = Pos - ML/2;
	asserta(Lo + ML <= L);
	for (uint k = 0; k < ML; ++k)
		{
		char c = m_Seq[Lo+k];
		s += c;
		}
	return s.c_str();
	}

void PDB::GetXYZ(uint Pos, double &x, double &y, double &z) const
	{
	asserta(Pos < SIZE(m_Xs));
	asserta(Pos < SIZE(m_Ys));
	asserta(Pos < SIZE(m_Zs));
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

void PDB::GetMotifTriangle(double &LAB, double &LBC, double &LAC) const
	{
	uint PosA = m_MotifPosVec[0];
	uint PosB = m_MotifPosVec[1];
	uint PosC = m_MotifPosVec[2];

	LAB = GetDist(PosA, PosB);
	LBC = GetDist(PosB, PosC);
	LAC = GetDist(PosA, PosC);
	}

uint PDB::GetSeqLength() const
	{
	return SIZE(m_Seq);
	}

void PDB::GetMotifSeqFromMidPos(uint MidPos, uint MotifLength,
  bool FailOnOverflow, string &MotifSeq) const
	{
	MotifSeq.clear();
	const int L = (int) GetSeqLength();
	const int Len2 = int(MotifLength)/2;
	for (int i = 0; i < (int) MotifLength; ++i)
		{
		int Pos = (int) MidPos - Len2 + i;
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

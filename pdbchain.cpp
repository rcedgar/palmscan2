#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"

void GetThreeFromOne(char aa, string &AAA);

void ReadLinesFromFile(const string &FileName, vector<string> &Lines)
	{
	Lines.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	while (ReadLineStdioFile(f, Line))
		{
		Lines.push_back(Line);
		}
	CloseStdioFile(f);
	}

/***
PDBChain ATOM record format
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

void PDBChain::LogMe(bool WithCoords) const
	{
	const size_t L = m_Xs.size();
	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	string Label;
	GetLabel(Label);

	Log("\n");
	Log(">%s\n", Label.c_str());
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

void PDBChain::ToCal(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	ToCal(f);
	CloseStdioFile(f);
	}

void PDBChain::ToCal(FILE* f) const
	{
	uint QL = GetSeqLength();
	ToCalSeg(f, 0, QL);
	}

void PDBChain::ToCalSeg(FILE *f, uint Pos, uint n) const
	{
	const size_t L = m_Xs.size();
	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	string Label;
	GetLabel(Label);

	fprintf(f, ">%s", Label.c_str());
	if (!m_MotifPosVec.empty())
		{
		asserta(m_MotifPosVec.size() == 3);

		uint PosA = m_MotifPosVec[0];
		uint PosB = m_MotifPosVec[1];
		uint PosC = m_MotifPosVec[2];
		asserta(PosA >= Pos);
		asserta(PosC + CL <= Pos + n);

		string A, B, C;
		fprintf(f, " A:%u:%s B:%u:%s C:%u:%s",
		  PosA - Pos + 1, GetMotifSeqNoFail(0, A), 
		  PosB - Pos + 1, GetMotifSeqNoFail(1, B), 
		  PosC - Pos + 1, GetMotifSeqNoFail(2, C));
		}
	fputc('\n', f);

	for (size_t i = Pos; i < Pos + n; ++i)
		{
		if (i < SIZE(m_Seq))
			fputc(m_Seq[i], f);
		else
			fputc('*', f);
		fprintf(f, "\t%.3f", m_Xs[i]);
		fprintf(f, "\t%.3f", m_Ys[i]);
		fprintf(f, "\t%.3f", m_Zs[i]);
		fputc('\n', f);
		}
	}

void PDBChain::ToPDB(const string &FileName) const
	{
	if (FileName == "")
		return;

	FILE *f = CreateStdioFile(FileName);

	const size_t L = m_Xs.size();
	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	char Chain = m_Chain;
	if (Chain == 0)
		Chain = 'A';

	fprintf(f, "TITLE %s\n", m_ChainLabel.c_str());
	if (!m_MotifPosVec.empty())
		{
		asserta(m_MotifPosVec.size() == 3);
		string A, B, C;
		fprintf(f, "REMARK PALMPRINT A:%u:%s B:%u:%s C:%u:%s\n",
		  m_MotifPosVec[0] + 1, GetMotifSeq(0, A), 
		  m_MotifPosVec[1] + 1, GetMotifSeq(1, B), 
		  m_MotifPosVec[2] + 1, GetMotifSeq(2, C));
		}

	for (uint i = 0; i < L; ++i)
		{
		char aa = m_Seq[i];
		string sAAA;
		GetThreeFromOne(aa, sAAA);
		const char *AAA = sAAA.c_str();

		fprintf(f, "ATOM  ");			//  1 -  6        Record name   "ATOM  "
		fprintf(f, "%5u", i+1);			//  7 - 11        Integer       serial       Atom serial number.
		fprintf(f, " ");				// 12
		fprintf(f, " CA ");				// 13 - 16        Atom          name         Atom name.
		fprintf(f, " ");				// 17             Character     altLoc       Alternate location indicator.
		fprintf(f, "%3.3s", AAA);		// 18 - 20        Residue name  resName      Residue name.
		fprintf(f, " ");				// 21
		fprintf(f, "%c", Chain);		// 22             Character     chainID      Chain identifier.
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
		fprintf(f, "  ");				// 79 - 80        LString(2)    charge       Charge on the atom.

		fprintf(f, "\n");
		}
	}

void PDBChain::FromLines(const string &Label, char ChainChar,
  const vector<string> &Lines)
	{
	Clear();

	vector<string> Fields;
	Split(Label, Fields, '/');
	uint FieldCount = SIZE(Fields);
	asserta(FieldCount > 0);
	string FixedLabel = Fields[FieldCount-1];
	if (EndsWith(FixedLabel, ".pdb"))
		{
		uint n = SIZE(FixedLabel);
		if (n > 4)
			FixedLabel = FixedLabel.substr(0, n-4);
		}

	m_Chain = ChainChar;
	m_ChainLabel = FixedLabel;
	const uint N = SIZE(Lines);
	m_Chain = 0;
	for (uint LineNr = 0; LineNr < N; ++LineNr)
		{
		const string &Line = Lines[LineNr];
		const size_t L = Line.size();

		if (strncmp(Line.c_str(), "REMARK PALMPRINT", 16) == 0)
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

		if (strncmp(Line.c_str(), "ATOM  ", 6) != 0)
			continue;

		string AtomName = Line.substr(12, 4);
		StripWhiteSpace(AtomName);
		if (AtomName != "CA")
			continue;

		char Chain = Line[21];
		if (m_Chain == 0)
			m_Chain = Chain;
		else if (Chain != m_Chain)
			Die("PDBChain::FromLines() two chains %c, %c", m_Chain, Chain);
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
	}

const char *PDBChain::GetMotifSeqNoFail(uint i, string &s) const
	{
	if (m_MotifPosVec.size() != 3)
		{
		s = "ERR1";
		return s.c_str();
		}

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
	for (uint k = 0; k < ML; ++k)
		{
		if (Pos+k < SIZE(m_Seq))
			{
			char c = m_Seq[Pos+k];
			s += c;
			}
		else
			s += "*";
		}
	return s.c_str();
	}

const char *PDBChain::GetMotifSeq(uint i, string &s) const
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

void PDBChain::GetPt(uint Pos, vector<double> &Pt) const
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));

	Resize3(Pt);
	Pt[X] = m_Xs[Pos];
	Pt[Y] = m_Ys[Pos];
	Pt[Z] = m_Zs[Pos];
	}

void PDBChain::SetPt(uint Pos, const vector<double> &Pt)
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));

	m_Xs[Pos] = Pt[X];
	m_Ys[Pos] = Pt[Y];
	m_Zs[Pos] = Pt[Z];
	}

void PDBChain::GetXYZ(uint Pos, double &x, double &y, double &z) const
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));
	x = m_Xs[Pos];
	y = m_Ys[Pos];
	z = m_Zs[Pos];
	}

double PDBChain::GetDist(uint Pos1, uint Pos2) const
	{
	double x1, y1, z1;
	double x2, y2, z2;
	GetXYZ(Pos1, x1, y1, z1);
	GetXYZ(Pos2, x2, y2, z2);
	double d = GetDist3D(x1, y1, z1, x2, y2, z2);
	return d;
	}

double PDBChain::GetDist2(uint Pos1, uint Pos2) const
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

void PDBChain::GetMotifDists(double &AB, double &BC, double &AC) const
	{
	vector<vector<double> > MotifCoords;
	GetMotifCoords(MotifCoords);
	AB = GetDist_Mxij(MotifCoords, A, B);
	BC = GetDist_Mxij(MotifCoords, B, C);
	AC = GetDist_Mxij(MotifCoords, A, C);
	}

void PDBChain::GetMotifDists2(double &AB, double &BC, double &AC) const
	{
	vector<vector<double> > MotifCoords;
	GetMotifCoords(MotifCoords);
	AB = GetDist2_Mxij(MotifCoords, A, B);
	BC = GetDist2_Mxij(MotifCoords, B, C);
	AC = GetDist2_Mxij(MotifCoords, A, C);
	}

void PDBChain::GetMotifCoords(vector<vector<double> > &MotifCoords) const
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

void PDBChain::GetLabel(string &Label) const
	{
	if (m_Chain == 0 || m_Chain == 'A')
		Label = m_ChainLabel;
	else
		Label = m_ChainLabel + "_" + m_Chain;
	}

uint PDBChain::GetSeqLength() const
	{
	return SIZE(m_Seq);
	}

void PDBChain::GetSubSeq(uint MotifStartPos, uint n,
  bool FailOnOverflow, string &MotifSeq) const
	{
	MotifSeq.clear();
	const uint L = GetSeqLength();
	for (uint i = 0; i < n; ++i)
		{
		uint Pos = MotifStartPos + i;
		if (Pos < 0 || Pos >= L)
			{
			if (FailOnOverflow)
				Die("'%s' GetMotifSeqFromMidPos overflow", m_ChainLabel.c_str());
			MotifSeq += '!';
			}
		else
			MotifSeq += m_Seq[Pos];
		}
	}

void PDBChain::ReadChainsFromLines(const string &Label,
  const vector<string> &Lines, vector<PDBChain *> &Chains)
	{
	Chains.clear();
	const uint N = SIZE(Lines);
	vector<string> ChainLines;
	char CurrChainChar = 0;
	bool AnyAtoms = false;
	for (uint i = 0; i < N; ++i)
		{
		const string &Line = Lines[i];

		if (strncmp(Line.c_str(), "ATOM  ", 6) == 0)
			{
			if (Line.size() < 57)
				continue;
			char ChainChar = Line[21];
			if (ChainChar != CurrChainChar)
				{
				if (AnyAtoms && !ChainLines.empty())
					{
					PDBChain *Chain = new PDBChain;
					Chain->FromLines(Label, CurrChainChar, ChainLines);
					ChainLines.clear();
					AnyAtoms = false;

					Chains.push_back(Chain);
					}
				CurrChainChar = ChainChar;
				}
			ChainLines.push_back(Line);
			AnyAtoms = true;
			}
		else if (
		  strncmp(Line.c_str(), "REMARK PALMPRINT", 6) == 0)
			ChainLines.push_back(Line);
		}
	if (!ChainLines.empty())
		{
		PDBChain *Chain = new PDBChain;
		Chain->FromLines(Label, CurrChainChar, ChainLines);
		ChainLines.clear();
		Chains.push_back(Chain);
		}
	}

void PDBChain::ReadChainsFromFile(const string &FileName,
  vector<PDBChain *> &Chains)
	{
	vector<string> Lines;
	ReadLinesFromFile(FileName, Lines);
	ReadChainsFromLines(FileName, Lines, Chains);
	}

void PDBChain::CopyTriForm(const vector<double> &t,
  const vector<vector<double> > &R,
  PDBChain &XChain) const
	{
	XChain.Clear();
	XChain.m_ChainLabel = m_ChainLabel;
	XChain.m_Seq = m_Seq;
	XChain.m_MotifPosVec = m_MotifPosVec;
	const uint N = SIZE(m_Xs);
	asserta(SIZE(m_Seq) == N);
	asserta(SIZE(m_Ys) == N);
	asserta(SIZE(m_Zs) == N);
	vector<double> Pt(3);
	vector<double> XPt(3);
	for (uint Pos = 0; Pos < N; ++Pos)
		{
		GetPt(Pos, Pt);
		XFormPt(Pt, t, R, XPt);
		XChain.m_Xs.push_back(XPt[X]);
		XChain.m_Ys.push_back(XPt[Y]);
		XChain.m_Zs.push_back(XPt[Z]);
		}
	}

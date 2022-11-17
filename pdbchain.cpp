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

	Log("\n");
	Log(">%s\n", m_Label.c_str());
	Log("%s\n", m_Seq.c_str());
	if (!m_MotifPosVec.empty())
		{
		asserta(m_MotifPosVec.size() == 3);
		string MotifA, MotifB, MotifC;
		GetMotifSeq(A, MotifA);
		GetMotifSeq(B, MotifB);
		GetMotifSeq(C, MotifC);

		Log("A:%u:%s  B:%u:%s  C:%u:%s\n",
		  m_MotifPosVec[0]+1, MotifA.c_str(),
		  m_MotifPosVec[1]+1, MotifB.c_str(),
		  m_MotifPosVec[2]+1, MotifC.c_str());
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
	if (f == 0)
		return;
	uint QL = GetSeqLength();
	ToCalSeg(f, 0, QL);
	}

void PDBChain::ToCalSeg(FILE *f, uint Pos, uint n) const
	{
	if (f == 0)
		return;
	if (n == 0)
		return;

	const size_t L = m_Xs.size();
	if (L == 0)
		return;

	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	if (m_MotifPosVec.empty())
		fprintf(f, ">%s\n", m_Label.c_str());
	else
		{
		asserta(m_MotifPosVec.size() == 3);

		vector<string> Fields;
		Split(m_Label, Fields, ' ');
		string NewLabel = Fields[0];

		uint PosA = m_MotifPosVec[0];
		uint PosB = m_MotifPosVec[1];
		uint PosC = m_MotifPosVec[2];
		asserta(PosA >= Pos);
		asserta(PosC + CL <= Pos + n);

		string MotifA, MotifB, MotifC;
		GetMotifSeq(A, MotifA);
		GetMotifSeq(B, MotifB);
		GetMotifSeq(C, MotifC);

		fprintf(f, ">%s A:%u:%s B:%u:%s C:%u:%s\n",
		  NewLabel.c_str(),
		  m_MotifPosVec[0]+1, MotifA.c_str(),
		  m_MotifPosVec[1]+1, MotifB.c_str(),
		  m_MotifPosVec[2]+1, MotifC.c_str());
		}

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
	fprintf(f, "TITLE %s\n", m_Label.c_str());

	uint AtomCount = SIZE(m_ATOMs);
	if (AtomCount > 0)
		{
		for (uint i = 0; i < AtomCount; ++i)
			{
			const vector<string> &v = m_ATOMs[i];
			for (uint j = 0; j < SIZE(v); ++j)
				{
				fputs(v[j].c_str(), f);
				fputc('\n', f);
				}
			}
		return;
		}

	const size_t L = m_Xs.size();
	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	const char Chain = 'A';

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

char PDBChain::FromPDBLines(const string &Label,
  const vector<string> &Lines, bool SaveAtoms)
	{
	Clear();

	m_Label = Label;
	const uint N = SIZE(Lines);
	char Chain = 0;
	uint ResidueCount = 0;
	uint CurrentResidueNumber = 0;
	for (uint LineNr = 0; LineNr < N; ++LineNr)
		{
		const string &Line = Lines[LineNr];
		const size_t L = Line.size();

		if (strncmp(Line.c_str(), "ATOM  ", 6) != 0)
			continue;

		if (SaveAtoms)
			{
		// 23 - 26 residue nr.
			uint ResidueNumber = 0;
			for (uint i = 22; i < 26; ++i)
				{
				char c = Line[i];
				if (!isspace(c) && !isdigit(c))
					Die("Invalid character in residue number field (cols 23-26): %s\n",
					  Line.c_str());
				if (isdigit(c))
					ResidueNumber = ResidueNumber*10 + (c - '0');
				}
			if (ResidueNumber != CurrentResidueNumber)
				{
				CurrentResidueNumber = ResidueNumber;
				++ResidueCount;
				m_ATOMs.resize(ResidueCount);
				}
			m_ATOMs.back().push_back(Line);
			}

		string AtomName = Line.substr(12, 4);
		StripWhiteSpace(AtomName);
		if (AtomName != "CA")
			continue;

		char LineChain = Line[21];
		if (Chain == 0)
			Chain = LineChain;
		else if (Chain != LineChain)
			Die("PDBChain::FromLines() two chains %c, %c",
			  Chain, LineChain);
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
	AppendChainToLabel(m_Label, Chain);
	return Chain;
	}

void PDBChain::GetSubSeq(uint Pos, uint n, string &s) const
	{
	s.clear();
	size_t L = m_Seq.size();
	asserta(Pos + n <= L);

	for (uint i = 0; i < n; ++i)
		{
		char c = m_Seq[Pos+i];
		s += c;
		}
	}

uint PDBChain::GetMotifLength(uint i)
	{
	asserta(i < 3);
	static const uint Lengths[3] = { AL, BL, CL };
	return Lengths[i];
	}

void PDBChain::GetMotifSeq(uint i, string &Seq) const
	{
	asserta(m_MotifPosVec.size() == 3);
	uint Pos = GetMotifPos(i);
	uint ML = GetMotifLength(i);
	GetSubSeq(Pos, ML, Seq);
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

uint PDBChain::GetSeqLength() const
	{
	return SIZE(m_Seq);
	}

uint PDBChain::GetMotifPos(uint MotifIndex) const
	{
	asserta(MotifIndex < 3);
	asserta(m_MotifPosVec.size() == 3);
	uint Pos = m_MotifPosVec[MotifIndex];
	return Pos;
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
				Die("'%s' GetMotifSeqFromMidPos overflow",
				  m_Label.c_str());
			MotifSeq += '!';
			}
		else
			MotifSeq += m_Seq[Pos];
		}
	}

void PDBChain::ChainsFromLines(const string &Label,
  const vector<string> &Lines, vector<PDBChain *> &Chains,
  bool SaveAtoms)
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
					char ChainChar = Chain->FromPDBLines(Label, ChainLines, SaveAtoms);
					if (ChainChar != 0)
						Chains.push_back(Chain);
					ChainLines.clear();
					AnyAtoms = false;
					}
				CurrChainChar = ChainChar;
				}
			ChainLines.push_back(Line);
			AnyAtoms = true;
			}
		}

	if (!ChainLines.empty() && !AnyAtoms)
		{
		PDBChain *Chain = new PDBChain;
		Chain->FromPDBLines(Label, ChainLines, SaveAtoms);
		ChainLines.clear();
		Chains.push_back(Chain);
		}
	}

void PDBChain::ReadChainsFromFile(const string &FileName,
  vector<PDBChain *> &Chains, bool SaveAtoms)
	{
	vector<string> Lines;
	ReadLinesFromFile(FileName, Lines);
	ChainsFromLines(FileName, Lines, Chains, SaveAtoms);
	}

void PDBChain::CopyTriForm(const vector<double> &t,
  const vector<vector<double> > &R,
  PDBChain &XChain) const
	{
	XChain.Clear();
	XChain.m_Label = m_Label;
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

void PDBChain::AppendChainToLabel(string &Label, char Chain)
	{
	if (Chain == 0)
		return;

	string _X = "_";
	_X += Chain;
	if (!EndsWith(Label, _X))
		{
		Label += "_";
		Label += Chain;
		}
	}

uint PDBChain::GetPalmPrintLength(uint PosA, uint PosC, uint L)
	{
	asserta(PosA < PosC);
	asserta(PosC + CL <= L);

	uint Hi = PosC + CL - 1;
	uint PPL = Hi - PosA + 1;
	return PPL;
	}

bool PDBChain::CheckMotifCoords(bool FailOnError) const
	{
	uint n = SIZE(m_MotifPosVec);
	if (n != 3)
		{
		if (FailOnError)
			Die("CheckMotifCoords(), m_MotifPosVec.size()=%d", n);
		return false;
		}

	uint QL = SIZE(m_Seq);
	uint PosA = m_MotifPosVec[A];
	uint PosB = m_MotifPosVec[B];
	uint PosC = m_MotifPosVec[C];

	if (PosA + AL >= PosB)
		{
		if (FailOnError)
			Die("CheckMotifCoords(), PosA+AL>PosB");
		return false;
		}

	if (PosB + BL >= PosC)
		{
		if (FailOnError)
			Die("CheckMotifCoords(), PosB+BL>PosC");
		return false;
		}

	if (PosC + CL > QL)
		{
		if (FailOnError)
			Die("CheckMotifCoords(), PosC+CL > QC");
		return false;
		}

	return true;
	}

bool PDBChain::CheckPPCMotifCoords(bool FailOnError) const
	{
	uint n = SIZE(m_MotifPosVec);
	if (n == 0)
		{
		if (FailOnError)
			Die("CheckPPCMotifCoords(), m_MotifPosVec empty");
		return false;
		}
	asserta(n == 3);

	uint QL = SIZE(m_Seq);
	uint PosA = m_MotifPosVec[A];
	uint PosB = m_MotifPosVec[B];
	uint PosC = m_MotifPosVec[C];

	if (PosA != 0)
		{
		if (FailOnError)
			Die("CheckPPCMotifCoords(), PosA=%u", PosA);
		return false;
		}

	if (PosA + AL >= PosB)
		{
		if (FailOnError)
			Die("CheckPPCMotifCoords(), PosA+AL>PosB");
		return false;
		}

	if (PosB + BL >= PosC)
		{
		if (FailOnError)
			Die("CheckPPCMotifCoords(), PosB+BL>PosC");
		return false;
		}

	if (PosC + CL != QL)
		{
		if (FailOnError)
			{
			LogMe();
			Die("CheckPPCMotifCoords(), PosC=%u +CL != QL=%u", PosC, QL);
			}
		return false;
		}

	return true;
	}

void PDBChain::GetPPC(uint PosA, uint PosB, uint PosC,
  PDBChain &PPC) const
	{
	PPC.Clear();

	uint L = GetSeqLength();
	uint PPL = GetPalmPrintLength(PosA, PosC, L);

	asserta(PosB > PosA + AL);
	asserta(PosC > PosB + AL);
	asserta(PosC + CL <= L);

	uint PPC_PosA = 0;
	uint PPC_PosB = PosB - PosA;
	uint PPC_PosC = PosC - PosA;

	PPC.m_MotifPosVec.clear();
	PPC.m_MotifPosVec.push_back(PPC_PosA);
	PPC.m_MotifPosVec.push_back(PPC_PosB);
	PPC.m_MotifPosVec.push_back(PPC_PosC);

	vector<vector<double> > MotifCoords(3);
	GetPt(PosA, MotifCoords[A]);
	GetPt(PosB, MotifCoords[B]);
	GetPt(PosC, MotifCoords[C]);

	vector<double> t;
	vector<vector<double> > R;
	GetTriForm(MotifCoords, t, R);

	vector<double> Pt;
	vector<double> XPt;
	for (uint i = 0; i < PPL; ++i)
		{
		uint Pos = PosA + i;
		char c = m_Seq[Pos];
		PPC.m_Seq += c;

		GetPt(Pos, Pt);
		XFormPt(Pt, t, R, XPt);

		PPC.m_Xs.push_back(XPt[X]);
		PPC.m_Ys.push_back(XPt[Y]);
		PPC.m_Zs.push_back(XPt[Z]);
		}

	string MotifA, MotifB, MotifC;
	GetSubSeq(PosA, AL, MotifA);
	GetSubSeq(PosB, BL, MotifB);
	GetSubSeq(PosC, CL, MotifC);

	string MotifCoordStr;
	Ps(MotifCoordStr, "A:%u:%s B:%u:%s C:%u:%s",
	  PPC_PosA + 1, MotifA.c_str(),
	  PPC_PosB + 1, MotifB.c_str(),
	  PPC_PosC + 1, MotifC.c_str());

	PPC.m_Label = m_Label + " " + MotifCoordStr;
	PPC.CheckPPCMotifCoords(true);
	}

void PDBChain::XFormATOM(string &ATOM, const vector<double> &t,
  const vector<vector<double> > &R) const
	{
	}

bool PDBChain::IsPDBAtomLine(const string &Line)
	{
	if (strncmp(Line.c_str(), "ATOM  ", 6) == 0)
		return true;
	else
		return false;
	}

char PDBChain::GetChainCharFromPDBAtomLine(const string &Line)
	{
	if (!IsPDBAtomLine(Line))
		return 0;
	if (Line.size() < 22)
		return 0;
	return Line[21];
	}

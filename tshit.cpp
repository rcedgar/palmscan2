#include "myutils.h"
#include "tshit.h"
#include "abcxyz.h"

float GetBlosum62Score(char a, char b);
void GetPalmSketch(const string& ss, uint PSL, string& Sketch);
void GetPalmSubSketch(const string& ss, uint Lo, uint n, 
  uint PSL, string& Sketch);

const uint V_SKETCH_LENGTH = 20;

static void MakeAnnotSketch(const string &QSketch, const string &RSketch,
  string &AnnotSketch)
	{
	AnnotSketch.clear();
	const uint L = SIZE(QSketch);
	asserta(SIZE(RSketch) == L);
	for (uint i = 0; i < L; ++i)
		{
		char q = QSketch[i];
		char r = RSketch[i];
		if (q == r  && q != ' ' && r != ' ')
			AnnotSketch += '|';
		else
			AnnotSketch += ' ';
		}
	}

// D:\src\py\conf_score_test.py
double GetConfidenceScore(double MotifRMSD)
	{
	const double T = 2.0;	// RMSD where results start to degrade
	double Ratio = MotifRMSD*MotifRMSD/(T*T);
	double ex = exp(Ratio);
	double Conf = (1.0 + ex)/(2.0*ex);
	return Conf;
	}

static char GetAnnotChar(char a, char b)
	{
	if (!isalpha(a) || !isalpha(b))
		return ' ';
	if (a == b)
		return '|';

	float Score = GetBlosum62Score(a, b);

	if (Score >= 0.5)
		return '+';
	else if (Score > 0)
		return '.';
	return ' ';
	}

static void MakeAnnot(const string &A, const string &B, 
  string &Annot)
	{
	const uint L = SIZE(A);
	assert(SIZE(B) == L);
	Annot.clear();
	for (uint i = 0; i < L; ++i)
		{
		char c = GetAnnotChar(A[i], B[i]);
		Annot += c;
		}
	}

void TSHit::WriteAln(FILE *f) const
	{
	asserta(m_Query != 0 && m_Ref != 0);
	
	string QLabel;
	string RLabel;
	m_Query->GetLabel(QLabel);
	m_Ref->GetLabel(RLabel);

	uint QPosA = m_QPosA;
	uint QPosB = m_QPosB;
	uint QPosC = m_QPosC;

	string QASeq, QBSeq, QCSeq;
	m_Query->GetSubSeq(QPosA, AL, true, QASeq);
	m_Query->GetSubSeq(QPosB, BL, true, QBSeq);
	m_Query->GetSubSeq(QPosC, CL, true, QCSeq);

	asserta(m_Ref->m_MotifPosVec.size() == 3);
	uint RPosA = m_Ref->m_MotifPosVec[A];
	uint RPosB = m_Ref->m_MotifPosVec[B];
	uint RPosC = m_Ref->m_MotifPosVec[C];

	string RASeq, RBSeq, RCSeq;
	m_Ref->GetSubSeq(RPosA, AL, true, RASeq);
	m_Ref->GetSubSeq(RPosB, BL, true, RBSeq);
	m_Ref->GetSubSeq(RPosC, CL, true, RCSeq);

	string AnnotA, AnnotB, AnnotC;
	MakeAnnot(QASeq, RASeq, AnnotA);
	MakeAnnot(QBSeq, RBSeq, AnnotB);
	MakeAnnot(QCSeq, RCSeq, AnnotC);

	string QRow = QASeq + "   " + QBSeq + "   " + QCSeq;
	string RRow = RASeq + "   " + RBSeq + "   " + RCSeq;
	string AnnotRow = AnnotA + "   " + AnnotB + "   " + AnnotC;

	QRow += "  " + QLabel;
	RRow += "  " + RLabel;

	string CoordStr;
	string CoordsA;
	string CoordsB;
	string CoordsC;

	Ps(CoordsA, "A:%u", QPosA+1);
	Ps(CoordsB, "B:%u", QPosB+1);
	Ps(CoordsC, "C:%u", QPosC+1);

	uint PPL = QPosC + CL - QPosA;
	Ps(CoordStr, "%-12.12s   %-14.14s   %-8.8s  Palmprint(%u-%u, %u aa)",
	  CoordsA.c_str(), CoordsB.c_str(), CoordsC.c_str(),
	  QPosA+1, QPosC+CL, PPL);

	fprintf(f, "\n");
	fprintf(f, "%s\n", CoordStr.c_str());
	fprintf(f, "%s\n", QRow.c_str());
	fprintf(f, "%s\n", AnnotRow.c_str());
	fprintf(f, "%s\n", RRow.c_str());
	}

void TSHit::WriteTsv(FILE *f) const
	{
	if (f == 0)
		return;

	static bool HdrDone = false;
	if (!HdrDone)
		{
		fprintf(f, "Query");
		fprintf(f, "\tRef");
		fprintf(f, "\tScore");
		fprintf(f, "\tRMSDm");
		fprintf(f, "\tPosA");
		fprintf(f, "\tPosB");
		fprintf(f, "\tPosC");
		fprintf(f, "\tMotifA");
		fprintf(f, "\tMotifB");
		fprintf(f, "\tMotifC");
		fprintf(f, "\n");
		HdrDone = true;
		}

	asserta(m_Query != 0);
	asserta(m_Ref != 0);

	string QLabel;
	string RLabel;
	m_Query->GetLabel(QLabel);
	m_Ref->GetLabel(RLabel);

	double MotifRMSD = sqrt(m_MotifRMSD2);
	double ConfScore = GetConfidenceScore(MotifRMSD);

	string QASeq, QBSeq, QCSeq;
	m_Query->GetSubSeq(m_QPosA, AL, true, QASeq);
	m_Query->GetSubSeq(m_QPosB, BL, true, QBSeq);
	m_Query->GetSubSeq(m_QPosC, CL, true, QCSeq);

	fprintf(f, "%s", QLabel.c_str());
	fprintf(f, "\t%s", RLabel.c_str());
	fprintf(f, "\t%.3f", ConfScore);
	fprintf(f, "\t%.3g", MotifRMSD);
	fprintf(f, "\t%u", m_QPosA + 1);
	fprintf(f, "\t%u", m_QPosB + 1);
	fprintf(f, "\t%u", m_QPosC + 1);
	fprintf(f, "\t%s", QASeq.c_str());
	fprintf(f, "\t%s", QBSeq.c_str());
	fprintf(f, "\t%s", QCSeq.c_str());
	fprintf(f, "\n");
	}

void TSHit::WritePalmprintFasta(FILE *f) const
	{
	if (f == 0)
		return;
	asserta(m_Query != 0);
	const string &Q = m_Query->m_Seq;
	const uint QL = SIZE(Q);
	uint PPLo = m_QPosA;
	uint PPHi = m_QPosC + CL - 1;
	uint PPL = PPHi - PPLo + 1;
	asserta(PPHi > PPLo && PPHi < QL);
	string Label;
	m_Query->GetLabel(Label);
	SeqToFasta(f, Label.c_str(), Q.c_str() + PPLo, PPL);
	}

void TSHit::WritePalmprintPDB(const string &FileNamePrefix) const
	{
	if (FileNamePrefix == "")
		return;
	asserta(m_Query != 0);
	vector<double> t;
	vector<vector<double> > R;

	vector<vector<double> > MotifCoords(3);
	m_Query->GetPt(m_QPosA, MotifCoords[A]);
	m_Query->GetPt(m_QPosB, MotifCoords[B]);
	m_Query->GetPt(m_QPosC, MotifCoords[C]);

	GetTriForm(MotifCoords, t, R);

	PDBChain XChain;
	m_Query->CopyTriForm(t, R, XChain);

	string Label;
	m_Query->GetLabel(Label);
	string FileName = FileNamePrefix + Label;
	if (!EndsWith(FileName, ".pdb"))
		FileName += ".pdb";
	XChain.ToPDB(FileName);
	}

void TSHit::WriteSketch(FILE* f) const
	{
	if (f == 0)
		return;

	const string& QSS = m_Query->m_SS;
	const string& RSS = m_Ref->m_SS;
	//fprintf(f, "QSS %s\n", QSS.c_str());
	//fprintf(f, "RSS %s\n", RSS.c_str());

	uint QL = m_Query->GetSeqLength();
	uint RL = m_Ref->GetSeqLength();

	asserta(m_QPosA + AL < m_QPosB);
	asserta(m_QPosB + BL < m_QPosC);
	asserta(m_QPosC + CL <= QL);

	asserta(m_RPosA + AL < m_RPosB);
	asserta(m_RPosB + BL < m_RPosC);
	asserta(m_RPosC + CL <= RL);

	string QASk = QSS.substr(m_QPosA, AL);
	string RASk = RSS.substr(m_RPosA, AL);

	string QBSk = QSS.substr(m_QPosB, BL);
	string RBSk = RSS.substr(m_RPosB, BL);

	string QCSk = QSS.substr(m_QPosC, CL);
	string RCSk = RSS.substr(m_RPosC, CL);

	uint QV1L = m_QPosB - (m_QPosA + AL);
	uint RV1L = m_RPosB - (m_RPosA + AL);

	uint QV2L = m_QPosC - (m_QPosB + BL);
	uint RV2L = m_RPosC - (m_RPosB + BL);

	string QV1Sk;
	string RV1Sk;
	GetPalmSubSketch(QSS, m_QPosA + 1, QV1L, V_SKETCH_LENGTH, QV1Sk);
	GetPalmSubSketch(RSS, m_RPosA + 1, RV1L, V_SKETCH_LENGTH, RV1Sk);

	string QV2Sk;
	string RV2Sk;
	GetPalmSubSketch(QSS, m_QPosB + 1, QV2L, V_SKETCH_LENGTH, QV2Sk);
	GetPalmSubSketch(RSS, m_RPosB + 1, RV2L, V_SKETCH_LENGTH, RV2Sk);

	string QSketch = QASk + "  " + QV1Sk + "  " + QBSk + "  " + QV2Sk + "  " + QCSk;
	string RSketch = RASk + "  " + RV1Sk + "  " + RBSk + "  " + RV2Sk + "  " + RCSk;
	string AnnotSketch;
	MakeAnnotSketch(QSketch, RSketch, AnnotSketch);

	fprintf(f, "\n");
	fprintf(f,
"_____A______  _________V1_________  ______B_______  _________V2_________  ___C____\n");

	fprintf(f, "%s  %s\n", QSketch.c_str(), m_Query->m_ChainLabel.c_str());
	fprintf(f, "%s\n", AnnotSketch.c_str());
	fprintf(f, "%s  %s\n", RSketch.c_str(), m_Ref->m_ChainLabel.c_str());
	}

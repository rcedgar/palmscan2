#include "myutils.h"
#include "structprof.h"
#include "chainreader.h"
#include "outputfiles.h"
#include "shapesearcher.h"
#include "motifsettings.h"
#include "abcxyz.h"

#define TRACE	1

void GetTriForm(vector<vector<double> > &MotifCoords,
  vector<double> &t, vector<vector<double> > &R);

static void CavityPic(const PDBChain &Chain,
  uint PosA, uint PosB, uint PosC);

void cmd_cavitypic()
	{
	const string &InputFN = g_Arg1;

	Shapes S;
	S.InitFromCmdLine();

	ShapeSearcher SS;
	SS.Init(S);

	ChainReader CR;
	CR.Open(InputFN);

	PDBChain Chain;
	PDBChain XChain;
	bool Done = false;
	while (CR.GetNext(Chain))
		{
		SS.SetQuery(Chain);
		SS.SearchABC();
		if (SS.m_ABCScore < SS.m_MinABCScore)
			continue;

		uint PosA = UINT_MAX;
		uint PosB = UINT_MAX;
		uint PosC = UINT_MAX;
		SS.GetPosABC(PosA, PosB, PosC);
		string A, B, C;
		SS.GetA(A);
		SS.GetB(B);
		SS.GetC(C);
		Log(">%s\n", Chain.m_Label.c_str());
		Log("A=%s B=%s C=%s\n", A.c_str(), B.c_str(), C.c_str());

		Chain.SetMotifPosVec(PosA, PosB, PosC);
		Chain.GetTriFormChain_DGD(XChain);
		CavityPic(XChain, PosA, PosB, PosC);
		}
	}

static double LoX;
static double LoY;
static double LoZ;
static double HiX;
static double HiY;
static double HiZ;

static void WriteSvgHeader(double Width, double Height)
	{
	if (g_fsvg == 0)
		return;

	fprintf(g_fsvg, "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"");
	fprintf(g_fsvg, " width=\"%.6g\" height=\"%.6g\">\n", Width, Height);
	}

static void WriteSvgFooter()
	{
	if (g_fsvg == 0)
		return;
	fprintf(g_fsvg, "</svg>\n");
	}

static void UpdateBoundingCube(const PDBChain &Chain, uint Pos,
  double &LoX, double &HiX, double &LoY, double &HiY, double &LoZ, double &HiZ)
	{
	double X = Chain.GetX(Pos);
	double Y = Chain.GetY(Pos);
	double Z = Chain.GetZ(Pos);
	LoX = min(X, LoX);
	HiX = max(X, HiX);

	LoY = min(Y, LoY);
	HiY = max(Y, HiY);

	LoZ = min(Z, LoZ);
	HiZ = max(Z, HiZ);
	}

// Out(r'<circle cx="%.6g" cy="%.6g" r="%.6g" stroke-width="%.6g" stroke="%s" fill="%s" />'

static void DrawCircle(FILE *f, double x, double y, double r,
  double stroke_width, const string &stroke_color, const string &fill_color)
	{
	if (f == 0)
		return;
	fprintf(f, "<circle cx=\"%.6g\" cy=\"%.6g\" r=\"%.6g\"", x, y, r);
	fprintf(f, " stroke-width=\"%.6g\" stroke=\"%s\" fill=\"%s\" />\n",
	  stroke_width, stroke_color.c_str(), fill_color.c_str());
	}

// Out(r'<text x="%.6g" y="%.6g" font-family="%s" font-size="%.6gpx" fill="%s" text-anchor="%s">%s</text>' % (x, y, font, size_px, color, anchor, text))
static void DrawText(FILE *f, double x, double y, const string &Text)
	{
	if (f == 0)
		return;
	fprintf(f, "<text x=\"%.6g\" y=\"%.6g\"", x, y);
	fprintf(f, " fill=\"black\"");
	fprintf(f, " font-size=\"1px\"");
	fprintf(f, " text-anchor=\"center\" dominant-baseline=\"middle\"");
	fprintf(f, ">%s</text>\n", Text.c_str());
	}

static void DrawAtom(FILE *f, double x, double y, double z,
  const string &ElName, const string &AtomName, const string &Color)
	{
	if (ElName == "H")
		return;
	double cx = x - LoX;
	double cy = y - LoY;
	double r = 1.0/(1.0 + fabs(z));
	if (z > 0)
		DrawCircle(f, cx, cy, r, 0.1, "black", Color.c_str());
	else
		DrawCircle(f, cx, cy, r, 0.2, "orange", Color.c_str());
	//DrawText(f, cx, cy, AtomName);
	}

static void DrawAtoms(FILE *f, const PDBChain &Chain,
  uint Pos, const string &Color)
	{
	vector<double> Xs;
	vector<double> Ys;
	vector<double> Zs;
	vector<string> ElementNames;
	vector<string> AtomNames;
	vector<string> Lines;
	Chain.GetResidueAtomsInfo(Pos, Xs, Ys, Zs,
	  ElementNames, AtomNames, Lines);
	const uint N = SIZE(Xs);
	for (uint i = 0; i < N; ++i)
		{
		const string &ElName = ElementNames[i];
		const string &AtomName = AtomNames[i];
		double z = Zs[i];
		if (fabs(z) > 2)
			continue;
		double r = 1.0/(1.0 + fabs(z));
		Log("Pos %4u %3.3s %3.3s (%10.1f, %10.1f, %10.1f) %s r=%.3g\n",
		  Pos, ElName.c_str(), AtomName.c_str(), Xs[i], Ys[i], Zs[i], Color.c_str(), r);
		DrawAtom(f, Xs[i], Ys[i], z, ElName, AtomName, Color);
		}
	}

static void CavityPic(const PDBChain &Chain,
  uint PosA, uint PosB, uint PosC)
	{
// Bounding cube
	LoX = Chain.GetX(PosA);
	LoY = Chain.GetY(PosA);
	LoZ = Chain.GetZ(PosA);
	HiX = LoX;
	HiY = LoY;
	HiZ = LoZ;

	for (uint i = 1; i < g_LA; ++i)
		UpdateBoundingCube(Chain, PosA + i, LoX, HiX, LoY, HiY, LoZ, HiZ);

	for (uint i = 0; i < g_LB; ++i)
		UpdateBoundingCube(Chain, PosB + i, LoX, HiX, LoY, HiY, LoZ, HiZ);

	for (uint i = 0; i < g_LC; ++i)
		UpdateBoundingCube(Chain, PosC + i, LoX, HiX, LoY, HiY, LoZ, HiZ);

	//Log("Bounding cube X %.1f - %.1f, Y %.1f - %.1f, Z %.1f - %.1f\n",
	//  LoX, HiX, LoY, HiY, LoZ, HiZ);

	const double MARGIN = 2.0;
	double Width = HiX - LoX + 2*MARGIN;
	double Height = HiY - LoY + 2*MARGIN;
	WriteSvgHeader(Width, Height);
	//DrawAtoms(g_fsvg, Chain, PosA + g_OffAd, "blue");
	//DrawAtoms(g_fsvg, Chain, PosB + g_OffBg, "green");
	//DrawAtoms(g_fsvg, Chain, PosC + g_OffCd, "red");

	//DrawAtoms(g_fsvg, Chain, PosB + 5, "limegreen");
	//DrawAtoms(g_fsvg, Chain, PosB + 6, "limegreen");
	//DrawAtoms(g_fsvg, Chain, PosB + 7, "limegreen");

	//DrawAtoms(g_fsvg, Chain, PosA + 9, "cyan");
	//DrawAtoms(g_fsvg, Chain, PosA + 10, "cyan");
	//DrawAtoms(g_fsvg, Chain, PosA + 11, "cyan");

	for (uint i = 0; i < g_LA; ++i)
		if (i == g_OffAd)
			DrawAtoms(g_fsvg, Chain, PosA + i, "cyan");
		else
			DrawAtoms(g_fsvg, Chain, PosA + i, "blue");

	for (uint i = 0; i < g_LB; ++i)
		if (i == g_OffBg)
			DrawAtoms(g_fsvg, Chain, PosB + i, "limegreen");
		else
			DrawAtoms(g_fsvg, Chain, PosB + i, "green");

	for (uint i = 0; i < g_LC; ++i)
		if (i == g_OffCd)
			DrawAtoms(g_fsvg, Chain, PosC + i, "magenta");
		else
			DrawAtoms(g_fsvg, Chain, PosC + i, "red");

	WriteSvgFooter();
	}

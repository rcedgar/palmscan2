#include "myutils.h"
#include "structprof.h"
#include "chainreader.h"
#include "outputfiles.h"
#include "shapesearcher.h"
#include "motifsettings.h"
#include "abcxyz.h"

#define TRACE	0

void GetTriForm(
  vector<vector<double> > &MotifCoords,
  vector<double> &t,
  vector<vector<double> > &R)
	{
#if TRACE
	Log("\n");
	Log("_________________________________________________________\n");
	Log("GetTriForm()\n");
	LogMx("MotifCoords", MotifCoords);
#endif

	t.resize(3);
	R.resize(3);

	vector<vector<double> > Basis;
	vector<double> CentroidCoords;
	GetTriangleBasis(MotifCoords, CentroidCoords, Basis);
	GetBasisR(Basis, R);

	t[X] = -CentroidCoords[X];
	t[Y] = -CentroidCoords[Y];
	t[Z] = -CentroidCoords[Z];

#if DEBUG
	{
	vector<vector<double> > XMotifCoords;
	XFormMx(MotifCoords, t, R, XMotifCoords);

	vector<double> XCentroidCoords;
	GetTriangleCentroid(XMotifCoords, XCentroidCoords);
#if TRACE
	LogMx("XMotifCoords", XMotifCoords);
	LogVec("XCentroidCoords", XCentroidCoords);
#endif

	assert(feq(XCentroidCoords[0], 0));
	assert(feq(XCentroidCoords[1], 0));
	assert(feq(XCentroidCoords[2], 0));

	assert(feq(XMotifCoords[A][Z], 0));
	assert(feq(XMotifCoords[B][Z], 0));
	assert(feq(XMotifCoords[C][Z], 0));

	assert(feq(XMotifCoords[B][Y], 0));
	}
#endif
	}

static void TestTriForm(
  double Ax, double Ay, double Az,
  double Bx, double By, double Bz,
  double Cx, double Cy, double Cz)
	{
	vector<vector<double> > MotifCoords(3);

	MotifCoords[A].push_back(Ax);
	MotifCoords[A].push_back(Ay);
	MotifCoords[A].push_back(Az);

	MotifCoords[B].push_back(Bx);
	MotifCoords[B].push_back(By);
	MotifCoords[B].push_back(Bz);

	MotifCoords[C].push_back(Cx);
	MotifCoords[C].push_back(Cy);
	MotifCoords[C].push_back(Cz);

	vector<double> t;
	vector<vector<double> > R;
	GetTriForm(MotifCoords, t, R);
	}

static void GetRandomMotifCoords(vector<vector<double> > &MotifCoords)
	{
	Resize3x3(MotifCoords);
	for (uint i = 0; i < 3; ++i)
		{
		for (uint j = 0; j < 3; ++j)
			{
			uint r = randu32()%1000000;
			double d = r/1e6;
			double Coord = 10*(d - 0.5);
			MotifCoords[i][j] = Coord;
			}
		}
	}

void cmd_test_triform()
	{
	if (1)
		TestTriForm(
		  0, 0, 0,
		  1, 1, 0,
		  1, -1, 0);

	if (1)
		TestTriForm(
		  0, 0, 0,
		  0, 1, 1,
		  0, 1, -1);

	if (1)
		TestTriForm(
		  -13.574, 7.392, 0.666,
		  5.691, 1.001, -4.517,
		  -7.191, 3.889, 2.399);

	for (uint Iter = 0; Iter < 100; ++Iter)
		{
		vector<vector<double> > MotifCoords;
		GetRandomMotifCoords(MotifCoords);

		vector<double> t;
		vector<vector<double> > R;
		GetTriForm(MotifCoords, t, R);
		}
	}

void cmd_triform()
	{
	const string &InputFN = opt_triform;

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

		Chain.SetMotifPosVec(PosA, PosB, PosC);
		Chain.GetTriFormChain_DGD(XChain);
		XChain.ToPDB(opt_pdb);
		Done = true;
		break;
		}
	if (!Done)
		Warning("Motifs not found");
	}

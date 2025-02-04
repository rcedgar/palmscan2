#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"

void LogCoords(const char *Name, coords_t c);
coords_t normalize(const coords_t &a);
coords_t cross_product(const coords_t &a, const coords_t &b);
double get_angle(const coords_t &a, const coords_t &b);
coords_t subtract(const coords_t &a, const coords_t &b);
coords_t rotate_around_vector(const coords_t &v,
							  const coords_t &axis, double theta);

void cmd_overlap()
	{
	PDBChain A, B;
	A.FromCal(g_Arg1);
	B.FromCal(opt_input1);

	A.ToPDB("inputa.pdb");
	B.ToPDB("inputb.pdb");

	const uint LA = A.GetSeqLength();
	const uint LB = B.GetSeqLength();

	const string Coords = string(opt_motif_coords);
	vector<string> Fields;
	Split(Coords, Fields, ',');
	asserta(SIZE(Fields) == 4);

	const uint LoA = StrToUint(Fields[0]);
	const uint HiA = StrToUint(Fields[1]);
	const uint LoB = StrToUint(Fields[2]);
	const uint HiB = StrToUint(Fields[3]);
	asserta(LoA < HiA && HiA < LA);
	asserta(LoB < HiB && HiB < LB);

	double CentroidAX = (A.m_Xs[HiA] + A.m_Xs[LoA])/2;
	double CentroidAY = (A.m_Ys[HiA] + A.m_Ys[LoA])/2;
	double CentroidAZ = (A.m_Zs[HiA] + A.m_Zs[LoA])/2;

	double CentroidBX = (B.m_Xs[HiB] + B.m_Xs[LoB])/2;
	double CentroidBY = (B.m_Ys[HiB] + B.m_Ys[LoB])/2;
	double CentroidBZ = (B.m_Zs[HiB] + B.m_Zs[LoB])/2;

	double dX = CentroidAX - CentroidBX;
	double dY = CentroidAY - CentroidBY;
	double dZ = CentroidAZ - CentroidBZ;

	ProgressLog("dX=%.1f, dY=%.1f, dZ=%.1f\n",
				dX, dY, dZ);

	PDBChain NewB;
	NewB.m_Label = B.m_Label + ".ov";
	NewB.m_Seq = B.m_Seq;
	for (uint PosB = 0; PosB < LB; ++PosB)
		{
		double newx = B.m_Xs[PosB] + dX;
		double newy = B.m_Ys[PosB] + dY;
		double newz = B.m_Zs[PosB] + dZ;

		NewB.m_Xs.push_back(newx);
		NewB.m_Ys.push_back(newy);
		NewB.m_Zs.push_back(newz);
		}

	double NewCentroidBX = (NewB.m_Xs[HiB] + NewB.m_Xs[LoB])/2;
	double NewCentroidBY = (NewB.m_Ys[HiB] + NewB.m_Ys[LoB])/2;
	double NewCentroidBZ = (NewB.m_Zs[HiB] + NewB.m_Zs[LoB])/2;

	asserta(feq(NewCentroidBX, CentroidAX));
	asserta(feq(NewCentroidBY, CentroidAY));
	asserta(feq(NewCentroidBZ, CentroidAZ));

	NewB.ToPDB("newb1.pdb");

// Now A and B have a common centroid
	coords_t A1 = A.GetCoords(LoA);
	coords_t A2 = A.GetCoords(HiA);
	coords_t av = subtract(A2, A1);

	coords_t B1 = NewB.GetCoords(LoB);
	coords_t B2 = NewB.GetCoords(HiB);
	coords_t bv = subtract(B2, B1);

	coords_t axis = normalize(cross_product(av, bv));

	double theta_rad = get_angle(av, bv);
	double theta_deg = degrees(theta_rad);
	Log("theta = %.1f\n", theta_deg);

	coords_t bv_rot = rotate_around_vector(bv, axis, theta_rad);
	coords_t bv_rot_minus = rotate_around_vector(bv, axis, -theta_rad);
	LogCoords("av", av);
	LogCoords("bv", bv);
	LogCoords("bv_rot", bv_rot);
	LogCoords("bv_rot_minus", bv_rot_minus);

	PDBChain NewB2;
	NewB2.m_Label = B.m_Label + ".ov";
	NewB2.m_Seq = B.m_Seq;
	for (uint PosB = 0; PosB < LB; ++PosB)
		{
		coords_t pt = NewB.GetCoords(PosB);
		coords_t pt_rot = rotate_around_vector(pt, axis, theta_rad);

		NewB2.m_Xs.push_back(pt.x);
		NewB2.m_Ys.push_back(pt.y);
		NewB2.m_Zs.push_back(pt.z);
		}
	NewB2.ToPDB("newb2.pdb");
	}

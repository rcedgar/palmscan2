#ifndef MY_VERSION
#define MY_VERSION	"2.0"
#endif

#define C(x)	STR_OPT(x)
#include "cmds.h"

#define F(x)	STR_OPT(x)
#include "ofiles.h"

STR_OPT(search_aa_top)
STR_OPT(search_nt_top)
STR_OPT(search_nt_all)
STR_OPT(search_nt_topx)
STR_OPT(search_pair_nt)
STR_OPT(search_pair_g)
STR_OPT(search_pair_aa)
STR_OPT(otumaps)
STR_OPT(alnqc)
STR_OPT(msaqc)
STR_OPT(msaqc3)
STR_OPT(logpssms)
STR_OPT(agreefastaout)
STR_OPT(qhitfastaout)
STR_OPT(rex)
STR_OPT(align_msas)
STR_OPT(test_pssms)
STR_OPT(ppp)
STR_OPT(tree)
STR_OPT(splitb)
STR_OPT(shapes)
STR_OPT(loaddir)

STR_OPT(train_cal)
STR_OPT(test_cal)

STR_OPT(ref)
STR_OPT(label)
STR_OPT(labels)
STR_OPT(fastaout)
STR_OPT(alnout)
STR_OPT(triout)
STR_OPT(calout)
STR_OPT(caloutx)
STR_OPT(pdboutx)
STR_OPT(chain)
STR_OPT(db)
STR_OPT(freqsout)
STR_OPT(relabel)
STR_OPT(fevout)
STR_OPT(bedout)
STR_OPT(model)
STR_OPT(ppout)
STR_OPT(ppout_nt)
STR_OPT(motifs_fastaout)
STR_OPT(motifs_fastaout2)
STR_OPT(motifs_fastaout3)
STR_OPT(nohit_fastaout)
STR_OPT(pssm_alnout)
STR_OPT(msa2)
STR_OPT(pssmdir)
STR_OPT(pdbout)
STR_OPT(xxx)
STR_OPT(centroids)
STR_OPT(exclude)
STR_OPT(include)
STR_OPT(seg_fasta_prefix)
STR_OPT(motif_coords)
STR_OPT(shapename)

STR_OPT(log)
STR_OPT(input)
STR_OPT(input1)
STR_OPT(output)
STR_OPT(output2)
STR_OPT(output3)
STR_OPT(psm)
STR_OPT(psm1)
STR_OPT(psm2)
STR_OPT(spec)
STR_OPT(logodds)

STR_OPT(fw_name)
STR_OPT(psm_fw)

STR_OPT(gapopen)
STR_OPT(gapext)
STR_OPT(frame)
STR_OPT(test)
STR_OPT(sep)
STR_OPT(fieldnrs)
STR_OPT(png)
STR_OPT(subchain)
STR_OPT(searchmfs)
STR_OPT(requiremfs)

UNS_OPT(threads,			8,			0,			UINT_MAX)
UNS_OPT(band,				16,			0,			UINT_MAX)
UNS_OPT(hspw,				0,			1,			UINT_MAX)
UNS_OPT(minhsp,				32,			1,			UINT_MAX)
UNS_OPT(iddrop,				8,			1,			UINT_MAX)
UNS_OPT(cluster_maxdiffs,	1,			0,			UINT_MAX)
UNS_OPT(tsv_topn,			32,			1,			UINT_MAX)
UNS_OPT(report_topn,		32,			1,			UINT_MAX)
UNS_OPT(maxx,				10,			1,			UINT_MAX)
UNS_OPT(refn,				10,			1,			UINT_MAX)
UNS_OPT(minflanklen,		10,			1,			UINT_MAX)
UNS_OPT(scoretype,			10,			1,			UINT_MAX)
UNS_OPT(minlen,				10,			1,			UINT_MAX)
UNS_OPT(lo,					10,			1,			UINT_MAX)
UNS_OPT(hi,					10,			1,			UINT_MAX)
UNS_OPT(pos,				10,			1,			UINT_MAX)

UNS_OPT(secs,				60,			1,			UINT_MAX)
UNS_OPT(flanks,				60,			1,			UINT_MAX)

UNS_OPT(maxseqlength,		500000000,		1,			UINT_MAX)
UNS_OPT(sfasta_buff_bytes,	512*1024*1024,1024,		UINT_MAX)
UNS_OPT(randseed,			0,			0,			UINT_MAX)
UNS_OPT(usort_w,			6,			1,			UINT_MAX)
UNS_OPT(mingap,				0,			1,			UINT_MAX)
UNS_OPT(maxgap,				999,		1,			UINT_MAX)
UNS_OPT(hiw,				45,			1,			UINT_MAX)
UNS_OPT(low,				16,			1,			UINT_MAX)
UNS_OPT(sample_size,		16,			1,			UINT_MAX)
UNS_OPT(topn,				16,			1,			UINT_MAX)

FLT_OPT(minscore,			10.0,		-9e9,		+9e9)
FLT_OPT(minppscore,			10.0,		-9e9,		+9e9)
FLT_OPT(minpalmscore,		10.0,		-9e9,		+9e9)
//FLT_OPT(maxlefppm,			10.0,		-9e9,		+9e9)
FLT_OPT(maxscore,			10.0,		-9e9,		+9e9)
FLT_OPT(motif_cluster_minscore,			0.0,		-9e9,		+9e9)
FLT_OPT(minscore1,			0.0,		-9e9,		+9e9)
FLT_OPT(minscore2,			0.0,		-9e9,		+9e9)
FLT_OPT(minscore_pair,		0.0,		-9e9,		+9e9)
FLT_OPT(stop_score,			-10.0,		-9e9,		+9e9)
FLT_OPT(minscore_fw,		0.0,		-9e9,		+9e9)
FLT_OPT(fastq_maxee,		1.0,		0.0,		9e9)
FLT_OPT(threshold,			1.0,		0.0,		9e9)
FLT_OPT(minselfscorenonpp,	1.0,		0.0,		9e9)
FLT_OPT(minselfscorepp,		1.0,		0.0,		9e9)

FLT_OPT(match,				1.0,		0.0,		DBL_MAX)
FLT_OPT(mismatch,			-2.0,		0.0,		DBL_MAX)
FLT_OPT(xdrop_u,			16.0,		0.0,		DBL_MAX)
FLT_OPT(xdrop_g,			32.0,		0.0,		DBL_MAX)
FLT_OPT(xdrop_nw,			16.0,		0.0,		DBL_MAX)
FLT_OPT(lopen,				10.0,		0.0,		DBL_MAX)
FLT_OPT(lext,				1.0,		0.0,		DBL_MAX)
FLT_OPT(mincscore,			0.0,		0.0,		DBL_MAX)
FLT_OPT(rmsd,				0.0,		0.0,		DBL_MAX)
FLT_OPT(train_pct,			0.0,		0.0,		DBL_MAX)
FLT_OPT(min_palm_score,			0.0,		0.0,		DBL_MAX)
FLT_OPT(radius,			0.0,		0.0,		DBL_MAX)

FLAG_OPT(trunclabels)
FLAG_OPT(notrunclabels)
FLAG_OPT(compilerinfo)
FLAG_OPT(quiet)
FLAG_OPT(logmemgrows)
FLAG_OPT(fulldp)
FLAG_OPT(logquery)
FLAG_OPT(verbose)
FLAG_OPT(gapped)
FLAG_OPT(coords)
FLAG_OPT(permuted)
FLAG_OPT(notpermuted)
FLAG_OPT(trace)
FLAG_OPT(self)
FLAG_OPT(leave_one_out)
FLAG_OPT(label_substr_match)
FLAG_OPT(pdb_ent)
FLAG_OPT(dgd)
FLAG_OPT(top_hit_only)
FLAG_OPT(refs)
FLAG_OPT(first_only)
FLAG_OPT(norows)
FLAG_OPT(calctps)
FLAG_OPT(compact)
FLAG_OPT(dotproduct)
FLAG_OPT(subchainc)
FLAG_OPT(delete_hetatm)
FLAG_OPT(misses)
FLAG_OPT(pml_savepse)
FLAG_OPT(traceabc)
FLAG_OPT(calibrate)
//FLAG_OPT(sensitive)
//FLAG_OPT(lowerrors)
FLAG_OPT(pmlrevmotifs)
FLAG_OPT(noappendchain)

#undef FLAG_OPT
#undef UNS_OPT
#undef FLT_OPT
#undef STR_OPT


const char *help_txt[] =
	{
	"CHAINSPEC",
	"    Command-line argument which specifies one or more structure files.",
	"    If CHAINSPEC has .cal or .ppc extension, it must be a CAL file.",
	"    If CHAINSPEC has .files extension, it must be a text file with",
	"      one filename or pathname per line, these must be PDB files.",
	"    Otherwise, it must be a filename or pathname for a file",
	"      in PDB format (typically .pdb or .ent extension).",
	"",
	"CAL",
	"    C-alpha format. Compact file format for one or more chains.",
	"    Each chain starts with FASTA-like label.",
	"    Following lines have four tab-separated fields: one-letter amino",
	"      acid code, then X, Y, Z coordinates for its C-alpha atom.",
	"",
	"PPC",
	"    Palmprint coordinates. Special case of CAL where",
	"      (1) chain is trimmed to a palmprint,",
	"      (2) ABC motif coordinates are specified in the label,",
	"      (3) structure is rotated to motif-anchored coordinates.",
	"    All PPCs are pre-aligned to each other, enabling very fast",
	"      pair-wise comparisons.",
	"",
	"MODELFILE",
	"    File in .ppm format containing PSSMs for one or more groups.",
	"",
	"Convert structure files to CAL format",
	"    palmscan2 -pdb2cal CHAINSPEC",
	"        -cal chains.cal",
	"",
	"Search FASTA for palmprints using PSSMs",
	"    palmscan2 -search_pssms seqs.fasta ",
	"        -model MODELFILE (required)",
	"        -report_pssms report.txt",
	"        -tsv hits.txt",
	"        -threads N (default min(number of cores, 10))",
	"",
	"Search CAL for palmprints using PSSMs",
	"    palmscan2 -search3d_pssms chains.cal",
	"        -model MODELFILE (required)",
	"        -ppc hits.ppc",
	"        -report_pssms report.txt",
	"        -tsv hits.txt",
	"        -threads N (default min(number of cores, 10))",
	"",
	"Search CAL for palmprints using PPCs",
	"    palmscan2 -search3d_cal chains.cal",
	"        -ref db.ppc (required)",
	"        -report_3d report.txt",
	"        -tsv hits.txt",
	"        -ppfa palmprints.fasta",
	"        -ppc palmprints.ppc",
	"        -threads N (default min(number of cores, 10))",
	"",
	"Align all vs. all PPCs",
	"    palmscan2 -search3d_ppc query.ppc",
	"        -ref db.ppc (required)",
	"        -report_3d report.txt",
	"        -tsv hits.txt (query, target, RMSD)",
	"        -threads N (default min(number of cores, 10))",
	"",
	"Cluster PPCs at given threshold (RMSD in Angstroms)",
	"    palmscan2 -cluster_ppc input.ppc",
	"	    -rmsd 2.0",
	"        -tsv input_rmsd2.tsv",
	"        -ppc centroids.rmsd2.ppc",
	"",
	"Classify PPCs as RdRp / not RdRp. Labels must start",
	"  with the five characters 'rdrp.' if they are RdRp.",
	"    palmscan2 -classify query.ppc",
	"        -ref db.ppc (required)",
	"        -tsv hits.txt",
	"        -threads N (default min(number of cores, 10))",
	"",
	"Remove chains with identical sequences",
	"    palmscan2 -remove_dupes CHAINSPEC",
	"        -cal chains.nodupes.cal",
	"",
	"Convert CAL, PPC or PDB to FASTA",
	"    palmscan2 -cal2fa CHAINSPEC",
	"        -fasta chains.fasta",
	"",
	"Extract chains with given labels",
	"    palmscan -getchains CHAINSPEC",
	"        -labels labels.txt (required, one label per line)",
	"        -cal subset.cal",
	};
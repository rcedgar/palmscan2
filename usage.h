
const char *usage_txt[] =
	{
	"3D search:",
	"    palmscan2 -search3d INPUT -tsv hits.tsv",
	"",
	"Convert PDB to CAL format:",
	"    palmscan2 -pdb2cal INPUT -cal chains.cal",
	"",
	"INPUT:",
	"    (1) PDB format file with one or more chains",
	"    (2) CAL format file with one or more chains",
	"    (3) Text file with one PDB or CAL pathname per line",
	"",
	"Options:",
	"    -tsv hits.tsv   Tab-separated file with hits",
	"    -misses         Include non-matches in tsv (default hits only)",
	"    -threads T      Run T threads (default nr CPU cores)",
	"    -pml hits.pml   Output pymol script to annotate chain(s)",
	"    -pml_save       Include save commands in pml script to generate",
	"                      one pse for each input chain",
	"    -loaddir DIR    Directory where pdb files are stored so that",
	"                      load commands in pml scripts work as-is",
	"    -searchmfs      Motifs to search for (default *)",
	"    -requiremfs     Required motifs (default ABC)",
	"    -minppscore     Minimum palmprint score 0..1 (default 0.5)",
	"    -minpalmscore   Minimum palm domain score 0..1 (default 0.69)",
	"    -minselfscorepp     Minimum ABC motif self-score (default 0.5)",
	"    -minselfscorenonpp  Minimum non-ABC self-score (default 0.5)",
	"",
	"Motifs",
	"    F1, F2, H, J, A, B, C and D are modeled.",
	"    The -searchmfs and -requiremfs options specify which",
	"    motifs should be searched and which should be required",
	"    to report a hit, respectively. ABC must be included.",
	"    Star (*) is wildcard, e.g. -searchmfs * means search",
	"    for all motifs.",
	"",
	"CAL format:",
	"    Sequence and C-alpha coordinates only.",
	"    More compact than PDB for faster search.",
	"    Use pdb2cal command to convert PDB to CAL.",
	};

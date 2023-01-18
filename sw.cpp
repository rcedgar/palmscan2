#include "myutils.h"

#if 0
// Textbook Smith-Waterman affine gap implementation.

#define	TRACE	0

static const char *LocalScoreToStr(double s)
	{
	static char str[16];
	if (MINUS_INFINITY == s)
		return "     *";
	sprintf(str, "%6.2f", s);
	return str;
	}

static void ListDP(const double *DPM_, const ProfPos *PA, const ProfPos *PB,
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	Log("        ");
	for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		char c = ' ';
		if (uPrefixLengthB > 0)
			c = ConsensusChar(PB[uPrefixLengthB - 1]);
		Log(" %4u:%c", uPrefixLengthB, c);
		}
	Log("\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			Log(" %s", LocalScoreToStr(DPM(uPrefixLengthA, uPrefixLengthB)));
		Log("\n");
		}
	}

double SW(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	assert(uLengthB > 0 && uLengthA > 0);

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;

// Allocate DP matrices
	const size_t LM = uPrefixCountA*uPrefixCountB;
	double *DPM_ = new double[LM];
	double *DPD_ = new double[LM];
	double *DPI_ = new double[LM];

	DPM(0, 0) = 0;
	DPD(0, 0) = MINUS_INFINITY;
	DPI(0, 0) = MINUS_INFINITY;

	DPM(1, 0) = MINUS_INFINITY;
	DPD(1, 0) = MINUS_INFINITY;
	DPI(1, 0) = MINUS_INFINITY;

	DPM(0, 1) = MINUS_INFINITY;
	DPD(0, 1) = MINUS_INFINITY;
	DPI(0, 1) = MINUS_INFINITY;

// Empty prefix of B is special case
	for (unsigned uPrefixLengthA = 2; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(uPrefixLengthA, 0) = MINUS_INFINITY;

	// D=LetterA+GapB, never optimal in local alignment with gap penalties
		DPD(uPrefixLengthA, 0) = MINUS_INFINITY;

	// I=GapA+LetterB, impossible with empty prefix
		DPI(uPrefixLengthA, 0) = MINUS_INFINITY;
		}

// Empty prefix of A is special case
	for (unsigned uPrefixLengthB = 2; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(0, uPrefixLengthB) = MINUS_INFINITY;

	// D=LetterA+GapB, impossible with empty prefix
		DPD(0, uPrefixLengthB) = MINUS_INFINITY;

	// I=GapA+LetterB, never optimal in local alignment with gap penalties
		DPI(0, uPrefixLengthB) = MINUS_INFINITY;
		}

	double scoreMax = MINUS_INFINITY;
	unsigned uPrefixLengthAMax = uInsane;
	unsigned uPrefixLengthBMax = uInsane;

// ============
// Main DP loop
// ============
	double scoreGapCloseB = MINUS_INFINITY;
	for (unsigned uPrefixLengthB = 1; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		const ProfPos &PPB = PB[uPrefixLengthB - 1];

		double scoreGapCloseA = MINUS_INFINITY;
		for (unsigned uPrefixLengthA = 1; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
			{
			const ProfPos &PPA = PA[uPrefixLengthA - 1];

			{
		// Match M=LetterA+LetterB
			double scoreLL = ScoreProfPos2(PPA, PPB);

			double scoreMM = DPM(uPrefixLengthA-1, uPrefixLengthB-1);
			double scoreDM = DPD(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapCloseA;
			double scoreIM = DPI(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapCloseB;

			double scoreBest;
			if (scoreMM >= scoreDM && scoreMM >= scoreIM)
				scoreBest = scoreMM;
			else if (scoreDM >= scoreMM && scoreDM >= scoreIM)
				scoreBest = scoreDM;
			else 
				{
				assert(scoreIM >= scoreMM && scoreIM >= scoreDM);
				scoreBest = scoreIM;
				}
			if (scoreBest < 0)
				scoreBest = 0;
			scoreBest += scoreLL;
			if (scoreBest > scoreMax)
				{
				scoreMax = scoreBest;
				uPrefixLengthAMax = uPrefixLengthA;
				uPrefixLengthBMax = uPrefixLengthB;
				}
			DPM(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			}

			{
		// Delete D=LetterA+GapB
			double scoreMD = DPM(uPrefixLengthA-1, uPrefixLengthB) +
			  PA[uPrefixLengthA-1].m_scoreGapOpen;
			double scoreDD = DPD(uPrefixLengthA-1, uPrefixLengthB);

			double scoreBest;
			if (scoreMD >= scoreDD)
				scoreBest = scoreMD;
			else
				{
				assert(scoreDD >= scoreMD);
				scoreBest = scoreDD;
				}
			DPD(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			}

		// Insert I=GapA+LetterB
			{
			double scoreMI = DPM(uPrefixLengthA, uPrefixLengthB-1) +
			  PB[uPrefixLengthB - 1].m_scoreGapOpen;
			double scoreII = DPI(uPrefixLengthA, uPrefixLengthB-1);

			double scoreBest;
			if (scoreMI >= scoreII)
				scoreBest = scoreMI;
			else 
				{
				assert(scoreII > scoreMI);
				scoreBest = scoreII;
				}
			DPI(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			}

			scoreGapCloseA = PPA.m_scoreGapClose;
			}
		scoreGapCloseB = PPB.m_scoreGapClose;
		}

#if TRACE
	Log("DPM:\n");
	ListDP(DPM_, PA, PB, uPrefixLengthA, uPrefixLengthB);
	Log("DPD:\n");
	ListDP(DPD_, PA, PB, uPrefixLengthA, uPrefixLengthB);
	Log("DPI:\n");
	ListDP(DPI_, PA, PB, uPrefixLengthA, uPrefixLengthB);
#endif

	assert(scoreMax == DPM(uPrefixLengthAMax, uPrefixLengthBMax));
	TraceBackSW(PA, uLengthA, PB, uLengthB, DPM_, DPD_, DPI_, 
	  uPrefixLengthAMax, uPrefixLengthBMax, Path);

#if	TRACE
	double scorePath = FastScorePath2(PA, uLengthA, PB, uLengthB, Path);
	Path.LogMe();
	Log("Score = %s Path = %s\n", LocalScoreToStr(scoreMax), LocalScoreToStr(scorePath));
#endif

	delete[] DPM_;
	delete[] DPD_;
	delete[] DPI_;

	return scoreMax;
	}
#endif // 0

#ifndef quarts_h
#define quarts_h

struct Quarts
	{
	unsigned Min;
	unsigned LoQ;
	unsigned Med;
	unsigned HiQ;
	unsigned Max;
	unsigned Total;
	double Avg;

	Quarts()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		}

	void Clear()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		}
	};

struct QuartsFloat
	{
	float Min;
	float LoQ;
	float Med;
	float HiQ;
	float Max;
	float Total;
	float Avg;
	float StdDev;

	QuartsFloat()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		StdDev = 0;
		}

	void Clear()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		StdDev = 0;
		}
	};

void GetQuarts(const vector<unsigned> &v, Quarts &Q);
void GetQuartsFloat(const vector<float> &v, QuartsFloat &Q);

#endif // quarts_h

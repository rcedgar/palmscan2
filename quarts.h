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

struct QuartsDouble
	{
	double Min;
	double LoQ;
	double Med;
	double HiQ;
	double Max;
	double Total;
	double Avg;
	double StdDev;

	QuartsDouble()
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
void GetQuartsDouble(const vector<double> &v, QuartsDouble &Q);

#endif // quarts_h

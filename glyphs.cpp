#include "myutils.h"
#include "alpha.h"

static vector<string> g_Paths;

static bool Init()
	{
#define X(c, path)	g_Paths.push_back(path);
#include "glyphs.h"
#undef X
	return true;
	}
static bool InitDone = Init();

const string &GetGlyphPath(uint Letter)
	{
	asserta(Letter < SIZE(g_Paths));
	return g_Paths[Letter];
	}

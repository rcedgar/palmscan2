#include "myutils.h"
#include "pdb.h"
#include "rdrpmodel.h"
#include "rdrpsearcher.h"

void cmd_scanpdb()
	{
	const string &QueryFN = opt_scanpdb;

	PDB Q;
	Q.FromFile(QueryFN);
	Q.LogMe();

	const string &QSeq = Q.m_Seq;
	const string &QLabel = Q.m_Label;

	asserta(optset_model);
	const string &ModelFileName = opt_model;
	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

	RdRpSearcher RS;
	RS.Init(Model);

	RdRpSearcher::InitOutput();
	RS.Search(QLabel, QSeq);
	RS.WriteOutput();
	RdRpSearcher::CloseOutput();
	}

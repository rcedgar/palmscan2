BINDIR := ../bin
OBJDIR := o
BINPATH := $(BINDIR)/palmscan2

CPPFLAGS := $(CPPFLAGS) -DNDEBUG -pthread

CXX = ccache g++
CXXFLAGS := $(CXXFLAGS) -O3 -fopenmp -ffast-math -msse -mfpmath=sse

UNAME_S := $(shell uname -s)
LDFLAGS := $(LDFLAGS) -O3 -fopenmp -pthread -lpthread
ifeq ($(UNAME_S),Linux)
    LDFLAGS += -static
endif

HDRS = \
  abcxyz.h \
  allocid.h \
  allocids.h \
  allocs.h \
  alnheuristics.h \
  alnparams.h \
  alpha.h \
  calreader.h \
  cmds.h \
  conncomp.h \
  diagbox.h \
  fastaseqsource.h \
  fastq.h \
  fastqrec.h \
  fileseqsource.h \
  gobuff.h \
  heuristics.h \
  hsp.h \
  linereader.h \
  lockobj.h \
  lockobjs.h \
  msaqc.h \
  msaqc2.h \
  mx.h \
  myopts.h \
  myutils.h \
  obj.h \
  objmgr.h \
  objtype.h \
  objtypes.h \
  ofiles.h \
  omplock.h \
  outputfiles.h \
  palmhit.h \
  pathinfo.h \
  pdbchain.h \
  ppcaligner.h \
  ppchit.h \
  pssm.h \
  pssms.h \
  pssmsearch.h \
  quarts.h \
  randomseqsource.h \
  rdrpmodel.h \
  rdrpsearcher.h \
  rphit.h \
  searchparams.h \
  seqdb.h \
  seqinfo.h \
  seqsource.h \
  sfasta.h \
  sort.h \
  timers.h \
  timing.h \
  tracebit.h \
  trisearcher.h \
  tshit.h \
  tshitmgr.h \
  usage.h \
  usearch.h \
  viterbi.h \
  xdpmem.h \
  xlat.h \

OBJS = \
  $(OBJDIR)/abcxyz.o \
  $(OBJDIR)/alignabc.o \
  $(OBJDIR)/alignpalm.o \
  $(OBJDIR)/align_msas.o \
  $(OBJDIR)/alloc.o \
  $(OBJDIR)/alpha.o \
  $(OBJDIR)/alpha2.o \
  $(OBJDIR)/blosum.o \
  $(OBJDIR)/blosum62.o \
  $(OBJDIR)/buildpsm.o \
  $(OBJDIR)/build_rdrp_model.o \
  $(OBJDIR)/cal2fa.o \
  $(OBJDIR)/calreader.o \
  $(OBJDIR)/classify.o \
  $(OBJDIR)/cluster_ppc.o \
  $(OBJDIR)/getchains.o \
  $(OBJDIR)/readppc.o \
  $(OBJDIR)/remove_dupes.o \
  $(OBJDIR)/cluster_cl.o \
  $(OBJDIR)/diagbox.o \
  $(OBJDIR)/fastaseqsource.o \
  $(OBJDIR)/fastqrec.o \
  $(OBJDIR)/fileseqsource.o \
  $(OBJDIR)/getseg.o \
  $(OBJDIR)/getss.o \
  $(OBJDIR)/help.o \
  $(OBJDIR)/invertmx.o \
  $(OBJDIR)/linereader.o \
  $(OBJDIR)/lockobj.o \
  $(OBJDIR)/logaln.o \
  $(OBJDIR)/model_strings.o \
  $(OBJDIR)/model_strings_rdrp.o \
  $(OBJDIR)/mx.o \
  $(OBJDIR)/myutils.o \
  $(OBJDIR)/ntmx.o \
  $(OBJDIR)/objmgr.o \
  $(OBJDIR)/one2three.o \
  $(OBJDIR)/output.o \
  $(OBJDIR)/outputfiles.o \
  $(OBJDIR)/palmscan2_main.o \
  $(OBJDIR)/getfilenames.o \
  $(OBJDIR)/palmsketch.o \
  $(OBJDIR)/pdb2cal.o \
  $(OBJDIR)/pdbchaincal.o \
  $(OBJDIR)/ppcaligner.o \
  $(OBJDIR)/readchains.o \
  $(OBJDIR)/search3d.o \
  $(OBJDIR)/pdbchain.o \
  $(OBJDIR)/psms.o \
  $(OBJDIR)/pssm.o \
  $(OBJDIR)/quarts.o \
  $(OBJDIR)/rdrpsearcher.o \
  $(OBJDIR)/search3d_cal.o \
  $(OBJDIR)/searchaatop.o \
  $(OBJDIR)/search_pssms.o \
  $(OBJDIR)/searchgroup.o \
  $(OBJDIR)/search3d_pssms.o \
  $(OBJDIR)/searchparams.o \
  $(OBJDIR)/search3d_ppc.o \
  $(OBJDIR)/searchtriangle.o \
  $(OBJDIR)/split_train_test.o \
  $(OBJDIR)/test_pssms.o \
  $(OBJDIR)/seqdb.o \
  $(OBJDIR)/seqinfo.o \
  $(OBJDIR)/seqsource.o \
  $(OBJDIR)/sfasta.o \
  $(OBJDIR)/sort.o \
  $(OBJDIR)/rdrpmodel.o \
  $(OBJDIR)/test.o \
  $(OBJDIR)/timing.o \
  $(OBJDIR)/tmscore.o \
  $(OBJDIR)/tracebackbitmem.o \
  $(OBJDIR)/triangle.o \
  $(OBJDIR)/triangles.o \
  $(OBJDIR)/triform.o \
  $(OBJDIR)/trisearcher.o \
  $(OBJDIR)/tshit.o \
  $(OBJDIR)/tshitmgr.o \
  $(OBJDIR)/tsoutput.o \
  $(OBJDIR)/validate.o \
  $(OBJDIR)/viterbifastmem.o \
  $(OBJDIR)/xbasis.o \

.PHONY: clean

$(BINPATH) : $(BINDIR)/ $(OBJDIR)/ $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(BINPATH)
	strip -d $(BINPATH)

$(OBJDIR)/ :
	mkdir -p $(OBJDIR)/

$(BINDIR)/ :
	mkdir -p $(BINDIR)/

$(OBJDIR)/%.o : %.cpp $(HDRS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR)/ $(BINPATH)
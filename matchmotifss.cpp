#include "myutils.h"
#include "denovosearcher.h"

/***
                                   1234567890123
                             012 3 4567890123456
            D                sss s s~~~hh~~~~hhh 
-------------------------------------------------------
                             012 3 4567890123456
  =A=  CVAT D VSDHDTFWPGWLR  sss s s~~~hh~~~~hhh 
  =A=  GFSY D TRCFDSTVTESDI  sss s ~~~~hh~~~thhh 
  =A=  GFSY D TRCFDSTVTESDI  sss s s~~~hh~~~thhh 
  =A=  AVSF D TKNWDTQVTSRDL  sss s s~~~hh~~~thhh 
  =A=  LFAF D YTGYDASLSPAWF  sss s s~~~hh~~~thhh 
  =A=  LATV D LSAASDSISLALC  sss s s~~~hh~~~thhh 
  =A=  HYDA D YSRWDSTQQRAVL  sss s s~~~hh~~~thhh 
  =A=  GIKV D IRDAYGNVKIPVL  sss s s~~~hh~~~~hhh 
  =A=  LATV D LSAASDSISLALC  sss s ~~~~hh~~~thhh 
  =A=  AFFG D VSSYDHSFSEEKI  sss s s~~~hh~~~thhh 
  =A=  AFFG D VSSYDHSFSEEKI  sss s s~~~hh~~~thhh 
***/
bool DeNovoSearcher::MatchAd_SS(uint Pos) const
	{
	asserta(Pos >= 3);
#define x(i, c) if (strchr(c, m_SS[Pos-3+i]) == 0) return false;
	x(0, "s")
	x(1, "s")
	x(2, "s")
	x(3, "s")
	x(4, "s~")
	x(5, "s~")
	x(6, "~")
	x(7, "~")
	x(8, "h")
	x(9, "h")
	x(10, "~")
	x(11, "~")
	x(12, "~th")
	x(13, "~th")
	x(14, "~th")
	x(15, "h")
	x(16, "h")
#undef x
	return true;
	}

/***
                          12345   12345678901
=B=  VGLSS G QGATDLMGTLL  ~~~~t ~ ~t~~hhhhhhh
=B=  RCRAS G VLTTSCGNTLT  ~t~~~ ~ ~thhhhhhhhh
=B=  RCRAS G VLTTSCGNTLT  ~~~~~ ~ ~t~hhhhhhhh
=B=  GQRGS G QPDTSAGNSML  ~~s~t ~ ~t~~hhhhhhh
=B=  GGMPS G CSGTSIFNSMI  ~~~s~ ~ ~t~~hhhhhhh
=B=  KISSM G NGYTFELESLI  t~~~t ~ ~tt~hhhhhhh
=B=  EGLPS G VPCTSQWNSIA  ~~~~~ ~ ~t~~hhhhhhh
=B=  HGLLQ G DPLSGCLCELY  ~~~st ~ ~thhhhhhhhh
=B=  KISSM G NGYTFELESLI  t~~~t ~ ~t~~hhhhhhh
***/
bool DeNovoSearcher::MatchBg_SS(uint Pos) const
	{
// Positions +5 to +11 must be h
	for (uint i = 5; i <= 11; ++i)
		if (m_SS[Pos+i] != 'h')
			return false;

// Must be a turn in -2 ... +2
	bool HasTurn = false;
	for (int i = -2; i <= 2; ++i)
		if (m_SS[int(Pos)+i] == 't')
			{
			HasTurn = true;
			break;
			}
	if (!HasTurn)
		return false;

// Must be four ~s in range -5 .. 0
	uint n = 0;
	for (int i = -5; i <= 0; ++i)
		if (m_SS[int(Pos)+i] == '~')
			++n;
	if (n < 4)
		return false;
	return true;
	}

/***
                               123456   123456 
=C=  IRQISKS D DAMLGWTKGRALV  ~ssss~t t ~sssss~~~thhh  >1hhs_A-cyst_A
=C=  CTMLVNG D DLVVICESAGTQE  sssss~t t ~sssssss~~thh  >1nb7_A-hepc_A
=C=  ARIHVCG D DGFLITEKGLGLK  ~ssss~h t ~ssss~~thhhhh  >2cjq_A-pest_A
=C=  LKMIAYG D DVIASYPHEVDAS  ~ssss~h t ~ssss~tt~s~~h  >2ijd_2-poli_1
=C=  SEVTVYG D DIILPSCAVPALR  h~~ss~h t ~ss~~th~hhhhh  >3avt_A-qbta_A
=C=  SLFSFYG D DEIVSTDIKLDPE  ~ssss~h t ~ssss~~~ss~~h  >3bso_A-noro_A
=C=  AFIHRTV D DYFFCSPHPHKVY  ~ssss~h t ~ssss~~~~hhhh  >3du5_A-tert_A
***/
bool DeNovoSearcher::MatchCd_SS(uint Pos) const
	{
// Must be turn at -1 0 or +1
	if (m_SS[Pos-1] != 't' && m_SS[Pos] != 't' && m_SS[Pos+1] != 't')
		return false;

// In the 5 positions before the D, must be two s's
	uint n = 0;
	for (uint i = 1; i <= 5; ++i)
		if (m_SS[Pos-i] == 's')
			++n;
	if (n < 2)
		return false;

// In the 5 positions after the D, must be two s's
	n = 0;
	for (uint i = 1; i <= 5; ++i)
		if (m_SS[Pos+i] == 's')
			++n;
	if (n < 2)
		return false;

	return true;
	}

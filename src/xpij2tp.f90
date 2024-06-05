! This is a Fortran90 wrapper of the Fortran77 code of
!      A. Almasy, S. Moch and A. Vogt, 1107.2263
! =====================================================================
! ..File: xpij2pt.f 
!
!                           __
! ..The parametrized 3-loop MS singlet splitting functions  P^(2)T 
!    for the evolution of unpolarized singlet fragmentation densities 
!    at mu_r = mu_f.  The expansion parameter is alpha_s/(4 pi).

! ..The distributions (in the mathematical sense) are given as in eq.
!   (B.27) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
!   The name-endings A, B, and C of the functions below correspond to
!   the kernel superscripts [2], [3], and [1] in that equation.
!
! ..The relative accuracy of these parametrisations, as well as of
!    the convolution results, is better than one part in thousand.

! ..The coefficients of 1/(1-x)_+, (ln x)/x and 1/x are exact (up
!    to a truncation of irrational coefficients).  Furthermore all
!    coefficients written as fractions (e.g., 160./27.D0) are exact.
!    The other terms at x < 1 have fitted to the exact results for x
!    between 10^-6 and 1 - 10^-6.  The coefficient of delta(1-x) of
!    P_gg^(2) have been very slightly adjusted using the low moments.
!
! ..References: S. Moch and A. Vogt,
!               Phys. Lett. B659 (2008) 290, arXiv:0708.3899  
!               A. Almasy, S. Moch and A. Vogt,
!               arXiv:1107.nnnn                            
!                           
module xpij2tp
  public::p2psta,  P2NSMTA, P2QGTA, P2GQTA, P2QQTA, P2GGTA, P2GGTB, P2GGTC
  ! =====================================================================
  !
  ! ..The (regular) pure-singlet splitting functions P_ps^(2)T.
  !    P_qq^(2)T is obtained by adding the non-singlet quantity P_NS+^(2)T
  !    A parametrization of the latter is provided in the file  xpns2pt.f.
  !
  double precision FUNCTION P2PSTA (Y, NF)
    !
    IMPLICIT REAL*8 (A-Z)
    INTEGER NF

    DL  = log(Y)
    DL1 = log(1.-Y)

    P2PST1 = - 256./(9.D0*Y)* DL**3 - 128./(9.D0*Y)* DL**2 &
         + 324.07/Y* DL + 479.87/Y&
         - 5.926* DL1**3 - 9.751* DL1**2 - 8.65* DL1 - 106.65&
         - 848.97* Y + 368.79* Y**2 - 61.284* Y**3&
         + 96.171* DL*DL1 + 656.49* DL + 425.14* DL**2 &
         + 47.322* DL**3 + 9.072* DL**4
    P2PST2 = - 128./(81.D0*Y) + 1.778* DL1**2 + 16.611* DL1 + 87.795&
         - 57.688* Y - 41.827* Y**2 + 25.628* Y**3 - 7.9934* Y**4&
         - 2.1031* DL*DL1 + 57.713* DL + 9.1682* DL**2 &
         - 1.9* DL**3 + 0.019122* DL**4&
         + 26.294* Y*DL - 7.8645* Y*DL**3

    P2PSTA = (1.-Y) * NF * ( P2PST1 + NF * P2PST2 )

  END FUNCTION P2PSTA

  ! ---------------------------------------------------------------------
  !
  !
  ! ..The quark-gluon splitting functions P_qg^(2)T.
  !   The nf^3 part is exact up to a truncation of zeta_2
  !
  double precision FUNCTION P2QGTA (Y, NF)

    IMPLICIT REAL*8 (A-Z)
    INTEGER  NF

    DL  = log(Y)
    DL1 = log(1.-Y)


    P2QG1 = - 64./Y* DL**3 - 64./Y* DL**2 + 675.83/Y* DL + 1141.7/Y&
         + 100./27.D0* DL1**4 + 350./9.D0* DL1**3&
         + 263.07* DL1**2 + 693.84* DL1 + 603.71&
         - 882.48* Y + 4723.2* Y**2 - 4745.8* Y**3 - 175.28* Y**4&
         + 1864.* DL + 1512.* DL**2 + 361.28* DL**3 &
         + 42.328* DL**4 - 1809.4* DL*DL1 - 107.59* Y*DL*DL1 &
         - 885.5* Y*DL**4
    P2QG2 = - 32./(27.D0*Y)* DL**2 - 3.1752/Y* DL - 2.8986/Y&
         - 100./27.D0* DL1**3 - 35.446* DL1**2 - 103.609* DL1&
         - 113.81 + 341.26* Y - 853.35* Y**2 + 492.1* Y**3 &
         + 14.803* Y**4 + 619.75* DL + 255.62* DL**2 &
         + 21.569* DL**3 + 966.96* DL*DL1 - 1.593*DL*DL1**2 &
         - 333.8* Y*DL**3 - 709.1* Y*DL*DL1 
    P2QG3 = 4./9.D0* (4. + 6.* (DL + DL1)&
         + (1. - 2.*Y + 2.*Y**2) * ( 3.8696 + 4.* (DL + DL1) &
         + 3.* (DL + DL1)**2 ) )

    P2QGTA = NF * ( P2QG1 + NF * P2QG2 + NF**2 * P2QG3 )

  END FUNCTION P2QGTA
  !
  ! ---------------------------------------------------------------------
  !
  !
  ! ..The gluon-quark splitting functions P_gq^(2)T
  !
  double precision  FUNCTION P2GQTA (Y, NF)

    IMPLICIT REAL*8 (A-Z)
    INTEGER NF

    DL  = log(Y)
    DL1 = log(1.-Y)

    P2GQ0 =  400./81.D0* DL1**4 + 520./27.D0* DL1**3 &
         - 220.13* DL1**2 - 152.60* DL1 + 272.85 - 7188.7* Y &
         + 5693.2* Y**2 + 146.98* Y**3 + 128.19* Y**4&
         - 30.062* DL**4 - 126.38* DL**3 - 0.71252* DL**2 &
         + 4.4136* DL - 1300.6* DL*DL1 - 71.23* DL*DL1**2 &
         + 543.8* Y*DL**3 &
         + 256./Y* DL**4 + 3712./(3.D0*Y)* DL**3&
         + 1001.89/Y* DL**2 + 4776.5/Y* DL + 5803.7/Y
    P2GQ1 =  80./81.D0* DL1**3 + 1040./81.D0* DL1**2 - 16.914* DL1&
         - 871.3 + 790.13* Y - 241.23* Y**2 + 43.252* Y**3&
         - 48.600 * DL**3 - 343.1* DL**2 - 492.* DL&
         + 55.048* DL*DL1 - 4.3465* Y*DL**3&
         + 6.0041/Y + 141.93/Y* DL &
         + 2912./(27.D0*Y)* DL**2 + 1280./(81.D0*Y)*DL**3

    P2GQTA = ( P2GQ0 + NF * P2GQ1 )

  END FUNCTION P2GQTA
  !
  ! ---------------------------------------------------------------------
  !
  !
  ! ..The regular piece of the gluon-gluon splitting function P_gg^(2)
  !
  FUNCTION P2GGTA (Y, NF)

    IMPLICIT REAL*8 (A-Z)
    INTEGER NF

    DL  = log(Y)
    DL1 = log(1.-Y)

    P2GGTA0 =   576./Y* DL**4 + 3168./Y* DL**3 + 3651.1/Y* DL**2 &
         + 10233./Y* DL + 14214.4/Y - 3590.1* DL1 - 28489.&
         + 7469.* Y + 30421.* Y**2 - 53017.* Y**3 + 19556.* Y**4&
         + 191.99* DL**4 + 3281.7* DL**3 + 13528.* DL**2 &
         + 12258.* DL - 186.4* DL*DL1 - 21328.* DL**2*DL1 &
         + 5685.8* Y*DL**3
    P2GGTA1 = + 448./(9.D0*Y)* DL**3 + 2368./(9.D0*Y)* DL**2 &
         - 5.47/Y* DL - 804.13/Y + 248.95 + 319.97* DL1&
         + 260.6* Y + 272.79* Y**2 + 2133.2* Y**3 - 926.87* Y**4&
         + 4.9934* DL + 482.94* DL**2 + 155.10* DL**3 &
         + 18.085* DL**4 + 485.18* Y*DL**3 + 1266.5* DL*DL1 &
         - 29.709* DL**2*DL1 + 87.771* DL*DL1**2
    P2GGTA2 =   32./(27.D0*Y)* DL**2 + 368./(81.D0*Y)* DL &
         + 472./(243.D0*Y) &
         - 77.190 + 153.27* Y - 106.03* Y**2 + 11.995* Y**3&
         - 5.0372* DL**3 - 44.8* DL**2 - 69.712* DL&
         - 115.01* DL*DL1 + 96.522* Y*DL*DL1 - 62.908* DL**2*DL1

    P2GGTA = P2GGTA0 + NF * ( P2GGTA1 + NF * P2GGTA2 )

  END FUNCTION P2GGTA
  !
  ! ---------------------------------------------------------------------
  !
  !
  ! ..The singular piece of the gluon-gluon splitting function P_gg^(2)T
  !    (identical to the spacelike case)
  !
  double precision FUNCTION P2GGTB (Y, NF)

    IMPLICIT REAL*8 (A-Z)
    INTEGER NF

    P2GGTB = (2643.521 - NF * 412.172 - NF**2 * 16./9.D0) / ( 1.-Y)

  END FUNCTION P2GGTB
  !
  ! ---------------------------------------------------------------------
  !
  !
  ! ..The 'local' piece of the gluon-gluon splitting function P_gg^(2)
  !   (as in the spacelike case, up to the smaller delta(1-x) adjustment)
  !
  double precision FUNCTION P2GGTC (Y, NF)

    IMPLICIT REAL*8 (A-Z)
    INTEGER NF

    DL1 = log(1.-Y)

    P2GGTC =       2643.521 * DL1 + 4425.448 + 0.003&
         - NF * ( 412.172 * DL1 +  528.720 - 0.001 )&
         + NF**2 * ( - 16./9.D0 * DL1 + 6.4630 - 0.0002)

  END FUNCTION P2GGTC

  ! =====================================================================
end module xpij2tp

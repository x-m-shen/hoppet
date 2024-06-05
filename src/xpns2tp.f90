! This is a Fortran90 wrapper of the Fortran77 code of
!      A. Almasy, S. Moch and A. Vogt, 1107.2263
! =====================================================================
! ..File: xpns2tp.f 
!
!                           __
! ..The parametrized 3-loop MS non-singlet splitting functions P^(2)T
!    for the evolution of unpolarized fragmentation densities at 
!    mu_r = mu_f.  The expansion parameter is alpha_s/(4 pi).
!
! ..The distributions (in the mathematical sense) are given as in eq.
!    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
!    The name-endings A, B, and C of the functions below correspond to
!    the kernel superscripts [2], [3], and [1] in that equation.
!
!  ..The relative accuracy of these parametrizations, as well as of
!    the convolution results, is better than one part in thousand.
!
! ..References: A. Mitov, S. Moch and A. Vogt,
!               Phys. Lett. B638 (2006) 61, hep-ph/0604053 (exact res.)
!               A. Almasy, S. Moch and A. Vogt,
!               arXiv:1107.nnnn                            (this code)
!
! ..Some (parts) of these functions are identical to the spacelike
!    results of S. Moch, J. Vermaseren and A. Vogt,
!               Nucl. Phys. B688 (2004) 101, hep-ph/040319
!      
! =====================================================================
module xpns2tp
  public::p2nspta,  P2NSMTA, P2NSTB, P2NSPTC, P2NSMTC, P2NSSTA
  !
  !
  ! ..This is the regular piece of P2_NS+.  The rational coefficients are
  !    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
  !    The N_f^2 part is exact is identical to the spacelike case.
  !
  double precision FUNCTION P2NSPTA (Y, NF)
    IMPLICIT REAL*8 (A - Z)
    INTEGER NF

    DL  = log(Y)
    DL1 = log(1.-Y)
    D81 = 1./81.D0

    P2NSPTA =   1658.7 - 707.67* DL1 + 1327.5* DL - 56.907* DL*DL1 &
         - 189.37* DL**2 - 519.37* DL1*DL**2 - 352./9.D0* DL**3 &
         + 128./81.D0* DL**4 - 4249.4* Y - 559.1* DL1*DL*Y &
         - 1075.3* Y**2 + 593.9* Y**3&
         + NF * (64./27.D0* DL**3 - 176./81.D0* DL**2 - 168.89* DL&
         - 198.10 + 466.29* Y + 181.18* Y**2 - 31.84* Y**3&
         + 5120./81.D0* DL1 - 50.758* DL*DL1 + 28.551* DL**2*DL1&
         - 39.113* Y*DL + 85.72* Y*DL*DL1 - 23.102* Y*DL**2*DL1)&
         + NF**2 * ( 32.* Y*DL/(1.-Y) * (3.* DL + 10.) + 64.&
         + (48.* DL**2 + 352.* DL + 384.) * (1.-Y) ) * D81

  END FUNCTION P2NSPTA
  !
  ! ---------------------------------------------------------------------
  !
  !
  ! ..This is the regular piece of P2_NS-.  The rational coefficients are 
  !    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
  !    The N_f^2 part is exact (and identical to that of P2_NS+). 
  !
  double precision FUNCTION P2NSMTA (Y, NF)
    IMPLICIT REAL*8 (A - Z)
    INTEGER NF

    DL  = log(Y)
    DL1 = log(1.-Y)
    D81 = 1./81.D0

    P2NSMTA =  - 140./81.D0* DL**4 - 1024./27.D0* DL**3 &
         - 38.298* DL**2 + 1625.5* DL - 707.94* DL1 + 1981.3&
         - 4885.7* Y - 577.42* Y**2 + 407.89* Y**3&
         + 1905.4* DL**2*DL1 + 1969.5* Y*DL**2*DL1 &
         + 4563.2* DL*DL1 - 34.683* Y*DL**4 &
         - 5140.6* Y*DL*DL1 - 437.03* Y*DL**3&
         + NF * ( 128./81.D0* DL**3 - 784./81.D0* DL**2 &
         - 188.99* DL - 217.84 + 511.92* Y + 209.19* Y**2 &
         - 85.786* Y**3 + 5120./81.D0* DL1 + 71.428* DL*DL1 &
         + 30.554* DL**2*DL1 + 92.453* Y*DL - 23.722* Y*DL*DL1 &
         - 18.975* Y*DL**2*DL1 )&
         + NF**2 * ( 32.* Y*DL/(1.-Y) * (3.* DL + 10.) + 64.&
         + (48.* DL**2 + 352.* DL + 384.) * (1.-Y) ) * D81

  END FUNCTION P2NSMTA
  !
  !
  ! ---------------------------------------------------------------------
  !
  !
  ! ..This is the singular piece of both P2_NS+ and P2_NS-. It is exact
  !    up to the truncation of the irrational coefficients, and identical
  !    to the spacelike result.
  !
  double precision FUNCTION P2NSTB (Y, NF)
    IMPLICIT REAL*8 (A-Z)
    INTEGER NF

    P2NSTB = ( 1174.898 - NF * 183.187 - NF**2 * 64./81.D0 ) / (1.-Y)

    RETURN
  END FUNCTION P2NSTB
  !
  ! ---------------------------------------------------------------------
  !
  !
  ! ..This is the 'local' piece of P2_NS+. Coefficients of delta(1-x) have 
  !    been shifted (minimally) relative to the exact (truncated) values.
  !
  double precision FUNCTION P2NSPTC (Y, NF)
    IMPLICIT REAL*8 (A - Z)
    INTEGER NF

    DL1 = log(1.-Y)

    P2NSPTC =       1174.898 * DL1 + 1295.624 + 0.001&
         - NF * ( 183.187 * DL1 + 173.938 - 0.003)&
         + NF**2 * ( - 64./81.D0 * DL1 + 1.13067 )

    RETURN
  END FUNCTION P2NSPTC
  !
  !
  ! ---------------------------------------------------------------------
  !
  ! ..This is the 'local' piece of P2_NS-. Coefficients of delta(1-x) have
  !    been shifted (minimally) relative to the exact (truncated) values.
  !
  double precision FUNCTION P2NSMTC (Y, NF)
    IMPLICIT REAL*8 (A - Z)
    INTEGER NF

    DL1 = log(1.-Y)

    P2NSMTC =       1174.898 * DL1 + 1295.624 - 0.002&
         - NF * ( 183.187 * DL1 + 173.938 - 0.0004 )&
         + NF**2 * ( - 64./81.D0 * DL1 + 1.13067 )

    RETURN
  END FUNCTION P2NSMTC
  !
  !
  ! ---------------------------------------------------------------------
  !
  !
  ! ..This is P2_NSS, the difference of P2_NSV and P2_NS-. It is identical
  !    to the spacelike result.
  !
  double precision FUNCTION P2NSSTA (Y, NF)

    IMPLICIT REAL*8 (A-Z)
    INTEGER NF

    D27 = 1./27.D0
    DL  = log(Y)
    Y1  = 1.- Y
    DL1 = log(Y1)

    P2NSSA = Y1* ( 151.49 + 44.51 * Y - 43.12 * Y**2 + 4.820 * Y**3 )&
         + 40.*D27 * DL**4 - 80.*D27 * DL**3 + 6.892 * DL**2&
         + 178.04 * DL + DL*DL1 * ( - 173.1 + 46.18 * DL )&
         + Y1*DL1 * ( - 163.9 / Y - 7.208 * Y )

    P2NSSTA  = NF * P2NSSA

    RETURN
  END FUNCTION P2NSSTA

end module xpns2tp
  ! =====================================================================

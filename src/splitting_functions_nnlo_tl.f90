! This is a HOPPET interface for
! - NNLO timelike splitting functions from A. Almasy, S. Moch and A. Vogt, 1107.2263, and
! - corrections to P2T_qg from H. Chen, T.Z. Yang, H.X. Zhu, Y.J Zhu 2006.10534
! <xmshen137@gmail.com>
module splitting_functions_nnlo_tl
  use types; use consts_dp
  use convolution_communicator
  use coefficient_functions; use qcd; use warnings_and_errors

  use xpns2tp
  use xpij2tp

 IMPLICIT NONE

  public :: P2Tnsp, P2Tnsm, P2Tnsv
  public :: P2Tqq, P2Tgg, P2Tqg, P2Tgq, P2Tgq_CYZZ

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NS_plus
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function P2Tnsp_Regular(x, nf) result(res)
    implicit none
    real(8), intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = p2nspta(x, nf)

  end function P2Tnsp_Regular

  function P2Tnsp_Singular(x, nf) result(res)
    implicit none
    real(8), intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = (1174.898d0 - nf * 183.187d0 - nf * nf * 64.0d0 / 81.0d0) / (1.0d0 - x)
  end function P2Tnsp_Singular

  function P2Tnsp_Delta(nf) result(res)
    implicit none
    integer, intent(in) :: nf
    real(8) :: res
    real(8) :: dl1

    dl1 = 0.0d0

    res = 1174.898d0 * dl1 + 1295.624d0 + 0.001d0 &
          - nf * (183.187d0 * dl1 + 173.938d0 - 0.003d0) &
          + nf * nf * (-64.0d0 / 81.0d0 * dl1 + 1.13067d0)
  end function P2Tnsp_Delta
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NS_minus
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function P2Tnsm_Regular(x, nf) result(res)
    implicit none
    real(8), intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = p2nsmta(x, nf)
  end function P2Tnsm_Regular

  function P2Tnsm_Singular(x, nf) result(res)
    implicit none
    real(8), intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = (1174.898d0 - nf * 183.187d0 - nf * nf * 64.0d0 / 81.0d0) / (1.0d0 - x)
  end function P2Tnsm_Singular

  function P2Tnsm_Delta(nf) result(res)
    implicit none
    integer, intent(in) :: nf
    real(8) :: res
    real(8) :: dl1

    dl1 = 0.0d0

    res = 1174.898d0 * dl1 + 1295.624d0 - 0.002d0 &
          - nf * (183.187d0 * dl1 + 173.938d0 - 0.0004d0) &
          + nf * nf * (-64.0d0 / 81.0d0 * dl1 + 1.13067d0)
  end function P2Tnsm_Delta
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NSS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function P2Tnss_Regular(x, nf) result(res)
    implicit none
    real(8), intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = p2nssta(x, nf)

  end function P2Tnss_Regular

!-----------------------------------------------------------------------------
! the singlet splitting functions
!-----------------------------------------------------------------------------
  function P2Tps_Regular(x, nf) result(res)
    implicit none
    real(8), intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = p2psta(x, nf)
  end function P2Tps_Regular
!-----------------------------------------------------------------------------
  function P2Tqg_Regular(x, nf) result(res)
    implicit none
    double precision, intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = 2 * nf * p2gqta(x, nf) ! swap q,g
  end function P2Tqg_Regular
!-----------------------------------------------------------------------------
  function P2Tgq_Regular(x, nf) result(res)
    implicit none
    double precision, intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = p2qgta(x, nf) / 2d0 / nf ! swap q,g

  end function P2Tgq_Regular

  function P2Tgq_Regular_CYZZ(x, nf) result(res)
    implicit none
    double precision, intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = p2qgta(x, nf) / 2d0 / nf ! swap q,g
    ! ref: H. Chen, T.Z. Yang, H.X. Zhu, Y.J Zhu 2006.10534
    res = res + pi*pi/3.0*(4.0/3.0 - 3.0)*(11.0 - 2.0*nf/3.0)*(-4.0 + 8*x + x*x + 6*log(x)*(1-2*x+2*x*x))
  end function P2Tgq_Regular_CYZZ
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function P2Tgg_Regular(x, nf) result(res)
    implicit none
    double precision, intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = p2ggta(x, nf)

  end function P2Tgg_Regular

  function P2Tgg_Singular(x, nf) result(res)
    implicit none
    double precision, intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res

    res = ( 2643.521 - nf * 412.172 - nf * nf * 16. / 9. ) / ( 1 - x )
  end function P2Tgg_Singular

  function P2Tgg_Local(x, nf) result(res)
    implicit none
    double precision, intent(in) :: x
    integer, intent(in) :: nf
    real(8) :: res
    double precision :: dl1

    dl1 = log(1 - x)
    res = 2643.521 * dl1 + 4425.448 + 0.003 - nf * ( 412.172 * dl1 + 528.720 - 0.001 ) + &
          nf * nf * ( - 16. / 9. * dl1 + 6.4630 - 0.0002 )
  end function P2Tgg_Local

  function P2Tgg_Delta(nf) result(res)
    implicit none
    integer, intent(in) :: nf
    real(8) :: res
    double precision :: dl1

    dl1 = 0.0
    res = 2643.521 * dl1 + 4425.448 + 0.003 - &
         nf * ( 412.172 * dl1 + 528.720 - 0.001 ) + &
         nf * nf * ( - 16. / 9. * dl1 + 6.4630 - 0.0002 )
  end function P2Tgg_Delta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! public
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function P2Tnsp(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 0.125_dp * P2Tnsp_regular(x, nf_int) + 0.125_dp * P2Tnsp_singular(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - 0.125_dp * P2Tnsp_singular(x, nf_int)
    case(cc_DELTA)
       res = 0.125_dp * P2Tnsp_delta(nf_int)
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function P2Tnsp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function P2Tnsm(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 0.125_dp * P2Tnsm_regular(x, nf_int) + 0.125_dp * P2Tnsm_singular(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - 0.125_dp * P2Tnsm_singular(x, nf_int)
    case(cc_DELTA)
       res = 0.125_dp * P2Tnsm_delta(nf_int)
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function P2Tnsm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function P2Tnss(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 0.125_dp * P2Tnss_regular(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function P2Tnss

  function P2Tnsv(y) result(res)
    implicit none
    real(dp), intent(in) :: y
    real(dp)             :: res

    res = P2Tnsm(y) + P2Tnss(y)
  end function P2Tnsv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! singlet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function P2Tps(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 0.125_dp * P2Tps_regular(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function P2Tps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function P2Tqq(y) result(res)
    implicit none
    real(dp), intent(in) :: y
    real(dp)             :: res

    res = P2Tnsp(y) + P2Tps(y)
  end function P2Tqq
  !----------------------------------------------------------------------
  function P2Tgg(y) result(res)
    implicit none
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 0.125_dp * P2Tgg_regular(x, nf_int) + 0.125_dp * P2Tgg_singular(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - 0.125_dp * P2Tgg_singular(x, nf_int)
    case(cc_DELTA)
       res = 0.125_dp * P2Tgg_delta(nf_int)
    end select

    if (cc_piece /= cc_DELTA) res = res * x

 end function P2Tgg
  !----------------------------------------------------------------------
  function P2Tqg(y) result(res)
    implicit none
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 0.125_dp * P2Tqg_regular(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function P2Tqg
  !----------------------------------------------------------------------
  function P2Tgq(y) result(res)
    implicit none
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 0.125_dp * P2Tgq_regular(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function P2Tgq

  function P2Tgq_CYZZ(y) result(res)
    implicit none
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 0.125_dp * P2Tgq_regular_CYZZ(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function P2Tgq_CYZZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module SPLITTING_FUNCTIONS_NNLO_TL

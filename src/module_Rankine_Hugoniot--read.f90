!------------------------------------------------------------------------------
!MODULE: Rankine_Hugoniot
!------------------------------------------------------------------------------
!> \brief Provide Rankine-Hugoniot relations for a multi-especies gas mixture.
!!
!! Resolution of pre- and post-shock conditions for an oblique shock
!! wave, with a angle of incidence \f$ \beta \f$ from the ground, and
!! spreads left (index 1) to right (index 2) at an average speed \f$
!! v_1 \f$. In the reference frame associated with shock, in STEADY
!! state, we calculate the flow thermodynamic state after the shock
!! passage \f$ (p_2, T_2, v_2) \f$ and the angle of deflection of the
!! velocity after impact \f$ \theta \f$.
!!
!! To resolve Rankine-Hugoniot's equations, we use the
!! Newton-Raphson's iterative method with the initial value \f$
!! \phi=T_2/T_1 \f$, calculated with Rankine-Hugoniot's relation for
!! an ideal gas and thermodynamics properties independent to
!! temperature. If this method fails, we then use the dichotomy
!! method.
!!
!!                          shock
!!                            ||
!!          <----             ||                <----
!!           (2)              ||                 (1)
!!      u2, T2, P2, rho2      ||           u1, T1, P1, rho1
!!
!! References:
!!
!! -# R. E. Mitchell, R. J. Kee. _A general-purpose computer code for
!!    predicting chemical kinetic behavior behind incident and
!!    reflected shocks_. Sandia Report, p. 13, (1982).
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
module Rankine_Hugoniot

  use parameters
  use adim

  implicit none

  integer (ip) , parameter , private        :: itmax   =  50        ! maximum number of iterations for the method
  real (dp)    , parameter , private        :: eps_phi =  1.0e-6_dp ! solution error tolerance
  real (dp)    , parameter , private        :: dphi    =  1.0e-5_dp ! small epsilon to add to the current value of phi
  real (dp)    , parameter , private        :: neg_tol = -1.0e-2_dp ! maximum negative tolerance to avoid square of negative numbers


contains


!> \brief Selector to update the boundary. This subroutine reads the keyword to select the corresponding boundary condition.

  subroutine RH_variables ( adi , p1 , T1 , rho1 , v1 , beta_ , &
                                  p2 , T2 , rho2 , v2 , theta_ )


    type (adi_type)  , intent (in)      :: adi    !< non-dimensional derived type
    real (dp)        , intent (in)      :: p1     !< pressure before shock
    real (dp)        , intent (in)      :: T1     !< temperature before shock
    real (dp)        , intent (in)      :: v1     !< velocity before shock
    real (dp)        , intent (in)      :: beta_  !< angle between wall and the shock
    real (dp)        , intent (inout)   :: rho1   !< density before shock
    real (dp)        , intent (inout)   :: p2     !< pressure after shock
    real (dp)        , intent (inout)   :: T2     !< temperature after shock
    real (dp)        , intent (inout)   :: v2     !< velocity after shock
    real (dp)        , intent (inout)   :: rho2   !< density after shock
    real (dp)        , intent (inout)   :: theta_ !< angle between wall and velocity after shock


    real (dp) :: u1 , cs1 , M1n
    real (dp) :: u2
    real (dp) :: phi , psi


    ! Initilize the state 1

    rho1   = P1 / T1
    u1     = v1 * sin (beta_)            ! normal component of velocity relative to the shock
    cs1    = sqrt ( adi % gamma * p1 / rho1 )  ! speed of sound
    M1n    = u1 / cs1                    ! normal Mach number

    phi    = RH_phi ( adi % gamma , M1n )      ! temperature ratio
    psi    = RH_psi ( adi % gamma , M1n )      ! pressure ration

    
    ! calculate conditions at state 2

    T2     = T1 * phi
    p2     = p1 * psi
    rho2   = P2 / T2
    u2     = rho1 * u1 / rho2                       ! normal component of velocity relative to the shock
    theta_ = beta_ - atan ( u2 * tan (beta_) / u1 )
    v2     = u2 / sin ( beta_ - theta_ )            ! velocity magnitude
    


  end subroutine RH_variables


!> \brief Calculate the initial value of the temperature ratio. This uses the Rankine-Hugoniot's relations for a pure air flow assuming that thermodynamics properties are independant to temperature.

  function RH_phi ( gamma1 , M1n )


    real (dp) , intent (in)  :: gamma1  !< ratio of specific heats
    real (dp) , intent (in)  :: M1n     !< normal component of the Mach relative to the shock
    
    real (dp)                :: RH_phi !< temperature ratio
    real (dp)                :: wrk


    wrk = 0.5_dp * ( gamma1 + 1.0_dp ) * M1n

    RH_phi = gamma1 * M1n * M1n - 0.5_dp * ( gamma1 - 1.0_dp )
    RH_phi = RH_phi * ( 1.0_dp + 0.5_dp * ( gamma1 - 1.0_dp ) * M1n * M1n )
    RH_phi = RH_phi / ( wrk * wrk )


  end function RH_phi


!> \brief Calculate the pressure ratio.

  function RH_psi ( gamma1 , M1n )

    real (dp) , intent (in)  :: gamma1  !< ratio of specific heats
    real (dp) , intent (in)  :: M1n     !< normal component of the Mach relative to the shock

    real (dp)                :: RH_psi  !< temperature ratio


    RH_psi = 2.0_dp * gamma1 * M1n * M1n - ( gamma1 - 1.0_dp )
    RH_psi = RH_psi / ( gamma1 + 1.0_dp )
   

  end function RH_psi


end module Rankine_Hugoniot

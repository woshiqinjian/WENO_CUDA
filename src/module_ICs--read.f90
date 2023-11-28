!------------------------------------------------------------------------------
! MODULE: BCs
!------------------------------------------------------------------------------
!> \brief Initial conditions.
!!
!! This module contains the definition of all the boundary conditions.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
module ICs

  use parameters
  use input
  use adim
  use BCs
  use tools
  use deriv

  implicit none


contains


!> \brief Selector to set the initial condition. This subroutine reads the keyword to select the corresponding initial condition.

  subroutine init_selector ( inp , adi , grid , T , v )


    type (inp_type)                               , intent (in)    :: inp  !< input derived type
    type (adi_type)                               , intent (inout) :: adi  !< non-dimensional derived type
    type (inp_grid)                               , intent (in)    :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

   !  real (dp) :: f_x


    if ( inp % init == init_rank ) then ! rank
       call ini_rank ( adi , T , v )
    else if ( inp % init == Ambient ) then ! ambient
       call ini_amb ( inp , adi , grid % x , grid % y , grid % z , T , v )
    else if ( inp % init == Vortex ) then ! vortex
       call ini_vortex ( inp , adi , grid % x , grid % y , grid % z , T , v)
    else if ( inp % init == vortexring ) then ! vortex ring
       call ini_vortexring ( inp , adi , grid % x , grid % y , grid % z , T , v )
    else if ( inp % init == shock ) then ! shock
       call ini_shock ( inp , adi , grid % x , grid % y , grid % z , T , v)
    else if ( inp % init == DMR ) then ! double Mach reflection
       call ini_DMR ( adi , grid % x , grid % y , T , v )
    else if ( inp % init == HIT ) then ! decaying HIT
       call ini_HIT ( inp , adi , grid , T , v )
    else if ( inp % init == Bogey ) then ! shear layer
       call ini_shear_layer_bogey ( inp , adi , grid , T , v )       
    else if ( inp % init == shear ) then ! shear layer
       call ini_shear_layer ( inp , adi , grid , T , v )
    else if ( inp % init == channel ) then
       call ini_channel_flow ( inp , adi , grid , T , v )
    ! else if ( inp % init == Diffusion ) then
    !    call ini_diffusion ( adi , thd , grid % x , T , W_i , cp , ha , v )
    ! else if ( inp % init == Fu ) then ! shear layer
    !    call ini_shear_layer_fu ( adi , thd , grid % y , T , W_i , cp , ha , v )
    ! else if ( inp % init == Miller ) then ! shear layer
    !    call ini_shear_layer_miller ( adi , thd , grid % y , T , W_i , cp , ha , v )
    ! else if ( inp % init == MIXChengMiller ) then ! shear layer
    !    call ini_shear_layer_mixchengmiller ( inp , adi , thd , grid % y , T , W_i , cp , ha , v )
    ! else if ( inp % init == Premix ) then
    !    call ini_premix ( adi , thd , T , W_i , cp , ha , v )
    ! else if ( inp % init == Jet ) then
    !    call ini_jet ( adi , thd , T , W_i , cp , ha , v )
    ! else if ( inp % init == DIFFMIXChengMiller ) then
    !    call ini_diffusion_mixchengmiller ( inp , adi , thd , grid % x , T , W_i , cp , ha , v )
    ! else if ( inp % init == Fedkiw ) then ! shock tube
    !    call ini_fedkiw ( adi , thd , grid % x , T , W_i , cp , ha , v )
    ! else if ( inp % init == DMR ) then ! double Mach reflection
    !    call ini_DMR ( adi , thd , grid % x , grid % y , T , W_i , cp , ha , v )
    ! else if ( inp % init == Poiseuille ) then ! poiseuille
    !    call ini_poiseuille ( adi , thd , grid % x , grid % y , grid % z , T , W_i , cp , ha , v )
    ! else if ( inp % init == Vortex ) then ! vortex
    !    call ini_vortex ( adi , thd , grid % x , grid % y , T , W_i , cp , ha , v )
    ! else if ( inp % init == Detonation ) then ! Detonation
    !    call ini_detonation ( inp , adi , thd , grid % x , grid % y , grid % z , T , W_i , cp , ha , v )
    else
       call end_cuda ('initalization ' // trim (inp % init) // ' not defined')
    end if


  end subroutine init_selector


  subroutine ini_rank ( adi , T , v )


    type (adi_type)                               , intent (in)    :: adi
    real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v


    integer (ip)  :: ok , i , j , k , l
    real    (dp)  :: P , Temp , rho


    allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
               T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate ini')    


    !Initialisation des matrices
    v = 1.0_dp
    T = 1.0_dp

    v(sx:ex,sy:ey,sz:ez,2) = 5.0_dp  !real (  rank , dp ) + 1.0_dp
    v(sx:ex,sy:ey,sz:ez,3) = 10.0_dp !real ( -rank , dp ) - 1.0_dp
    v(sx:ex,sy:ey,sz:ez,4) = 0.0_dp
    v(sx:ex,sy:ey,sz:ez,5) = ( 1.0_dp / ( adi % gamma - 1.0_dp ) ) &
            + 0.5_dp * v(sx:ex,sy:ey,sz:ez,1) * (5.0_dp * 5.0_dp + 10.0_dp * 10.0_dp )
    
    
  end subroutine ini_rank

!> \brief Ambient IC.0 The field is at rest and in ambient conditions (U=0, T=300K, P=1atm, in air).

  subroutine ini_amb ( inp , adi , x , y , z , T , v )


    type (inp_type)                               , intent (in)    :: inp  !< input derived type
    type (adi_type)                               , intent (in)    :: adi  !< non-dimensional derived type
    real (dp) , allocatable , dimension (:)       , intent (in)    :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

    integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
    integer (ip)                                           :: ok , i , j , k , l

    real (dp)  :: ux, vy, wz, rho , Press , Ma , cs


    allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
               T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate ini_amb')


    v   = 0.0_dp
    T   = 0.0_dp

    
    ux = 0.0_dp
    vy = 0.0_dp
    wz = 0.0_dp
    
    Rho = 1.0_dp
    Ma  = 0.1_dp
    cs  = 1.0_dp / Ma
    

    ! Application of (T, p, rho, u, v, w) in the calculation domain
    do k = sz , ez ! z-direction
       do j = sy , ey ! y-direction
          do i = sx , ex ! x-direction
             
             T(i,j,k) = cs * cs / adi % gamma
             Press    = rho * T(i,j,k)

             v (i,j,k,1) = rho
             v (i,j,k,2) = rho * ux
             v (i,j,k,3) = rho * vy
             v (i,j,k,4) = rho * wz
             v (i,j,k,5) = Press / ( adi % gamma - 1.0_dp ) + 0.5_dp * rho * (ux*ux + vy*vy + wz*wz)

             
          end do
       end do
    end do

    
    call upd_prim_var_domain ( adi , T , v )
    call comm_cons (v)
    call upd_extrapolation(v)
    call upd_prim_var_ghost ( adi , T , v )
    

  end subroutine ini_amb

!> \brief Vortex advection IC.

  subroutine ini_vortex ( inp , adi , x , y , z , T , v)
    
    type (inp_type)                               , intent (in)    :: inp  !< input derived type
    type (adi_type)                               , intent (in)    :: adi  !< non-dimensional derived type
    real (dp) , allocatable , dimension (:)       , intent (in)    :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

    integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
    integer (ip)                                           :: ok , i , j , k


    real    (dp) :: x_0 , y_0 , circ , r2 , wrk1 , wrk2 , wrk3 , wrk4
    real    (dp) :: lamda0 , rho , ux , vy , wz , Press , eint , u0 , v0

    allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
               T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate ini_vortex')
    
    
    x_0 = 5.0_dp
    y_0 = 5.0_dp

    u0 = 1.0_dp
    v0 = 0.0_dp

    circ = 5.0_dp
    
    wrk1  = circ / (pi+pi)
    wrk2  = - ( (adi % gamma - 1.0_dp) * circ * circ ) / ( 8.0_dp * adi % gamma * pi * pi )

    wz = 0.0_dp
    lamda0 = 1.0_dp
         
    do k = sz , ez 
       do j = sy , ey
          do i = sx , ex
             
             r2 = ( x(i) - x_0 ) * ( x(i) - x_0 ) + &
                  ( y(j) - y_0 ) * ( y(j) - y_0 )

             wrk3 = exp ( 0.5_dp * ( 1.0_dp - r2 ) ) ! this requires r2 to be larger than one (i.e. physical)
             wrk4 = exp ( 1.0_dp - r2 )

             T(i,j,k) = 1.0_dp + wrk2 * wrk4

             rho = T(i,j,k)**(1.0_dp/(adi % gamma-1.0_dp))

             ux  = u0 - wrk1 * wrk3 * ( y (j) - y_0 )
             vy  = v0 + wrk1 * wrk3 * ( x (i) - x_0 )

             Press = rho * T(i,j,k)

             
             v(i,j,k,1) = rho
             v(i,j,k,2) = rho * ux
             v(i,j,k,3) = rho * vy
             v(i,j,k,4) = rho * wz
             v(i,j,k,5) = Press / (adi % gamma - 1.0_dp) + 0.5_dp * rho * (ux*ux + vy*vy + wz*wz)

                          
          end do
       end do
    end do


    call upd_prim_var_domain ( adi , T , v )
    call comm_cons (v)
    call upd_extrapolation (v)
    call upd_prim_var_ghost ( adi , T , v )
        
  end subroutine ini_vortex

!> \brief Vortex advection IC.

  subroutine ini_vortexring ( inp , adi , x , y , z , T , v )
    
    type (inp_type)                               , intent (in)    :: inp  !< input derived type
    type (adi_type)                               , intent (in)    :: adi  !< non-dimensional derived type
    real (dp) , allocatable , dimension (:)       , intent (in)    :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

    integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
    integer (ip)                                           :: ok , i , j , k


    real    (dp) :: x_0 , y_0 , z_0 , circ_1 , circ_2 ,omga
    real    (dp) :: lx , lc, sinalpha , cosalpha , sintheta , costheta
    real    (dp) :: rho , Press , ux , vy , wz 

    allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
               T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               stat = ok )
               
    if ( ok > 0 ) call abort_mpi ('error allocate ini_vortexring')
    
    !< center of ring
    x_0 = 1.0_dp
    y_0 = 3.0_dp
    z_0 = 3.0_dp
    !< parameter of the ring
    circ_1 = 2.0_dp
    circ_2 = 0.1_dp
    !< flow information
    omga  = 4.0_dp * pi !< rotation of vortexring
    rho   = 1.4_dp
    Press = 1.0_dp
    
    do k = sz , ez 
       do j = sy , ey
          do i = sx , ex
             
             lx = sqrt( ( z(k) - z_0 )**2 + ( y(j) - y_0 )**2 )
             lc = sqrt( ( x(i) - x_0 )**2 + ( lx - circ_1 )**2 )

             if ( lc <= circ_2 ) then
               sinalpha = ( z(k) - z_0 )/lx
               cosalpha = ( y(j) - y_0 )/lx
               sintheta = ( x(i) - x_0 )/lc
               costheta = ( lx - circ_1 )/lc

               ux = -0.5_dp * omga * lc * costheta
               vy =  0.5_dp * omga * lc * sintheta * cosalpha
               wz =  0.5_dp * omga * lc * sintheta * sinalpha
             else
               ux = 0.0_dp
               vy = 0.0_dp
               wz = 0.0_dp
             end if
             
             v(i,j,k,1) = rho
             v(i,j,k,2) = rho * ux
             v(i,j,k,3) = rho * vy
             v(i,j,k,4) = rho * wz
             v(i,j,k,5) = Press / (adi % gamma - 1.0_dp) + 0.5_dp * rho * (ux**2 + vy**2 + wz**2)
            
          end do
       end do
    end do

    call upd_prim_var_domain ( adi , T , v )
    call comm_cons (v)
    call upd_extrapolation (v)
    call upd_prim_var_ghost ( adi , T , v )
        
  end subroutine ini_vortexring

!> \brief Vortex advection IC.

  subroutine ini_shock ( inp , adi , x , y , z , T , v)
    
    type (inp_type)                               , intent (in)    :: inp  !< input derived type
    type (adi_type)                               , intent (in)    :: adi  !< non-dimensional derived type
    real (dp) , allocatable , dimension (:)       , intent (in)    :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

    integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
    integer (ip)                                           :: ok , i , j , k

    real    (dp) :: rho , ux , vy , wz , Press , Ma , cs

    allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
               T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate ini_vortex')
    

    v = 0.0_dp
    T = 0.0_dp

    ux = 1.0_dp
    vy = 0.0_dp
    wz = 0.0_dp

    rho = 1.0_dp
    Ma  = 5.0_dp
    cs  = ux / Ma
 
         
    do k = sz , ez 
       do j = sy , ey
          do i = sx , ex
             
             T(i,j,k) = cs * cs / adi % gamma
             Press    = rho * T(i,j,k)
                          
             v(i,j,k,1) = rho
             v(i,j,k,2) = rho * ux
             v(i,j,k,3) = rho * vy
             v(i,j,k,4) = rho * wz
             v(i,j,k,5) = Press / (adi % gamma - 1.0_dp) + 0.5_dp * rho * (ux*ux + vy*vy + wz*wz)

                          
          end do
       end do
    end do


    call upd_prim_var_domain ( adi , T , v )
    call comm_cons (v)
    call upd_extrapolation (v)
    call upd_prim_var_ghost ( adi , T , v )
        
  end subroutine ini_shock


!> \brief Double Mach reflection IC.
!!
!! Initialization of the 2D test case of the double Mach reflection
!! (DMR) on a wall with an incident _unsteady_ shock wave arriving
!! bias on a horizontal plane wall.
!!
!!                          shock
!!                            ||
!!                          --||--> Us
!!           (2)              ||                 (1)
!!                            ||        u1=000, T1=300K, P1=1atm
!!
!! References :
!!
!! -# Vasilev et al. "The wall-jetting effect in Mach reflection:
!!    Navier–Stokes simulations," Journal of Fluid Mechanics.
!! -# http://ufrmeca.univ-lyon1.fr/~buffat/COURS/AERO_HTML/node49.html
!! -# A.Kourta, 'Dynamique des gaz, Supersoniques linearises et ondes
!!    instationnaires, Fascicule de cours', Polytech Orléans, 2012.
!! -# R.E.Mitchell, R.J.Kee, 'A general-purpose computer code for
!!    predicting chemical kinetic behavior behind incident and
!!    reflected shocks', Sandia Report, 1982, p.13.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  subroutine ini_DMR ( adi , x , y , T , v )


    type (adi_type)                               , intent (in)    :: adi
    real (dp) , allocatable , dimension (:)       , intent (in)    :: x
    real (dp) , allocatable , dimension (:)       , intent (in)    :: y
    real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v


    integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)

    integer (ip) :: ok , i , j , k , l

    real (dp)  :: Ms, Us ! shock Mach number & shock velocity
    real (dp)  :: u1, v1, w1, vnorme1, p1, T1, rho1, cs1   ! pre-shock variables (RIGHT)                  
    real (dp)  :: u2, v2, w2, vnorme2, p2, T2, rho2, cs2   ! post-shock variables (LEFT)

    real (dp)             :: x0      ,& ! abscisse of the shock impact
                             thetaw  ,& ! angle of the angled wedges
                             beta    ,& ! angle=(wall, shock)
                             slope   ,& ! = tan(beta)
                             slope_i ,& ! = 1/tan(beta)
                             theta      ! angle(wall, velocity magnitude after shock)

    allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
               T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate ini')

    v   = 0.0_dp
    T   = 0.0_dp

    ! !====== Parameters of the DMR problem ======
    ! ! Angled wedge
    ! x0      = 1.0_dp / 6.0_dp
    ! thetaw  = 30.0_dp * pi / 180.0_dp
    ! beta    = pi/2.0_dp - thetaw
    ! slope   = tan(beta)
    ! slope_i = 1.0_dp/slope
    
    ! Ms = 10.0_dp  ! shock Mach number
    
    ! ! Pre-shock conditions
    ! u1 = 0.0_dp
    ! v1 = 0.0_dp
    ! w1 = 0.0_dp
    
    ! vnorme1 = sqrt(u1*u1 + v1*v1 + w1*w1)
    
    ! rho1 = 1.4_dp
    ! P1   = 1.0_dp
    ! T1   = P1 / rho1
    
    ! cs1  = sqrt ( adi % gamma * p1 / rho1 )

    ! Us = Ms * cs1 !'cs1' because of wave speed is measured from the speed of sound upstream of shock
    ! vnorme1 = Us

    ! call RH_variables (adi , p1, T1, rho1, vnorme1, pi/2.0_dp, &
    !                          p2, T2, rho2, vnorme2, theta)

    ! write (*,*) vnorme1 , vnorme2
    
    ! ! Post-shock conditions
    ! u2 =  vnorme2 * sin(beta)
    ! v2 = -vnorme2 * cos(beta)
    ! w2 =  0.0_dp

    ! vnorme2 = sqrt(u2*u2 + v2*v2 + w2*w2)
    
    ! !rho2 = 8.0_dp
    ! !P2   = 116.5_dp
    ! !T2   = P2 / rho2
    
    ! cs2  = sqrt ( adi % gamma * p2 / rho2 )

    ! End of calculation of pre- and post-shock conditions

    ! Application of (T, p, rho, u, v, w) in the calculation domain
    do k = sz , ez ! z-direction
       do j = sy , ey ! y-direction
          do i = sx , ex ! x-direction

             !if ( x(i) <= (x0+y(j)*slope_i) .and. y(j) >= (slope*(x(i)-x0)) ) then

             if ( x(i) < (1.0_dp / 6.0_dp + y(j) / sqrt(3.0_dp) ) ) then
                ! POST-shock conditions (left)

                rho2      =  8.0_dp
                P2        =  116.5_dp
                T (i,j,k) =  P2 / rho2
                u2        =  8.25_dp * cos(pi / 6.0_dp)
                v2        = -8.25_dp * sin(pi / 6.0_dp)
                w2        =  0.0_dp
                                
                v (i,j,k,1) = rho2
                v (i,j,k,2) = rho2 * u2
                v (i,j,k,3) = rho2 * v2
                v (i,j,k,4) = rho2 * w2
                v (i,j,k,5) = P2 / (adi % gamma - 1.0_dp) + 0.5_dp * rho2 * (u2*u2 + v2*v2 + w2*w2)

             else
                ! PRE-shock conditions (right)

                rho1      = 1.4_dp
                P1        = 1.0_dp
                T (i,j,k) = P1 / rho1
                u1        = 0.0_dp
                v1        = 0.0_dp
                w1        = 0.0_dp                               

                v (i,j,k,1) = rho1
                v (i,j,k,2) = rho1 * u1
                v (i,j,k,3) = rho1 * v1
                v (i,j,k,4) = rho1 * w1
                v (i,j,k,5) = P1 / (adi % gamma - 1.0_dp) + 0.5_dp * rho1 * (u1*u1 + v1*v1 + w1*w1)

             end if

          end do
       end do
    end do

    call upd_prim_var_domain ( adi , T , v )
    call comm_cons (v)
    call upd_extrapolation (v)
    call upd_prim_var_ghost ( adi , T , v )
    

  end subroutine ini_DMR  

!> \brief Vortex advection IC.

  subroutine ini_HIT ( inp , adi , grid , T , v)
    
    type (inp_type)                               , intent (in)    :: inp  !< input derived type
    type (adi_type)                               , intent (inout) :: adi  !< non-dimensional derived type
    type (inp_grid)                               , intent (in)    :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

    integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)

    real (dp) , allocatable , dimension (:,:,:) :: ux , vy , wz , dux_x , dvy_y , dwz_z
    real (dp) , allocatable , dimension (:)     :: dummy , dx_i2
    
    integer (ip) :: ok , i , j , k , nx_HIT

    real    (dp) :: rho , Press , Mt , cs
    real    (dp) :: KE , kp_HIT , Re_taylor , KE_0

    real    (dp) :: dudx2 , mean_rho , u_rms , mean_ux , lamda1
    real    (dp) :: dvdy2 , mean_vy , lamda2
    real    (dp) :: dwdz2 , mean_wz , lamda3

    allocate ( v     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
               T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               ux    ( 1-ng:ntx+ng , 1-ng:nty+ng , 1-ng:ntz+ng )         , &
               vy    ( 1-ng:ntx+ng , 1-ng:nty+ng , 1-ng:ntz+ng )         , &
               wz    ( 1-ng:ntx+ng , 1-ng:nty+ng , 1-ng:ntz+ng )         , &
               dux_x ( 1-ng:ntx+ng , 1-ng:nty+ng , 1-ng:ntz+ng )         , &
               dvy_y ( 1-ng:ntx+ng , 1-ng:nty+ng , 1-ng:ntz+ng )         , &
               dwz_z ( 1-ng:ntx+ng , 1-ng:nty+ng , 1-ng:ntz+ng )         , &
               dx_i2 ( 1-ng:ntx+ng )                                     , &
               dummy ( 1-ng:ntx+ng )                                     , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate ini_vortex')
    

    v = 0.0_dp
    T = 0.0_dp

    nx_HIT    = 256_ip 
    kp_HIT    = 4.0_dp
    Mt        = 0.488_dp
    KE_0      = 0.5_dp
    Re_taylor = 175.0_dp
 
    KE       = 0.0_dp

    rho   = 1.0_dp
    cs    = sqrt(2.0_dp * KE_0) / Mt  ! sqrt(2 x KE)/Mt
    press = rho * cs * cs / adi % gamma

    open(unit_HIT , file = 'HIT_0032_4.0' , form = 'formatted')

    do k = 1 , ntz
       do j = 1 , nty
          do i = 1 , ntx

             read(unit_HIT,'(3(F25.20,2X))') ux(i,j,k) , vy(i,j,k) , wz(i,j,k)

             KE = KE + 0.5_dp * ( ux(i,j,k)*ux(i,j,k) + &
                                  vy(i,j,k)*vy(i,j,k) + &
                                  wz(i,j,k)*wz(i,j,k)   )

          end do
       end do      
    end do

    close(unit_HIT)


    do k = sz , ez 
       do j = sy , ey
          do i = sx , ex
             
             T(i,j,k) = cs * cs / adi % gamma
             Press    = rho * T(i,j,k)
                          
             v(i,j,k,1) = rho
             v(i,j,k,2) = rho * ux(i,j,k)
             v(i,j,k,3) = rho * vy(i,j,k)
             v(i,j,k,4) = rho * wz(i,j,k)
             v(i,j,k,5) = Press / (adi % gamma - 1.0_dp) + 0.5_dp * rho * ( ux(i,j,k)*ux(i,j,k) + &
                                                                            vy(i,j,k)*vy(i,j,k) + &
                                                                            wz(i,j,k)*wz(i,j,k)   )
                                      
          end do
       end do
    end do

    call upd_prim_var_domain ( adi , T , v )
    call comm_cons (v)
    !call upd_boundaries ( inp , adi , 0.0_dp , 0.0_dp , grid , T , v )
    call upd_prim_var_ghost ( adi , T , v )

    !u_rms    = 0.0_dp
    !mean_rho = 1.0_dp
    
    !dudx2    = 0.0_dp
    !dvdy2    = 0.0_dp
    !dwdz2    = 0.0_dp   
    
    !mean_ux  = 0.0_dp
    !mean_vy  = 0.0_dp
    !mean_wz  = 0.0_dp
    
    !dux_x    = 0.0_dp
    !dvy_y    = 0.0_dp
    !dwz_z    = 0.0_dp
    
    !dummy    = 1.0_dp
    
    
    !dx_i2 = 256.0_dp / ( acos(-1.0_dp) + acos(-1.0_dp) )
    !call dx_IC ( dx_i2 , ux , dux_x )
    !call dy_IC ( dx_i2 , vy , dvy_y )
    !call dz_IC ( dx_i2 , wz , dwz_z )
        
    !do k = 1 , ntz
    !   do j = 1 , nty
    !      do i = 1 , ntx

    !         mean_ux = mean_ux + ux(i,j,k)
    !         mean_vy = mean_vy + vy(i,j,k)
    !         mean_wz = mean_wz + wz(i,j,k)
                         
    !         u_rms = u_rms + ux(i,j,k)*ux(i,j,k) + &
    !                         vy(i,j,k)*vy(i,j,k) + &
    !                         wz(i,j,k)*wz(i,j,k)
             
    !         dudx2 = dudx2 + dux_x(i,j,k) * dux_x(i,j,k)
    !         dvdy2 = dvdy2 + dvy_y(i,j,k) * dvy_y(i,j,k)
    !         dwdz2 = dwdz2 + dwz_z(i,j,k) * dwz_z(i,j,k)
             
    !      end do
    !   end do
    !end do    

    !mean_ux = mean_ux / ( ntx * nty * ntz )
    !mean_vy = mean_vy / ( ntx * nty * ntz )
    !mean_wz = mean_wz / ( ntx * nty * ntz )
        
    !dudx2  = dudx2 / ( ntx * nty * ntz )
    !dvdy2  = dvdy2 / ( ntx * nty * ntz )
    !dwdz2  = dwdz2 / ( ntx * nty * ntz )
    
    !u_rms  = u_rms / ( ntx * nty * ntz )
    !u_rms  = sqrt ( u_rms / 3.0_dp )
    
    !lamda1 = u_rms / sqrt(dudx2)
    !lamda2 = u_rms / sqrt(dvdy2)
    !lamda3 = u_rms / sqrt(dwdz2)

    !adi % Re = Re_taylor / ( u_rms * mean_rho * lamda1 )

    if(rank==rank_default) then 
       
       write(*,*) 'Initial conditions : Homogeneous isotropic turbulence'
       write(*,*)
       write(*,*) 'Mesh Size', ntx, nty, ntz
       write(*,*) 'Total kinetic energy = ' , KE / (ntx * nty * ntz)
       write(*,*) 'Pressure - Density : ', press , rho
       write(*,*) 'Speed of sound : ' , cs
       write(*,*) 'Reynolds : ' , adi % Re
       write(*,*)
       !write(*,*) 'U_rms   : ' , u_rms
       !write(*,*) 'U_mean  : ' , mean_ux , mean_vy , mean_wz
       !write(*,*) 'dudx    : ' , dudx2 , dvdy2 , dwdz2
       !write(*,*) 'Lambda1 : ' , lamda1 , lamda2 , lamda3
       
    end if  

    deallocate ( ux , vy , wz , dux_x , dvy_y , dwz_z , dx_i2 , dummy )
        
  end subroutine ini_HIT


!> \brief "Bogey" shear layer IC.

  subroutine ini_shear_layer_bogey ( inp , adi , grid , T , v )

    type (inp_type)                               , intent (in)    :: inp  !< input derived type
    type (adi_type)                               , intent (in)    :: adi  !< non-dimensional derived type
    type (inp_grid)                               , intent (in)    :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
    integer (ip) :: ok , i , j , k 
    real (dp)    :: u1_ , u2_ , Ma
    real (dp)    :: Deltau , Sigmau , d2O0_i , dw0_i
    real (dp)    :: P , rho , ux , vy , wz


    allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
               T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate ini')


    v   = 0.0_dp
    T   = 0.0_dp

    u1_ = 1.0_dp
    u2_ = 0.5_dp

    Ma = 0.01_dp

    Deltau = u1_ - u2_
    Sigmau = u1_ + u2_

    dw0_i  = 1.6e-3_dp
    d2O0_i = 0.25_dp * dw0_i

    ! invert
    d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
    dw0_i  = 1.0_dp / dw0_i

    ux = 0.0_dp
    vy = 0.0_dp
    wz = 0.0_dp


    do k = sz , ez
       do j = sy , ey
          do i = sx , ex

             T (i,j,k) = ( (0.5_dp / (2.0_dp * Ma)) * (0.5_dp / (2.0_dp * Ma)) ) / adi % gamma
             
             rho = 1.22_dp
             P   = rho * T(i,j,k)

             
             ux = 0.5_dp * ( Sigmau + Deltau * tanh ( grid % y(j) * d2O0_i ) )


             v (i,j,k,1) = rho
             v (i,j,k,2) = rho * ux
             v (i,j,k,3) = rho * vy
             v (i,j,k,4) = rho * wz
             v (i,j,k,5) = P / (adi % gamma - 1.0_dp) + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
             
          end do
       end do
    end do

    
    call upd_prim_var_domain ( adi , T , v )
    call comm_cons (v)
    call upd_boundaries ( inp , adi , 0.0_dp , 0.0_dp , grid , T , v )
    call upd_prim_var_ghost ( adi , T , v )


  end subroutine ini_shear_layer_bogey


  subroutine ini_shear_layer ( inp , adi , grid , T , v )

   type (inp_type)                               , intent (in)    :: inp
   type (adi_type)                               , intent (inout) :: adi
   type (inp_grid)                               , intent (in)    :: grid 
   real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T 
   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v

   real (dp)    :: Press , y_m , d_s , alph
   real (dp)    :: rho , ux , vy , wz 
   real (dp)    :: Mc , rho_1 , rho_2 , c_1 , c_2
   real (dp)    :: eta , geta , dgdy ,k_x , k_z , aij , alph_i , beta_j
   integer (ip) :: i , j , k , ok , n , q , Np


   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
              stat = ok )

   if ( ok > 0 ) call abort_mpi ('error allocate ini')

   rho_0   = 1.0_dp
   delta_0 = 0.2_dp
   
   y_m    = Ly / 2.0_dp
   d_s    = 1.0_dp
   alph   = (d_s - 1.0_dp) / (d_s + 1.0_dp) 

   Mc     = 1.2_dp 
   rho_1  = rho_0 * ( 1.0_dp + alph )
   rho_2  = rho_0 * ( 1.0_dp - alph )

   Press  = 1.0_dp
   
   c_1    = sqrt( adi % gamma * Press / rho_1 )
   c_2    = sqrt( adi % gamma * Press / rho_2 )

   delta_U = Mc * ( c_1 + c_2 )

   Np = 75

   do k = sz , ez 
      do j = sy ,ey 
         do i =  sx , ex

            eta = ( grid % y(j) - y_m ) / delta_0

            rho = rho_0 * (1.0_dp + alph * tanh( eta / 2.0_dp )) 
            ux  = 0.5_dp * delta_U * tanh( eta / 2.0_dp )
            vy  = 0.0_dp
            wz  = 0.0_dp

            geta = exp( -5.0_dp * eta**2 ) * ( sin(eta) + 10.0_dp * eta * cos(eta) )
            dgdy = -10.0_dp * eta * geta + exp( -5.0_dp * eta**2 ) * (11.0_dp * cos(eta) - 10.0_dp * eta * sin(eta))
            dgdy = dgdy / delta_0

            do q = 1 , Np
               do n = 1 , Np 
                  alph_i = 2.0_dp * pi * rand()
                  beta_j = 2.0_dp * pi * rand()
               
                  k_x = 2.0_dp * pi * n / Lx
                  k_z = 2.0_dp * pi * q / Lz
                  aij = 0.15_dp / ( 4.0_dp * pi**2 * k_x * k_z )
                  ux  = ux - aij * sin(k_x * grid % x(i) + alph_i) * cos(k_z * grid % z(k) + beta_j) * k_x * geta
                  vy  = vy + aij * cos(k_x * grid % x(i) + alph_i) * cos(k_z * grid % z(k) + beta_j) * dgdy
                  wz  = wz - aij * cos(k_x * grid % x(i) + alph_i) * sin(k_z * grid % z(k) + beta_j) * k_z * geta
               end do
            end do

            v(i,j,k,1) = rho 
            v(i,j,k,2) = rho * ux
            v(i,j,k,3) = rho * vy 
            v(i,j,k,4) = rho * wz 
            v(i,j,k,5) = Press / ( adi % gamma - 1.0_dp ) + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )

         end do
      end do
   end do

   call upd_prim_var_domain ( adi , T , v )
   call comm_cons (v)
   call upd_boundaries ( inp , adi , 0.0_dp , 0.0_dp , grid , T , v )
   call upd_prim_var_ghost ( adi , T , v )


  end subroutine ini_shear_layer


  subroutine ini_channel_flow ( inp , adi , grid , T , v )

    type (inp_type)                               , intent (in)    :: inp
    type (adi_type)                               , intent (in)    :: adi
    type (inp_grid)                               , intent (in)    :: grid 
    real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T 
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v

    real    (dp) :: rho , ux , vy , wz , pre , u_max , eta , yy
    real    (dp) :: m_trg_tmp , gm13 
    integer (ip) :: i , j , k , ok 

    allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
               T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate ini')

    gm13 = adi%gm1 / 3.0_dp 
    
    rho   = 1.0_dp
    u_max = 1.5_dp
    
    !m_trg_tmp = 0.0_dp

    do i = sx , ex
       do j = sy , ey
          do k = sz , ez

             call random_seed()
             call random_number(eta)
          
             eta = 2.0_dp * ( eta - 0.5_dp )
             yy  = grid%y(j) 

             T(i,j,k) = 1.0_dp! + gm13 * adi % pr * adi % ma**2 * u_max**2 * (1.0_dp - yy**4) 

             ux = u_max * (1.0_dp - yy**2) * ( 1.0_dp + 0.2_dp * eta )
             vy = 0.0_dp
             wz = 0.0_dp
             
             pre = rho * T(i,j,k) / (adi%gamma * adi%ma * adi%ma)
                        
             v(i,j,k,1) = rho 
             v(i,j,k,2) = rho * ux
             v(i,j,k,3) = rho * vy
             v(i,j,k,4) = rho * wz
             v(i,j,k,5) = pre / adi%gm1 + 0.5_dp * rho * ( ux * ux + vy * vy + wz * wz )

            !  m_trg_tmp = m_trg_tmp + v(i,j,k,2)

          end do 
       end do
    end do

   !  call mpi_allreduce ( m_trg_tmp , m_trg , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpicode )
    
   !  m_trg = 1.0_dp * Ly * Lz 

    call upd_prim_var_domain ( adi , T , v )
    call comm_cons (v)
    call upd_boundaries ( inp , adi , 0.0_dp , 0.0_dp , grid , T , v )
    call upd_prim_var_ghost ( adi , T , v )

  end subroutine ini_channel_flow
 

  
  ! subroutine ini_psr ( adi , thd , T , W_i , cp , ha , v )


  !   type (adi_type) , intent (in)                                  :: adi
  !   type (thd_type) , intent (in)                                  :: thd
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v


  !   integer (ip)                                                   :: ok , i , j , k , l
  !   real (dp)                                                      :: P , Temp , rho , Ws_i , hm

  !   real (dp)                                                      :: Ya (nrv+npv+nvv) , has (nrv) , ux , vy , wz

  !   real (dp)                                                      :: t1 , t2

  !   ! real (dp) , allocatable , dimension (:,:,:)                    :: mu , ct
  !   ! real (dp) , allocatable , dimension (:,:,:,:)                  :: Xa , dm


  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini')


  !   ! allocate ( mu   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !   !            ct   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !   !            Xa   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !   !            dm   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !   !            stat = ok )
  !   ! if ( ok > 0 ) call abort_mpi ('error allocate ini 2')


  !   P    = thd % P0  / adi % P_ref
  !   Temp = 1000.0_dp / adi % T_ref

  !   Ya (:)  = 0.0_dp
  !   Ya (1)  = 1.999669132986E-02_dp
  !   Ya (3)  = 9.522160086615E-01_dp
  !   Ya (11) = 2.778730000868E-02_dp

  !   ux = 0.0_dp
  !   vy = 0.0_dp
  !   wz = 0.0_dp

  !   call Wmix_i_scalar ( thd , Ya , Ws_i )
  !   call ha_scalar ( thd , temp , has )

  !   rho = P / ( Temp * Ws_i )

  !   hm = 0.0_dp
  !   do l = 1 , nrv
  !      hm = hm + has (l) * Ya (l)
  !   end do

  !   do k = sz-ng , ez+ng
  !      do j = sy-ng , ey+ng
  !         do i = sx-ng , ex+ng

  !            T (i,j,k)   = Temp ! initial guess for the temperature (mandatory)

  !            v (i,j,k,1) = rho
  !            v (i,j,k,2) = rho * ux
  !            v (i,j,k,3) = rho * vy
  !            v (i,j,k,4) = rho * wz
  !            v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )

  !            do l = 1 , nrv+npv+nvv
  !               v (i,j,k,niv+l) = rho * Ya (l)
  !            end do

  !         end do
  !      end do
  !   end do


  !   call upd_prim_var_domain ( thd , T , W_i , cp , ha , v )
  !   call comm_cons (v)
  !   call upd_boundaries ( adi , thd , T , W_i , cp , ha , v )
  !   call upd_prim_var_ghost ( thd , T , W_i , cp , ha , v )


  !   write (*,*) v (sx,sy,sz,:)

  !   read (*,*)


  ! end subroutine ini_psr


!> \brief "Fedkiw" IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine ini_fedkiw ( adi , thd , x , T , W_i , cp , ha , v )


  !   type (adi_type) , intent (in)                                  :: adi !< non-dimensional derived type
  !   type (thd_type) , intent (in)                                  :: thd !< thermodynamic derived type
  !   real (dp) , allocatable , dimension (:) , intent (in)          :: x   !< x-coordinate array
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T   !< temperature
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i !< inverted molar mass
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp  !< heat capacity
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha  !< especies enthalpy
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v   !< conserved variables array


  !   integer (ip) , dimension (ndimmax+ndimmax) , parameter         :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
  !   integer (ip)                                                   :: ok , i , j , k , l
  !   real (dp)                                                      :: P , Temp , rho , Ws_i , hm

  !   real (dp)                                                      :: Ya (nrv+npv+nvv) , has (nrv) , ux , vy , wz


  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini')


  !   v   = 0.0_dp
  !   ha  = 0.0_dp
  !   T   = 0.0_dp
  !   W_i = 0.0_dp
  !   cp  = 0.0_dp


  !   do k = sz , ez
  !      do j = sy , ey
  !         do i = sx , ex


  !            if ( x (i) * adi % L_ref < 5.0e-2_dp ) then

  !               ! H O H2 O2 OH H2O HO2 H2O2 AR (in chem_fedkiw.inp)
  !               Ya (:)  = 0.0_dp ! Initialization
  !               Ya (3)  = 0.2_dp ! H2
  !               Ya (4)  = 0.1_dp ! O2
  !               Ya (9)  = 0.7_dp ! Ar
  !               call X_to_Y ( thd , Ya )

  !               ux = 0.0_dp
  !               vy = 0.0_dp
  !               wz = 0.0_dp

  !               P    = 8000.0_dp / adi % P_ref
  !               Temp = 400.0_dp  / adi % T_ref

  !               call Wmix_i_scalar ( thd , Ya , Ws_i )
  !               call ha_scalar ( thd , temp , has )

  !               hm = 0.0_dp
  !               do l = 1 , nrv
  !                  hm = hm + has (l) * Ya (l)
  !               end do

  !               rho = P / ( Temp * Ws_i )


  !               T (i,j,k)   = Temp ! initial guess for the temperature (mandatory)

  !               v (i,j,k,1) = rho
  !               v (i,j,k,2) = rho * ux
  !               v (i,j,k,3) = rho * vy
  !               v (i,j,k,4) = rho * wz
  !               v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )

  !               do l = 1 , nrv+npv+nvv
  !                  v (i,j,k,niv+l) = rho * Ya (l)
  !               end do

  !            else

  !               Ya (:)  = 0.0_dp ! Initialization
  !               Ya (3)  = 0.2_dp ! H2
  !               Ya (4)  = 0.1_dp ! O2
  !               Ya (9)  = 0.7_dp ! Ar
  !               call X_to_Y ( thd , Ya )

  !               ux = 0.0_dp
  !               vy = 0.0_dp
  !               wz = 0.0_dp

  !               P    = 80000.0_dp / adi % P_ref
  !               Temp = 1200.0_dp  / adi % T_ref

  !               call Wmix_i_scalar ( thd , Ya , Ws_i )
  !               call ha_scalar ( thd , temp , has )

  !               hm = 0.0_dp
  !               do l = 1 , nrv
  !                  hm = hm + has (l) * Ya (l)
  !               end do

  !               rho = P / ( Temp * Ws_i )


  !               T (i,j,k)   = Temp ! initial guess for the temperature (mandatory)

  !               v (i,j,k,1) = rho
  !               v (i,j,k,2) = rho * ux
  !               v (i,j,k,3) = rho * vy
  !               v (i,j,k,4) = rho * wz
  !               v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )

  !               do l = 1 , nrv+npv+nvv
  !                  v (i,j,k,niv+l) = rho * Ya (l)
  !               end do

  !            end if


  !         end do
  !      end do
  !   end do


  !   call comm_cons (v)
  !   call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
  !   call upd_extrapolation (v)
  !   do l = 1 , ndim+ndim
  !      call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
  !   end do


  ! end subroutine ini_fedkiw


!> \brief Double Mach reflection IC.
!!
!! Initialization of the 2D test case of the double Mach reflection
!! (DMR) on a wall with an incident _unsteady_ shock wave arriving
!! bias on a horizontal plane wall.
!!
!!                          shock
!!                            ||
!!                          --||--> Us
!!           (2)              ||                 (1)
!!                            ||        u1=000, T1=300K, P1=1atm
!!
!! References :
!!
!! -# Vasilev et al. "The wall-jetting effect in Mach reflection:
!!    Navier–Stokes simulations," Journal of Fluid Mechanics.
!! -# http://ufrmeca.univ-lyon1.fr/~buffat/COURS/AERO_HTML/node49.html
!! -# A.Kourta, 'Dynamique des gaz, Supersoniques linearises et ondes
!!    instationnaires, Fascicule de cours', Polytech Orléans, 2012.
!! -# R.E.Mitchell, R.J.Kee, 'A general-purpose computer code for
!!    predicting chemical kinetic behavior behind incident and
!!    reflected shocks', Sandia Report, 1982, p.13.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine ini_DMR ( adi , thd , x , y , T , W_i , cp , ha , v )


  !   type (adi_type) , intent (in)                                  :: adi
  !   type (thd_type) , intent (in)                                  :: thd
  !   real (dp) , allocatable , dimension (:) , intent (in)          :: x
  !   real (dp) , allocatable , dimension (:) , intent (in)          :: y
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v


  !   integer (ip) , dimension (ndimmax+ndimmax) , parameter         :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
  !   integer (ip)                                                   :: ok , i , j , k , l
  !   real (dp)                                                      :: Ws_i ! mixing molar mass (Left & Right)

  !   real (dp)  :: Ya (nrv+npv+nvv) ! mixing composition

  !   real (dp)  :: Ms, Us ! shock Mach number & shock velocity
  !   real (dp)  :: u1, v1, w1, vnorme1, p1, T1, rho1, cs1, hm1, &   ! pre-shock variables (RIGHT)
  !                 cp1, gamma1
  !   real (dp)  :: u2, v2, w2, vnorme2, p2, T2, rho2, hm2      ! post-shock variables (LEFT)

  !   real (dp)             :: x0      ,& ! abscisse of the shock impact
  !                            thetaw  ,& ! angle of the angled wedges
  !                            beta    ,& ! angle=(wall, shock)
  !                            slope   ,& ! = tan(beta)
  !                            slope_i ,& ! = 1/tan(beta)
  !                            theta      ! angle(wall, velocity magnitude after shock)

  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini')

  !   v   = 0.0_dp
  !   ha  = 0.0_dp
  !   T   = 0.0_dp
  !   W_i = 0.0_dp
  !   cp  = 0.0_dp

  !   !====== Parameters of the DMR problem ======
  !       ! Species : O2   N2 (in chem.inp)
  !       Ya (1)  = 0.00_dp ! 02 (0.23_dp for air)
  !       Ya (2)  = 1.00_dp ! N2 (0.77_dp for air)

  !       ! Angled wedge
  !       x0      = 0.0_dp / adi % L_ref
  !       thetaw  = 36.0_dp * pi/180.0_dp
  !       beta    = pi/2.0_dp - thetaw
  !       slope   = tan(beta)
  !       slope_i = 1.0_dp/slope

  !       Ms = 4.5_dp  ! shock Mach number

  !       ! Pre-shock conditions
  !       u1 = 0.0_dp
  !       v1 = 0.0_dp
  !       w1 = 0.0_dp
  !       vnorme1 = sqrt(u1*u1 + v1*v1 + w1*w1)
  !       P1 = 101325.0_dp / adi % P_ref
  !       T1 = 300.0_dp    / adi % T_ref
  !   !===========================================

  !   call Wmix_i_scalar ( thd , Ya , Ws_i )

  !   ! Calculate of pre-shock conditions (RIGHT)
  !       rho1 = P1 / ( T1 * Ws_i )

  !       call cp_scalar ( thd , T1 , Ya , cp1 )
  !       gamma1 = thd % gam2 * Ws_i / cp1
  !       gamma1 = 1.0_dp / ( 1.0_dp - gamma1 )
  !       cs1    = sqrt ( gamma1 * p1 / rho1 )

  !   Us = Ms * cs1 !'cs1' because of wave speed is measured from the speed of sound upstream of shock

  !   ! Calculate of post-shock conditions (LEFT)
  !       ! Using RH subroutine, by assuming the shock is STEADY, in shock point of view
  !           vnorme1 = Us ! in shock point of view (STEADY shock)
  !           call RH_variables (thd, Ya, p1, T1, rho1, vnorme1, hm1, pi/2.0_dp, &
  !                              p2, T2, rho2, vnorme2, hm2, theta)
  !           ! where vnorme2 is now the velocity relative to the reference related to the STEADY shock
  !           vnorme2 = Us - vnorme2 ! conversion to the reference related to the laboratory (UNSTEADY shock)
  !           vnorme1 = 0.0_dp ! put back the rest state upstream of the shock (UNSTEADY shock)
  !       ! END of using RH subroutine

  !       u2 =  vnorme2 * sin(beta) ! velocity relative to cartesian coord.
  !       v2 = -vnorme2 * cos(beta) ! v2 < 0 because flow direction downwards
  !       w2 = 0.0_dp

  !   ! End of calculation of pre- and post-shock conditions

  !   ! Application of (T, p, rho, u, v, w, hm) in the calculation domain
  !   do k = sz , ez ! z-direction
  !      do j = sy , ey ! y-direction
  !         do i = sx , ex ! x-direction

  !            if ( x(i) <= (x0+y(j)*slope_i) .and. y(j) >= (slope*(x(i)-x0)) ) then
  !               ! POST-shock conditions (left)

  !               T (i,j,k)   = T2 ! initial guess for the temperature (mandatory)

  !               v (i,j,k,1) = rho2
  !               v (i,j,k,2) = rho2 * u2
  !               v (i,j,k,3) = rho2 * v2
  !               v (i,j,k,4) = rho2 * w2
  !               v (i,j,k,5) = rho2 * hm2 - P2 + 0.5_dp * rho2 * (u2*u2 + v2*v2 + w2*w2)

  !               do l = 1 , nrv+npv+nvv
  !                  v (i,j,k,niv+l) = rho2 * Ya (l)
  !               end do

  !            else
  !               ! PRE-shock conditions (right)

  !               T (i,j,k)   = T1 ! initial guess for the temperature (mandatory)

  !               v (i,j,k,1) = rho1
  !               v (i,j,k,2) = rho1 * u1
  !               v (i,j,k,3) = rho1 * v1
  !               v (i,j,k,4) = rho1 * w1
  !               v (i,j,k,5) = rho1 * hm1 - P1 + 0.5_dp * rho1 * (u1*u1 + v1*v1 + w1*w1)

  !               do l = 1 , nrv+npv+nvv
  !                  v (i,j,k,niv+l) = rho1 * Ya (l)
  !               end do

  !            end if

  !         end do
  !      end do
  !   end do


  !   call comm_cons (v)
  !   call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
  !   call upd_extrapolation (v)
  !   do l = 1 , ndim+ndim
  !      call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
  !   end do


  ! end subroutine ini_DMR



!> \brief Poiseuille IC.
!!
!! Poiseuille flow with white noise at ambient conditions (T=300K, P=1atm).
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
!   subroutine ini_poiseuille ( adi , thd , x , y , z , T , W_i , cp , ha , v )


!     type (adi_type) , intent (in)                                  :: adi
!     type (thd_type) , intent (in)                                  :: thd
!     real (dp) , allocatable , dimension (:) , intent (in)          :: x
!     real (dp) , allocatable , dimension (:) , intent (in)          :: y
!     real (dp) , allocatable , dimension (:) , intent (in)          :: z
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v

!     integer (ip) , dimension (ndimmax+ndimmax) , parameter         :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
!     integer (ip)                                                   :: ok , i , j , k , l
!     real (dp)                                                      :: Ws_i

!     real (dp)               :: Ya (nrv+npv+nvv) , has (nrv)
!     real (dp)               :: ux, vy, wz, Pamb, Tamb, rho, hm, Ma, cp0, gamma0
!     real (dp)               :: ux_inf , wrky , wrkz , eps , noise
!     real (dp)               :: alpha , beta


!     allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
!                ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
!                T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                stat = ok )
!     if ( ok > 0 ) then
!        write (*,*) 'error allocate ini_amb'
!        call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
!     end if


!     v   = 0.0_dp
!     ha  = 0.0_dp
!     T   = 0.0_dp
!     W_i = 0.0_dp
!     cp  = 0.0_dp

!     !======== Parameters of the problem =========
!         ! Species : in chem.inp
!         Ya(:) = 0.0_dp
!         if      (thd % Nspc == 1) then
!             Ya (1) = 1.000_dp ! N2
!         else if (thd % Nspc == 2) then ! air only
!             Ya (1) = 0.233_dp ! O2 (0.233_dp for air)
!             Ya (2) = 0.767_dp ! N2 (0.767_dp for air)
!         else
!             write (*,*) 'module_ICs : Problem with species initialisation'
!             call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
!         end if

!         ux = 0.0_dp
!         vy = 0.0_dp
!         wz = 0.0_dp

!         Pamb = 101325.0_dp / adi % P_ref
!         Tamb = 300.0_dp    / adi % T_ref
!         Ma   = 0.25_dp

!         ! Noise amplitude in percent
!         noise = 0.80_dp ! in %

!         alpha = 0.04_dp
!         beta = 1.0_dp
!     !===========================================

!     noise = noise / 100_dp

!     call Wmix_i_scalar ( thd , Ya , Ws_i )
!     rho = Pamb / ( Tamb * Ws_i )

!     call cp_scalar ( thd , Tamb , Ya , cp0 )
!     gamma0 = thd % gam2 * Ws_i / cp0
!     gamma0 = 1.0_dp / ( 1.0_dp - gamma0 )
!     ux = Ma * sqrt ( gamma0 * Pamb / rho )

!     ux_inf = ux

!     call ha_scalar ( thd , Tamb , has )
!     hm = 0.0_dp
!     do l = 1 , nrv
!        hm = hm + has (l) * Ya (l)
!     end do

!     ! Application of (T, p, rho, u, v, w, hm) in the calculation domain
!     do k = sz , ez ! z-direction
!        do j = sy , ey ! y-direction
!           do i = sx , ex ! x-direction

!                 T (i,j,k)   = Tamb ! initial guess for the temperature (mandatory)

! !                wrky = ( y(j) - Ly ) / Ly ! (-y,+y) = (adiabaticwall,symmetry)
! !                wrkz = ( z(k) - Lz ) / Lz ! (-z,+z) = extrapolation/periodic
! !                ux = ux_inf * ( 1.0_dp - wrky * wrky )

!                 wrky = 2.0_dp * y(j) / Ly ! (-y,+y) = (adiabaticwall,adiabaticwall)
!                 wrkz = 2.0_dp * z(k) / Lz ! (-z,+z) = (adiabaticwall,adiabaticwall)
!                 ux = ux_inf * sqrt ( ( 1.0_dp - wrky * wrky ) * ( 1.0_dp - wrkz * wrkz ) )

! !                ux = ux_inf

! !                wrky = erf ( ( abs(y(j)) - beta * 0.5_dp * Ly ) / ( alpha * Ly ) )
! !                wrkz = erf ( ( abs(z(k)) - beta * 0.5_dp * Lz ) / ( alpha * Lz ) )
! !                ux = 0.5_dp * ux_inf * sqrt ( ( 1.0_dp - wrky ) * ( 1.0_dp - wrkz ) )

!                 ! add white noise
!                 call random_number (eps)
!                 eps = noise * ( eps + eps - 1.0_dp )
!                 ux = ux * ( 1.0_dp + eps ) 

!                 v (i,j,k,1) = rho
!                 v (i,j,k,2) = rho * ux
!                 v (i,j,k,3) = rho * vy
!                 v (i,j,k,4) = rho * wz
!                 v (i,j,k,5) = rho * hm - Pamb + 0.5_dp * rho * (ux*ux + vy*vy + wz*wz)

!                 do l = 1 , nrv+npv+nvv
!                    v (i,j,k,niv+l) = rho * Ya (l)
!                 end do

!           end do
!        end do
!     end do


!     call comm_cons (v)
!     call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
!     call upd_extrapolation (v)
!     do l = 1 , ndim+ndim
!        call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
!     end do


!   end subroutine ini_poiseuille


!   subroutine ini_vicquelin_periodic ( adi , thd , x , y , T , W_i , cp , ha , v )


!     type (adi_type) , intent (in)                                  :: adi
!     type (thd_type) , intent (in)                                  :: thd
!     real (dp) , allocatable , dimension (:) , intent (in)          :: x
!     real (dp) , allocatable , dimension (:) , intent (in)          :: y
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v


!     integer (ip) , dimension (ndimmax+ndimmax) , parameter         :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)

!     real (dp)                                                      :: Tmax , Tmin , &
!                                                                       x0 , d0

!     integer (ip)                                                   :: ok , i , j , k , l

!     real (dp)                                                      :: P , rho , Ws_i , hm , zeta

!     real (dp)                                                      :: Ya (nrv+npv+nvv) , has (nrv) , ux , vy , wz


!     allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
!                ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
!                T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                stat = ok )
!     if ( ok > 0 ) call abort_mpi ('error allocate ini')


!     v   = 0.0_dp
!     ha  = 0.0_dp
!     T   = 0.0_dp
!     W_i = 0.0_dp
!     cp  = 0.0_dp


!     Tmax = 320.0_dp / adi % T_ref
!     Tmin = 320.0_dp / adi % T_ref
!     x0   = 2.5e-3_dp / adi % L_ref
!     d0   = 2.5e-4_dp / adi % L_ref
!     P    = thd % P0 / adi % P_ref
!     ux   = 10.0_dp / adi % u_ref
!     vy   = 0.0_dp / adi % u_ref
!     wz   = 0.0_dp / adi % u_ref




!     do k = sz , ez
!        do j = sy , ey
!           do i = sx , ex


!              zeta =  1.0_dp - exp ( - ( x(i) - x0 ) * ( x(i) - x0 ) / ( d0 * d0 ) )

!              T (i,j,k) = ( Tmin - Tmax ) * zeta + Tmax

!              Ya (1) = ( 0.195_dp - 0.142_dp ) * zeta + 0.142_dp
!              Ya (2) = ( 0.591_dp - 0.758_dp ) * zeta + 0.758_dp
!              Ya (3) = ( 0.0_dp   - 0.1_dp )   * zeta + 0.1_dp
!              Ya (4) = ( 0.214_dp - 0.0_dp )   * zeta + 0.0_dp

!              call Wmix_i_scalar ( thd , Ya , Ws_i )
!              rho = P / ( T (i,j,k) * Ws_i )

!              call ha_scalar ( thd , T (i,j,k) , has )
!              hm = 0.0_dp
!              do l = 1 , nrv
!                 hm = hm + has (l) * Ya (l)
!              end do

!              v (i,j,k,1) = rho
!              v (i,j,k,2) = rho * ux
!              v (i,j,k,3) = rho * vy
!              v (i,j,k,4) = rho * wz
!              v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
!              do l = 1 , nrv+npv+nvv
!                 v (i,j,k,niv+l) = rho * Ya (l)
!              end do


!           end do
!        end do
!     end do


!     call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
!     call comm_cons (v)
!     call upd_boundaries ( adi , thd , y , T , W_i , cp , ha , v )

!     do l = 1 , ndim+ndim
!        call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
!     end do


!   end subroutine ini_vicquelin_periodic


!   subroutine ini_vicquelin_temp_val ( adi , thd , time , x , y , z , T , W_i , cp , ha , v )


!     type (adi_type) , intent (in)                                  :: adi
!     type (thd_type) , intent (in)                                  :: thd
!     real (dp) , intent (in)                                        :: time
!     real (dp) , allocatable , dimension (:) , intent (in)          :: x
!     real (dp) , allocatable , dimension (:) , intent (in)          :: y
!     real (dp) , allocatable , dimension (:) , intent (in)          :: z
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v


!     integer (ip) , dimension (ndimmax+ndimmax) , parameter         :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)

!     real (dp)                                                      :: Tmax , Tmin , x0 , d0 , Lx_i

!     integer (ip)                                                   :: ok , i , j , k , l

!     real (dp)                                                      :: P , rho , Ws_i , hm , zeta

!     real (dp)                                                      :: Ya (nrv+npv+nvv) , has (nrv) , ux , vy , wz


!     allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
!                ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
!                T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                stat = ok )
!     if ( ok > 0 ) call abort_mpi ('error allocate ini')


!     v   = 0.0_dp
!     ha  = 0.0_dp
!     T   = 0.0_dp
!     W_i = 0.0_dp
!     cp  = 0.0_dp


!     Tmax = 1350.0_dp / adi % T_ref
!     Tmin = 320.0_dp / adi % T_ref
!     x0   = 2.5e-3_dp / adi % L_ref
!     d0   = 5.0e-4_dp / adi % L_ref
!     P    = thd % P0 / adi % P_ref
! !    ux   = 10.0_dp / adi % u_ref
!     vy   = 0.0_dp / adi % u_ref
!     wz   = 0.0_dp / adi % u_ref
!     Lx_i = adi % L_ref / 5.0e-2_dp


!     do k = sz , ez
!        do j = sy , ey
!           do i = sx , ex


!              zeta = 0.0_dp !1.0_dp - x(i) * Lx_i !- exp ( - ( x(i) - x0 ) * ( x(i) - x0 ) / ( d0 * d0 ) )

!              ux   = 0.01_dp * 734.6_dp * exp ( - ( x(i) - x0 ) * ( x(i) - x0 ) / ( d0 * d0 ) ) / adi % u_ref

!              T (i,j,k) = ( Tmin - Tmax ) * zeta + Tmax

!              ! Ya (1) = ( 0.195_dp - 0.142_dp ) * zeta + 0.142_dp
!              ! Ya (2) = ( 0.591_dp - 0.758_dp ) * zeta + 0.758_dp
!              ! Ya (3) = ( 0.0_dp   - 0.1_dp )   * zeta + 0.1_dp
!              ! Ya (4) = ( 0.214_dp - 0.0_dp )   * zeta + 0.0_dp

!              Ya (1) = ( 0.195_dp - 0.142_dp ) * zeta + 0.142_dp
!              Ya (2) = ( 0.214_dp - 0.0_dp )   * zeta + 0.0_dp
!              Ya (3) = ( 0.0_dp   - 0.1_dp )   * zeta + 0.1_dp
!              Ya (4) = ( 0.591_dp - 0.758_dp ) * zeta + 0.758_dp

!              call Wmix_i_scalar ( thd , Ya , Ws_i )

!              P   = thd % P0 / adi % P_ref + 734.6_dp * 0.24_dp * ux / ( adi % rho_ref * adi % u_ref )

!              rho = ( 0.24_dp / adi % rho_ref ) * ( 1.0_dp + ux * adi % u_ref / 734.6_dp )

!              T (i,j,k) = P / ( rho * Ws_i )

!              call ha_scalar ( thd , T (i,j,k) , has )
!              hm = 0.0_dp
!              do l = 1 , nrv
!                 hm = hm + has (l) * Ya (l)
!              end do

!              v (i,j,k,1) = rho
!              v (i,j,k,2) = rho * ux
!              v (i,j,k,3) = rho * vy
!              v (i,j,k,4) = rho * wz
!              v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
!              do l = 1 , nrv+npv+nvv
!                 v (i,j,k,niv+l) = rho * Ya (l)
!              end do


!           end do
!        end do
!     end do


!     call comm_cons (v)
!     call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
!     call upd_extrapolation (v)
!     do l = 1 , ndim+ndim
!        call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
!     end do


!   end subroutine ini_vicquelin_temp_val3


!> \brief Diffusion IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine ini_diffusion ( adi , thd , x , T , W_i , cp , ha , v )


  !   type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
  !   type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
  !   real (dp) , allocatable , dimension (:) , intent (in)          :: x    !< x-coordinate array
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


  !   integer (ip) , dimension (ndimmax+ndimmax) , parameter         :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
  !   real (dp)                                                      :: Tmax , Tmin , &
  !                                                                     x0 , d0
  !   integer (ip)                                                   :: ok , i , j , k , l
  !   real (dp)                                                      :: P , rho , Ws_i , hm , zeta
  !   real (dp)                                                      :: Ya (nrv+npv+nvv) , has (nrv) , ux , vy , wz


  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini')


  !   v   = 0.0_dp
  !   ha  = 0.0_dp
  !   T   = 0.0_dp
  !   W_i = 0.0_dp
  !   cp  = 0.0_dp


  !   Tmax = 1350.0_dp / adi % T_ref
  !   Tmin = 320.0_dp / adi % T_ref
  !   x0   = 2.5e-2_dp / adi % L_ref
  !   d0   = 2.5e-3_dp / adi % L_ref
  !   P    = thd % P0 / adi % P_ref
  !   ux   = 0.0_dp / adi % u_ref
  !   vy   = 0.0_dp / adi % u_ref
  !   wz   = 0.0_dp / adi % u_ref


  !   do k = sz , ez
  !      do j = sy , ey
  !         do i = sx , ex


  !            zeta =  1.0_dp - 0.5_dp * exp ( - ( x(i) - x0 ) * ( x(i) - x0 ) / ( d0 * d0 ) )

  !            T (i,j,k) = ( Tmin - Tmax ) * zeta + Tmax

  !            Ya (1) = ( 0.195_dp - 0.142_dp ) * zeta + 0.142_dp
  !            Ya (2) = ( 0.214_dp - 0.0_dp )   * zeta + 0.0_dp
  !            Ya (3) = ( 0.0_dp   - 0.1_dp )   * zeta + 0.1_dp
  !            Ya (4) = ( 0.591_dp - 0.758_dp ) * zeta + 0.758_dp

  !            call Wmix_i_scalar ( thd , Ya , Ws_i )
  !            rho = P / ( T (i,j,k) * Ws_i )

  !            call ha_scalar ( thd , T (i,j,k) , has )
  !            hm = 0.0_dp
  !            do l = 1 , nrv
  !               hm = hm + has (l) * Ya (l)
  !            end do

  !            v (i,j,k,1) = rho
  !            v (i,j,k,2) = rho * ux
  !            v (i,j,k,3) = rho * vy
  !            v (i,j,k,4) = rho * wz
  !            v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
  !            do l = 1 , nrv+npv+nvv
  !               v (i,j,k,niv+l) = rho * Ya (l)
  !            end do


  !         end do
  !      end do
  !   end do


  !   call comm_cons (v)
  !   call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
  !   call upd_extrapolation (v)
  !   do l = 1 , ndim+ndim
  !      call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
  !   end do


  ! end subroutine ini_diffusion


!> \brief Air shear layer IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine ini_shear_layer ( adi , thd , y , T , W_i , cp , ha , v )


  !   type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
  !   type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
  !   real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


  !   integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
  !   integer (ip) :: ok , i , j , k , l
  !   real (dp)    :: u1_ , u2_ , cs1 , cs2 , uc , ym2
  !   real (dp)    :: Deltau , Sigmau , d2O0_i , dw0_i
  !   real (dp)    :: P , rho , Ws_i , hm , ux , vy , wz
  !   real (dp)    :: Ya (nrv+npv+nvv) , has (nrv)


  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini')


  !   v   = 0.0_dp
  !   ha  = 0.0_dp
  !   T   = 0.0_dp
  !   W_i = 0.0_dp
  !   cp  = 0.0_dp


  !   u1_  = 1780.0_dp / adi % u_ref ; cs1 = 1771.89_dp / adi % u_ref
  !   u2_  = 1420.0_dp / adi % u_ref ; cs2 = 714.59_dp / adi % u_ref

  !   Deltau = u1_ - u2_
  !   Sigmau = u1_ + u2_

  !   d2O0_i = 8.88e-5_dp / adi % L_ref
  !   dw0_i  = 4.0_dp * d2O0_i

  !   uc = ( cs1 * u1_ + cs2 * u2_ ) / ( cs1 + cs2 )

  !   ! invert
  !   d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
  !   dw0_i  = 1.0_dp / dw0_i

  !   ym2 = 0.0_dp

  !   vy = 0.0_dp
  !   wz = 0.0_dp


  !   do k = sz , ez
  !      do j = sy , ey
  !         do i = sx , ex


  !            if ( y (j) >= ym2 ) then


  !               T (i,j,k) = 545.0_dp / adi % T_ref
  !               P         = 112.0e3_dp / adi % P_ref

  !               Ya (:)  = 0.0_dp
  !               Ya (1)  = 1.0_dp

  !               call Wmix_i_scalar ( thd , Ya , Ws_i )
  !               rho = P / ( T (i,j,k) * Ws_i )

  !               call ha_scalar ( thd , T (i,j,k) , has )
  !               hm = 0.0_dp
  !               do l = 1 , nrv
  !                  hm = hm + has (l) * Ya (l)
  !               end do

  !               ux = 0.5_dp * ( Sigmau + Deltau * tanh ( ( y(j) - ym2 ) * d2O0_i ) )


  !               v (i,j,k,1) = rho
  !               v (i,j,k,2) = rho * ux
  !               v (i,j,k,3) = rho * vy
  !               v (i,j,k,4) = rho * wz
  !               v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
  !               do l = 1 , nrv+npv+nvv
  !                  v (i,j,k,niv+l) = rho * Ya (l)
  !               end do


  !            else


  !               T (i,j,k) = 1250.0_dp / adi % T_ref
  !               P         = 107.0e3_dp / adi % P_ref

  !               Ya (:)  = 0.0_dp
  !               Ya (2)  = 0.201_dp
  !               Ya (3)  = 0.255_dp
  !               Ya (4)  = 0.544_dp
  !               call X_to_Y ( thd , Ya )

  !               call Wmix_i_scalar ( thd , Ya , Ws_i )
  !               rho = P / ( T (i,j,k) * Ws_i )

  !               call ha_scalar ( thd , T (i,j,k) , has )
  !               hm = 0.0_dp
  !               do l = 1 , nrv
  !                  hm = hm + has (l) * Ya (l)
  !               end do

  !               ux = 0.5_dp * ( Sigmau + Deltau * tanh ( ( y(j) - ym2 ) * d2O0_i ) )


  !               v (i,j,k,1) = rho
  !               v (i,j,k,2) = rho * ux
  !               v (i,j,k,3) = rho * vy
  !               v (i,j,k,4) = rho * wz
  !               v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
  !               do l = 1 , nrv+npv+nvv
  !                  v (i,j,k,niv+l) = rho * Ya (l)
  !               end do


  !            end if


  !         end do
  !      end do
  !   end do


  !   call comm_cons (v)
  !   call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
  !   call upd_extrapolation (v)
  !   do l = 1 , ndim+ndim
  !      call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
  !   end do


  ! end subroutine ini_shear_layer







!> \brief "Fu" shear layer IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine ini_shear_layer_fu ( adi , thd , y , T , W_i , cp , ha , v )


  !   type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
  !   type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
  !   real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


  !   integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
  !   integer (ip) :: ok , i , j , k , l
  !   real (dp)    :: u1_ , u2_
  !   real (dp)    :: Deltau , Sigmau , d2O0_i , dw0_i
  !   real (dp)    :: P , rho , Ws_i , hm , ux , vy , wz
  !   real (dp)    :: Ya (nrv+npv+nvv) , has (nrv)


  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini')


  !   v   = 0.0_dp
  !   ha  = 0.0_dp
  !   T   = 0.0_dp
  !   W_i = 0.0_dp
  !   cp  = 0.0_dp


  !   u1_ = 647.653_dp / adi % u_ref
  !   u2_ = 374.957_dp / adi % u_ref

  !   Deltau = u1_ - u2_
  !   Sigmau = u1_ + u2_

  !   dw0_i  = 3.53e-5_dp / adi % L_ref
  !   d2O0_i = 0.25_dp * dw0_i

  !   ! invert
  !   d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
  !   dw0_i  = 1.0_dp / dw0_i

  !   ux = 0.0_dp
  !   vy = 0.0_dp
  !   wz = 0.0_dp


  !   do k = sz , ez
  !      do j = sy , ey
  !         do i = sx , ex


  !            T (i,j,k) = 288.0_dp / adi % T_ref
  !            P         = 1.0e5_dp / adi % P_ref

  !            Ya (1) = 0.21_dp
  !            Ya (2) = 0.79_dp
  !            call X_to_Y ( thd , Ya )
  !            call Wmix_i_scalar ( thd , Ya , Ws_i )
  !            rho = P / ( T (i,j,k) * Ws_i )

  !            call ha_scalar ( thd , T (i,j,k) , has )
  !            hm = 0.0_dp
  !            do l = 1 , nrv
  !               hm = hm + has (l) * Ya (l)
  !            end do

  !            ux = 0.5_dp * ( Sigmau + Deltau * tanh ( y(j) * d2O0_i ) )


  !            v (i,j,k,1) = rho
  !            v (i,j,k,2) = rho * ux
  !            v (i,j,k,3) = rho * vy
  !            v (i,j,k,4) = rho * wz
  !            v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
  !            do l = 1 , nrv+npv+nvv
  !               v (i,j,k,niv+l) = rho * Ya (l)
  !            end do


  !         end do
  !      end do
  !   end do


  !   call comm_cons (v)
  !   call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
  !   call upd_extrapolation (v)
  !   do l = 1 , ndim+ndim
  !      call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
  !   end do


  ! end subroutine ini_shear_layer_fu


!> \brief "Mahle" shear layer IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
!   subroutine ini_shear_layer_mahle ( adi , thd , y , T , W_i , cp , ha , v )


!     type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
!     type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
!     real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< x-coordinate array
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


!     integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
!     integer (ip) :: ok , i , j , k , l
!     real (dp)    :: u1_ , u2_
!     real (dp)    :: Deltau , Sigmau , d2O0_i , dw0_i
!     real (dp)    :: P , rho , Ws_i , hm , ux , vy , wz
!     real (dp)    :: Ya (nrv+npv+nvv) , has (nrv)


!     allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
!                ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
!                T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                stat = ok )
!     if ( ok > 0 ) call abort_mpi ('error allocate ini')


!     v   = 0.0_dp
!     ha  = 0.0_dp
!     T   = 0.0_dp
!     W_i = 0.0_dp
!     cp  = 0.0_dp


!     u1_ = 648.4_dp / adi % u_ref
!     u2_ = 380.7_dp / adi % u_ref

!     Deltau = u1_ - u2_
!     Sigmau = u1_ + u2_

!     dw0_i  = 3.59e-5_dp / adi % L_ref
!     d2O0_i = 0.25_dp * dw0_i

!     ! invert
!     d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
!     dw0_i  = 1.0_dp / dw0_i

!     ux = 0.0_dp
!     vy = 0.0_dp
!     wz = 0.0_dp


!     do k = sz , ez
!        do j = sy , ey
!           do i = sx , ex


! !             if ( y (j) >= 0.0_dp ) then


!                 T (i,j,k) = 288.0_dp / adi % T_ref
!                 P         = 1.0e5_dp / adi % P_ref

!                 ! Ya (:) = 0.0_dp
!                 ! Ya (1) = 1.0_dp

!                 Ya (1) = 0.5_dp * ( 1.0_dp + 1.0_dp * tanh ( y(j) * d2O0_i ) )
!                 Ya (2) = 0.5_dp * ( 1.0_dp - 1.0_dp * tanh ( y(j) * d2O0_i ) )

!                 call Wmix_i_scalar ( thd , Ya , Ws_i )
!                 rho = P / ( T (i,j,k) * Ws_i )

!                 call ha_scalar ( thd , T (i,j,k) , has )
!                 hm = 0.0_dp
!                 do l = 1 , nrv
!                    hm = hm + has (l) * Ya (l)
!                 end do

!                 ux = 0.5_dp * ( Sigmau + Deltau * tanh ( y(j) * d2O0_i ) )


!                 v (i,j,k,1) = rho
!                 v (i,j,k,2) = rho * ux
!                 v (i,j,k,3) = rho * vy
!                 v (i,j,k,4) = rho * wz
!                 v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
!                 do l = 1 , nrv+npv+nvv
!                    v (i,j,k,niv+l) = rho * Ya (l)
!                 end do


!              ! else


!              !    T (i,j,k) = 288.0_dp / adi % T_ref
!              !    P         = 1.0e5_dp / adi % P_ref

!              !    Ya (:)  = 0.0_dp
!              !    Ya (2)  = 1.0_dp
!              !    call X_to_Y ( thd , Ya )

!              !    call Wmix_i_scalar ( thd , Ya , Ws_i )
!              !    rho = P / ( T (i,j,k) * Ws_i )

!              !    call ha_scalar ( thd , T (i,j,k) , has )
!              !    hm = 0.0_dp
!              !    do l = 1 , nrv
!              !       hm = hm + has (l) * Ya (l)
!              !    end do

!              !    ux = 0.5_dp * ( Sigmau + Deltau * tanh ( y(j) * d2O0_i ) )


!              !    v (i,j,k,1) = rho
!              !    v (i,j,k,2) = rho * ux
!              !    v (i,j,k,3) = rho * vy
!              !    v (i,j,k,4) = rho * wz
!              !    v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
!              !    do l = 1 , nrv+npv+nvv
!              !       v (i,j,k,niv+l) = rho * Ya (l)
!              !    end do


!              ! end if


!           end do
!        end do
!     end do


!     call comm_cons (v)
!     call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
!     call upd_extrapolation (v)
!     do l = 1 , ndim+ndim
!        call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
!     end do


!   end subroutine ini_shear_layer_mahle


  ! subroutine ini_shear_layer_cheng ( adi , thd , time , x , y , z , T , W_i , cp , ha , v )


  !   type (adi_type) , intent (in)                                  :: adi
  !   type (thd_type) , intent (in)                                  :: thd
  !   real (dp) , intent (in)                                        :: time
  !   real (dp) , allocatable , dimension (:) , intent (in)          :: x
  !   real (dp) , allocatable , dimension (:) , intent (in)          :: y
  !   real (dp) , allocatable , dimension (:) , intent (in)          :: z
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v


  !   integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)

  !   integer (ip) :: ok , i , j , k , l
  !   real (dp)    :: u1_ , u2_
  !   real (dp)    :: Deltau , Sigmau , d2O0_i , dw0_i
  !   real (dp)    :: P , rho , Ws_i , hm , ux , vy , wz
  !   real (dp)    :: Ya (nrv+npv+nvv) , has (nrv)


  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini')


  !   v   = 0.0_dp
  !   ha  = 0.0_dp
  !   T   = 0.0_dp
  !   W_i = 0.0_dp
  !   cp  = 0.0_dp


  !   u1_ = 1908.0_dp / adi % u_ref
  !   u2_ = 928.41_dp / adi % u_ref

  !   Deltau = u1_ - u2_
  !   Sigmau = u1_ + u2_

  !   dw0_i  = 1.283e-4_dp / adi % L_ref
  !   d2O0_i = 0.25_dp * dw0_i

  !   ! invert
  !   d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
  !   dw0_i  = 1.0_dp / dw0_i

  !   ux = 0.0_dp
  !   vy = 0.0_dp
  !   wz = 0.0_dp


  !   do k = sz , ez
  !      do j = sy , ey
  !         do i = sx , ex


  !            if ( y (j) >= 0.0_dp ) then


  !               T (i,j,k) = 545.0_dp / adi % T_ref
  !               P         = 109.5e3_dp / adi % P_ref

  !               Ya (:) = 0.0_dp
  !               Ya (1) = 1.0_dp

  !               call Wmix_i_scalar ( thd , Ya , Ws_i )
  !               rho = P / ( T (i,j,k) * Ws_i )

  !               call ha_scalar ( thd , T (i,j,k) , has )
  !               hm = 0.0_dp
  !               do l = 1 , nrv
  !                  hm = hm + has (l) * Ya (l)
  !               end do

  !               ux = 0.5_dp * ( Sigmau + Deltau * tanh ( y(j) * d2O0_i ) )


  !               v (i,j,k,1) = rho
  !               v (i,j,k,2) = rho * ux
  !               v (i,j,k,3) = rho * vy
  !               v (i,j,k,4) = rho * wz
  !               v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
  !               do l = 1 , nrv+npv+nvv
  !                  v (i,j,k,niv+l) = rho * Ya (l)
  !               end do


  !            else


  !               T (i,j,k) = 1250.0_dp / adi % T_ref
  !               P         = 109.5e3_dp / adi % P_ref

  !               Ya (:)  = 0.0_dp
  !               Ya (2)  = 0.201_dp
  !               Ya (3)  = 0.255_dp
  !               Ya (4)  = 0.544_dp
  !               call X_to_Y ( thd , Ya )

  !               call Wmix_i_scalar ( thd , Ya , Ws_i )
  !               rho = P / ( T (i,j,k) * Ws_i )

  !               call ha_scalar ( thd , T (i,j,k) , has )
  !               hm = 0.0_dp
  !               do l = 1 , nrv
  !                  hm = hm + has (l) * Ya (l)
  !               end do

  !               ux = 0.5_dp * ( Sigmau + Deltau * tanh ( y(j) * d2O0_i ) )


  !               v (i,j,k,1) = rho
  !               v (i,j,k,2) = rho * ux
  !               v (i,j,k,3) = rho * vy
  !               v (i,j,k,4) = rho * wz
  !               v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
  !               do l = 1 , nrv+npv+nvv
  !                  v (i,j,k,niv+l) = rho * Ya (l)
  !               end do


  !            end if


  !         end do
  !      end do
  !   end do


  !   call comm_cons (v)
  !   call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
  !   call upd_extrapolation (v)
  !   do l = 1 , ndim+ndim
  !      call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
  !   end do


  ! end subroutine ini_shear_layer_cheng


!> \brief "Cheng" shear layer IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
!   subroutine ini_shear_layer_cheng_cont ( adi , thd , y , T , W_i , cp , ha , v )


!     type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
!     type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
!     real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


!     integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
!     integer (ip) :: ok , i , j , k , l
!     real (dp)    :: u1_ , u2_ , t1_ , t2_ , p1_ , p2_
!     real (dp)    :: Deltau , Sigmau , d2O0_i , dw0_i
!     real (dp)    :: P , rho , Ws_i , hm , ux , vy , wz
!     real (dp)    :: Ya (nrv+npv+nvv) , has (nrv)


!     allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
!                ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
!                T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                stat = ok )
!     if ( ok > 0 ) call abort_mpi ('error allocate ini')


!     v   = 0.0_dp
!     ha  = 0.0_dp
!     T   = 0.0_dp
!     W_i = 0.0_dp
!     cp  = 0.0_dp


!     u1_ = 2462.74_dp / adi % u_ref
!     u2_ = 1213.96_dp / adi % u_ref

!     t1_ = 875.0_dp  / adi % T_ref
!     t2_ = 1750.0_dp / adi % T_ref

!     p1_ = 109.5e3_dp / adi % P_ref
!     p2_ = 109.5e3_dp / adi % P_ref


!     Deltau = u1_ - u2_
!     Sigmau = u1_ + u2_

!     dw0_i  = 1.839e-4_dp / adi % L_ref
!     d2O0_i = 0.25_dp * dw0_i

!     ! invert
!     d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
!     dw0_i  = 1.0_dp / dw0_i

!     ux = 0.0_dp
!     vy = 0.0_dp
!     wz = 0.0_dp


!     do k = sz , ez
!        do j = sy , ey
!           do i = sx , ex


!              ux        = 0.5_dp * ( Sigmau    + Deltau    * tanh ( y(j) * d2O0_i ) )
!              T (i,j,k) = 0.5_dp * ( (t1_+t2_) + (t1_-t2_) * tanh ( y(j) * d2O0_i ) )
!              P         = 0.5_dp * ( (p1_+p2_) + (p1_-p2_) * tanh ( y(j) * d2O0_i ) )

!              Ya (:)  = 0.0_dp
! !             Ya (2)  = 0.5_dp * ( (0.071_dp) + (0.071_dp) * tanh ( y(j) * d2O0_i ) )
!              Ya (2)  = 0.5_dp * ( (1.000_dp) + (1.000_dp) * tanh ( y(j) * d2O0_i ) )
!              Ya (4)  = 0.5_dp * ( (0.245_dp) - (0.245_dp) * tanh ( y(j) * d2O0_i ) )
!              Ya (6)  = 0.5_dp * ( (0.175_dp) - (0.175_dp) * tanh ( y(j) * d2O0_i ) )
!              Ya (9)  = 0.5_dp * ( (0.580_dp) - (0.580_dp) * tanh ( y(j) * d2O0_i ) )
! !             Ya (9)  = 0.5_dp * ( (1.509_dp) + (0.349_dp) * tanh ( y(j) * d2O0_i ) )
! !             call X_to_Y ( thd , Ya )

!              call Wmix_i_scalar ( thd , Ya , Ws_i )
!              rho = P / ( T (i,j,k) * Ws_i )

!              call ha_scalar ( thd , T (i,j,k) , has )
!              hm = 0.0_dp
!              do l = 1 , nrv
!                 hm = hm + has (l) * Ya (l)
!              end do


!              v (i,j,k,1) = rho
!              v (i,j,k,2) = rho * ux
!              v (i,j,k,3) = rho * vy
!              v (i,j,k,4) = rho * wz
!              v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
!              do l = 1 , nrv+npv+nvv
!                 v (i,j,k,niv+l) = rho * Ya (l)
!              end do


!           end do
!        end do
!     end do


!     call comm_cons (v)
!     call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
!     call upd_extrapolation (v)
!     do l = 1 , ndim+ndim
!        call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
!     end do


!   end subroutine ini_shear_layer_cheng_cont


!> \brief "Miller" shear layer IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine ini_shear_layer_miller ( adi , thd , y , T , W_i , cp , ha , v )


  !   type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
  !   type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
  !   real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


  !   integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
  !   integer (ip) :: ok , i , j , k , l
  !   real (dp)    :: u1_ , u2_ , t1_ , t2_ , p1_ , p2_
  !   real (dp)    :: Deltau , Sigmau , d2O0_i , dw0_i
  !   real (dp)    :: P , rho , Ws_i , hm , ux , vy , wz
  !   real (dp)    :: Ya (nrv+npv+nvv) , has (nrv)


  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini')


  !   v   = 0.0_dp
  !   ha  = 0.0_dp
  !   T   = 0.0_dp
  !   W_i = 0.0_dp
  !   cp  = 0.0_dp


  !   u1_ = 168.00_dp / adi % u_ref
  !   u2_ = 953.11_dp / adi % u_ref

  !   t1_ = 269.0_dp  / adi % T_ref
  !   t2_ = 1475.0_dp / adi % T_ref

  !   p1_ = 94232.25_dp / adi % P_ref
  !   p2_ = 94232.25_dp / adi % P_ref

  !   Deltau = u1_ - u2_
  !   Sigmau = u1_ + u2_

  !   dw0_i  = 4.71e-5_dp / adi % L_ref
  !   d2O0_i = 0.25_dp * dw0_i

  !   ! invert
  !   d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
  !   dw0_i  = 1.0_dp / dw0_i

  !   ux = 0.0_dp
  !   vy = 0.0_dp
  !   wz = 0.0_dp


  !   do k = sz , ez
  !      do j = sy , ey
  !         do i = sx , ex


  !            ux        = 0.5_dp * ( Sigmau    + Deltau    * tanh ( y(j) * d2O0_i ) )
  !            T (i,j,k) = 0.5_dp * ( (t1_+t2_) + (t1_-t2_) * tanh ( y(j) * d2O0_i ) )
  !            P         = 0.5_dp * ( (p1_+p2_) + (p1_-p2_) * tanh ( y(j) * d2O0_i ) )

  !            Ya (:)  = 0.0_dp
  !            Ya (1)  = 0.5_dp * ( (0.100_dp) + (0.100_dp) * tanh ( y(j) * d2O0_i ) )
  !            Ya (2)  = 0.5_dp * ( (0.230_dp) - (0.230_dp) * tanh ( y(j) * d2O0_i ) )
  !            Ya (3)  = 0.5_dp * ( (0.250_dp) - (0.250_dp) * tanh ( y(j) * d2O0_i ) )
  !            Ya (4)  = 0.5_dp * ( (1.420_dp) + (0.380_dp) * tanh ( y(j) * d2O0_i ) )
  !            call X_to_Y ( thd , Ya )

  !            call Wmix_i_scalar ( thd , Ya , Ws_i )
  !            rho = P / ( T (i,j,k) * Ws_i )

  !            call ha_scalar ( thd , T (i,j,k) , has )
  !            hm = 0.0_dp
  !            do l = 1 , nrv
  !               hm = hm + has (l) * Ya (l)
  !            end do


  !            v (i,j,k,1) = rho
  !            v (i,j,k,2) = rho * ux
  !            v (i,j,k,3) = rho * vy
  !            v (i,j,k,4) = rho * wz
  !            v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
  !            do l = 1 , nrv+npv+nvv
  !               v (i,j,k,niv+l) = rho * Ya (l)
  !            end do


  !         end do
  !      end do
  !   end do


  !   call comm_cons (v)
  !   call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
  !   call upd_extrapolation (v)
  !   do l = 1 , ndim+ndim
  !      call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
  !   end do


  ! end subroutine ini_shear_layer_miller


!> \brief "Cheng" & "Miller" combination shear layer IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
!   subroutine ini_shear_layer_mixchengmiller ( inp , adi , thd , y , T , W_i , cp , ha , v )


!     type (inp_type) , intent (in)                                  :: inp  !< input derived type
!     type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
!     type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
!     real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


!     integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
!     integer (ip) :: ok , i , j , k , l
!     real (dp)    :: u1_ , u2_ , t1_ , t2_ , p1_ , p2_
!     real (dp)    :: Deltau , Sigmau , d2O0_i , dw0_i
!     real (dp)    :: P , rho , Ws_i , hm , ux , vy , wz
!     real (dp)    :: Ya (nrv+npv+nvv) , has (nrv)


!     allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
!                ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
!                T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                stat = ok )
!     if ( ok > 0 ) call abort_mpi ('error allocate ini')


!     v   = 0.0_dp
!     ha  = 0.0_dp
!     T   = 0.0_dp
!     W_i = 0.0_dp
!     cp  = 0.0_dp


!     u1_ = 669.14_dp / adi % u_ref
!     u2_ = 2185.35_dp / adi % u_ref

!     t1_ = 545.0_dp  / adi % T_ref
!     t2_ = 1475.0_dp / adi % T_ref

!     p1_ = 94232.25_dp / adi % P_ref
!     p2_ = 94232.25_dp / adi % P_ref

!     Deltau = u1_ - u2_
!     Sigmau = u1_ + u2_

!     dw0_i  = 6.28e-5_dp / adi % L_ref
!     d2O0_i = 0.25_dp * dw0_i

!     ! invert
!     d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
!     dw0_i  = 1.0_dp / dw0_i

!     ux = 0.0_dp
!     vy = 0.0_dp
!     wz = 0.0_dp


!     do k = sz , ez
!        do j = sy , ey
!           do i = sx , ex


!              ux        = 0.5_dp * ( Sigmau    + Deltau    * tanh ( y(j) * d2O0_i ) )
!              T (i,j,k) = 0.5_dp * ( (t1_+t2_) + (t1_-t2_) * tanh ( y(j) * d2O0_i ) )
!              P         = 0.5_dp * ( (p1_+p2_) + (p1_-p2_) * tanh ( y(j) * d2O0_i ) )

!              do l = 1,nrv-1
!                 Ya (l) = 0.5_dp * ( inp % Y1f(l)+inp % Y0o(l) +                           &
!                                   (inp % Y1f(l)-inp % Y0o(l)) * tanh ( y(j) * d2O0_i ) )
!              end do
!              Ya (nrv) = 1.0_dp - sum ( Ya (1:nrv-1) )
!              do l = nrv+npv,nv-niv
!                 Ya (l) = 0.5_dp * ( inp % Y1f(l)+inp % Y0o(l) +                           &
!                                   (inp % Y1f(l)-inp % Y0o(l)) * tanh ( y(j) * d2O0_i ) )
!              end do
!              if ( npv > 0 ) then
!                 do l = nrv+1,nrv+npv
!                    Ya (l) = 0.5_dp * ( inp % Y1f(l)+inp % Y0o(l) +                        &
!                                      (inp % Y1f(l)-inp % Y0o(l)) * tanh ( y(j) * d2O0_i ) )
!                 end do
!              end if
! !             call X_to_Y ( thd , Ya )

!              call Wmix_i_scalar ( thd , Ya , Ws_i )
!              rho = P / ( T (i,j,k) * Ws_i )

!              call ha_scalar ( thd , T (i,j,k) , has )
!              hm = 0.0_dp
!              do l = 1 , nrv
!                 hm = hm + has (l) * Ya (l)
!              end do


!              v (i,j,k,1) = rho
!              v (i,j,k,2) = rho * ux
!              v (i,j,k,3) = rho * vy
!              v (i,j,k,4) = rho * wz
!              v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
!              do l = 1 , nrv+npv+nvv
!                 v (i,j,k,niv+l) = rho * Ya (l)
!              end do


!           end do
!        end do
!     end do


!     call comm_cons (v)
!     call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
!     call upd_extrapolation (v)
!     do l = 1 , ndim+ndim
!        call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
!     end do


!   end subroutine ini_shear_layer_mixchengmiller


!> \brief Jet shear layer IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine ini_jet ( adi , thd , T , W_i , cp , ha , v )


  !   type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
  !   type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


  !   integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
  !   integer (ip) :: ok , i , j , k , l
  !   real (dp)    :: Ya (nrv+npv+nvv) , has (nrv)
  !   real (dp)    :: u1 , v1 , wz
  !   real (dp)    :: P , rho , Ws_i , hm


  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )        , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )        , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )        , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini')


  !   v   = 0.0_dp
  !   ha  = 0.0_dp
  !   T   = 0.0_dp
  !   W_i = 0.0_dp
  !   cp  = 0.0_dp

  !   u1  = 20.0_dp / adi % u_ref
  !   v1  = 0.00_dp / adi % u_ref
  !   wz  = 0.00_dp / adi % u_ref


  !   do k = sz , ez
  !      do j = sy , ey
  !         do i = sx , ex

  !            T (i,j,k)= 300.0_dp    / adi % T_ref
  !            P        = 101325.0_dp / adi % P_ref

  !            Ya (:)   = 0.0_dp
  !            Ya (1)   = 0.233_dp
  !            Ya (2)   = 1.0_dp - Ya (1)

  !            ! call X_to_Y ( thd , Ya )

  !            call Wmix_i_scalar ( thd , Ya , Ws_i )
  !            rho = P / ( T(i,j,k) * Ws_i )

  !            call ha_scalar ( thd , T (i,j,k) , has )
  !            hm = 0.0_dp
  !            do l = 1 , nrv
  !               hm = hm + has (l) * Ya (l)
  !            end do

  !            v (i,j,k,1) = rho
  !            v (i,j,k,2) = rho * u1
  !            v (i,j,k,3) = rho * v1
  !            v (i,j,k,4) = rho * wz
  !            v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( u1*u1 + v1*v1 + wz*wz )
  !            do l = 1 , nrv+npv+nvv
  !               v (i,j,k,niv+l) = rho * Ya (l)
  !            end do


  !         end do
  !      end do
  !   end do


  !   call comm_cons (v)
  !   call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
  !   call upd_extrapolation (v)
  !   do l = 1 , ndim+ndim
  !      call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
  !   end do


  !  end subroutine ini_jet


!> \brief Premixed IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine ini_premix ( adi , thd , T , W_i , cp , ha , v )


  !   type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
  !   type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


  !   integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
  !   integer (ip) :: ok , i , l
  !   real (dp)    :: P , dummy , rho , Ws_i , hm , ux , vy , wz
  !   real (dp)    :: Ya (nrv+npv+nvv) , has (nrv)


  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini')


  !   v   = 0.0_dp
  !   ha  = 0.0_dp
  !   T   = 0.0_dp
  !   W_i = 0.0_dp
  !   cp  = 0.0_dp


  !   P = thd % P0 / adi % P_ref
  !   ux = 0.0_dp ; vy = 0.0_dp ; wz = 0.0_dp


  !   open ( unit = 501 , file = 'post_premix.tmp' , status = 'old' , iostat = ok )
  !   if ( ok /= 0 ) then
  !      call abort_mpi ('error opening ')
  !   end if
  !   read (501,*)
  !   do i = sx , ex
  !      read (501,*) dummy , T (i,sy,sz) , dummy , dummy , v (i,sy,sz,2)
  !      T (i,sy,sz) = T (i,sy,sz) / adi % T_ref
  !      v (i,sy,sz,2) = v (i,sy,sz,2) / ( 100.0_dp * adi % u_ref )
  !   end do
  !   close (501)

  !   open ( unit = 501 , file = 'post_premix.spc' , status = 'old' , iostat = ok )
  !   if ( ok /= 0 ) then
  !      call abort_mpi ('error opening ')
  !   end if
  !   read (501,*)
  !   read (501,*)
  !   read (501,*)
  !   do i = sx , ex
  !      read (501,*) dummy , dummy , v (i,sy,sz,niv+1:niv+nrv)
  !   end do
  !   close (501)


  !   do i = sx , ex

  !      do l = 1,nrv-1
  !         Ya (l) = v (i,sy,sz,niv+l)
  !      end do
  !      Ya (nrv) = 1.0_dp - sum ( Ya (1:nrv-1) )
  !      ! call X_to_Y ( thd , Ya )

  !      call Wmix_i_scalar ( thd , Ya , Ws_i )
  !      rho = P / ( T (i,sy,sz) * Ws_i )

  !      call ha_scalar ( thd , T (i,sy,sz) , has )
  !      hm = 0.0_dp
  !      do l = 1 , nrv
  !         hm = hm + has (l) * Ya (l)
  !      end do

  !      ux = v (i,sy,sz,2)

  !      T (i,sy:ey,sz:ez) = T (i,sy,sz)

  !      v (i,sy:ey,sz:ez,1) = rho
  !      v (i,sy:ey,sz:ez,2) = rho * ux
  !      v (i,sy:ey,sz:ez,3) = rho * vy
  !      v (i,sy:ey,sz:ez,4) = rho * wz
  !      v (i,sy:ey,sz:ez,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
  !      do l = 1 , nrv+npv+nvv
  !         v (i,sy:ey,sz:ez,niv+l) = rho * Ya (l)
  !      end do

  !   end do


  !   call comm_cons (v)
  !   call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
  !   call upd_extrapolation (v)
  !   do l = 1 , ndim+ndim
  !      call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
  !   end do


  ! end subroutine ini_premix


!> \brief "Cheng" & "Miller" combination diffusion IC.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
!   subroutine ini_diffusion_mixchengmiller ( inp , adi , thd , x , T , W_i , cp , ha , v )


!     type (inp_type) , intent (in)                                  :: inp  !< input derived type
!     type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
!     type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
!     real (dp) , allocatable , dimension (:) , intent (in)          :: x    !< x-coordinate array
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
!     real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


!     integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
!     integer (ip) :: ok , i , j , k , l
!     real (dp)    :: u1_ , u2_ , t1_ , t2_ , p1_ , p2_
!     real (dp)    :: Deltau , Sigmau , d2O0_i , dw0_i
!     real (dp)    :: P , rho , Ws_i , hm , ux , vy , wz
!     real (dp)    :: Ya (nrv+npv+nvv) , has (nrv)


!     allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
!                ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
!                T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
!                stat = ok )
!     if ( ok > 0 ) call abort_mpi ('error allocate ini')


!     v   = 0.0_dp
!     ha  = 0.0_dp
!     T   = 0.0_dp
!     W_i = 0.0_dp
!     cp  = 0.0_dp


!     u1_ = 669.14_dp  / adi % u_ref
!     u2_ = 1151.57_dp / adi % u_ref

!     t1_ = 545.0_dp  / adi % T_ref
!     t2_ = 1475.0_dp / adi % T_ref

!     p1_ = 94232.25_dp / adi % P_ref
!     p2_ = 94232.25_dp / adi % P_ref

!     Deltau = u1_ - u2_
!     Sigmau = u1_ + u2_

!     dw0_i  = 1.975e-4_dp / adi % L_ref
!     d2O0_i = 0.25_dp * dw0_i

!     ! invert
!     d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
!     dw0_i  = 1.0_dp / dw0_i

!     ux = 0.0_dp
!     vy = 0.0_dp
!     wz = 0.0_dp


!     do k = sz , ez
!        do j = sy , ey
!           do i = sx , ex


! !             ux        = 0.5_dp * ( Sigmau    + Deltau    * tanh ( x(i) * d2O0_i ) )
!              T (i,j,k) = 0.5_dp * ( (t1_+t2_) + (t1_-t2_) * tanh ( x(i) * d2O0_i ) )
!              P         = 0.5_dp * ( (p1_+p2_) + (p1_-p2_) * tanh ( x(i) * d2O0_i ) )

!              do l = 1,nrv-1
!                 Ya (l) = 0.5_dp * ( inp % Y1f(l)+inp % Y0o(l) +                           &
!                                    (inp % Y1f(l)-inp % Y0o(l)) * tanh ( x(i) * d2O0_i ) )
!              end do
!              Ya (nrv) = 1.0_dp - sum ( Ya (1:nrv-1) )
! !             call X_to_Y ( thd , Ya )

!              call Wmix_i_scalar ( thd , Ya , Ws_i )
!              rho = P / ( T (i,j,k) * Ws_i )

!              call ha_scalar ( thd , T (i,j,k) , has )
!              hm = 0.0_dp
!              do l = 1 , nrv
!                 hm = hm + has (l) * Ya (l)
!              end do


!              v (i,j,k,1) = rho
!              v (i,j,k,2) = rho * ux
!              v (i,j,k,3) = rho * vy
!              v (i,j,k,4) = rho * wz
!              v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
!              do l = 1 , nrv+npv+nvv
!                 v (i,j,k,niv+l) = rho * Ya (l)
!              end do


!           end do
!        end do
!     end do


!     call comm_cons (v)
!     call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
!     call upd_extrapolation (v)
!     do l = 1 , ndim+ndim
!        call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
!     end do


!   end subroutine ini_diffusion_mixchengmiller


!> \brief Detonation IC.
!!
!! Detonation mix at ambient conditions (T=300K, P=1atm).
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine ini_detonation ( inp , adi , thd , x , y , z , T , W_i , cp , ha , v )

  !   type (inp_type)                               , intent (in)    :: inp
  !   type (adi_type)                               , intent (in)    :: adi
  !   type (thd_type)                               , intent (in)    :: thd
  !   real (dp) , allocatable , dimension (:)       , intent (in)    :: x
  !   real (dp) , allocatable , dimension (:)       , intent (in)    :: y
  !   real (dp) , allocatable , dimension (:)       , intent (in)    :: z
  !   real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: T
  !   real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: W_i
  !   real (dp) , allocatable , dimension (:,:,:)   , intent (inout) :: cp
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha
  !   real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v

  !   integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
  !   integer (ip)                                           :: ok , i , j , k , l
  !   real    (dp)                                           :: Ws_i

  !   real (dp) :: Ya (nrv+npv) , has (nrv)
  !   real (dp) :: ux, vy, wz, Pamb, Tamb, rho, hm, cp0, gamma0
  !   real (dp) :: ux_inf , wrky , wrkz , eps , noise
  !   real (dp) :: alpha , beta

  !   allocate ( v    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
  !              ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
  !              T    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              W_i  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              cp   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
  !              stat = ok )
  !   if ( ok > 0 ) then
  !      write (*,*) 'error allocate ini_amb'
  !      call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
  !   end if

  !   v   = 0.0_dp
  !   ha  = 0.0_dp
  !   T   = 0.0_dp
  !   W_i = 0.0_dp
  !   cp  = 0.0_dp

  !   ! Application of (T, p, rho, u, v, w, hm) in the calculation domain
  !   do k = sz , ez ! z-direction
  !      do j = sy , ey ! y-direction
  !         do i = sx , ex ! x-direction

  !            if ( x (i) * adi % L_ref < 5.0e-3_dp ) then

  !               !H2 H O2 O OH HO2 H2O2 H2O AR N2
  !               !Ya (:)  = 0.000_dp ! Initialization
  !               !Ya (2)  = 0.296_dp ! H2
  !               !Ya (4)  = 0.148_dp ! O2
  !               !Ya (10) = 0.556_dp ! N2

  !               !H2 H O2 O OH HO2 H2O2 H2O AR
  !               Ya (:) = 0.000_dp ! Initialization
  !               Ya (1) = 0.200_dp ! H2
  !               Ya (3) = 0.100_dp ! O2
  !               Ya (9) = 0.700_dp ! AR


                
  !               call X_to_Y ( thd , Ya )

  !               ux = 0.0_dp
  !               vy = 0.0_dp
  !               wz = 0.0_dp

  !               Pamb = 25.0_dp * 101325.0_dp / adi % P_ref
  !               Tamb =             2000.0_dp / adi % T_ref

  !               call Wmix_i_scalar ( thd , Ya , Ws_i )
  !               call ha_scalar ( thd , Tamb , has )

  !               hm = 0.0_dp
  !               do l = 1 , nrv
  !                  hm = hm + has (l) * Ya (l)
  !               end do

  !               rho = Pamb / ( Tamb * Ws_i )


  !               T (i,j,k)   = Tamb ! initial guess for the temperature (mandatory)

  !               v (i,j,k,1) = rho
  !               v (i,j,k,2) = rho * ux
  !               v (i,j,k,3) = rho * vy
  !               v (i,j,k,4) = rho * wz
  !               v (i,j,k,5) = rho * hm - Pamb + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )

  !               do l = 1 , nrv+npv
  !                  v (i,j,k,niv+l) = rho * Ya (l)
  !               end do

  !            else

  !               !Ya (:)  = 0.000_dp ! Initialization
  !               !Ya (2)  = 0.296_dp ! H2
  !               !Ya (4)  = 0.148_dp ! O2
  !               !Ya (10) = 0.556_dp ! N2

  !               Ya (:) = 0.000_dp ! Initialization
  !               Ya (1) = 0.200_dp ! H2
  !               Ya (3) = 0.100_dp ! O2
  !               Ya (9) = 0.700_dp ! AR


                
  !               call X_to_Y ( thd , Ya )

  !               ux = 0.0_dp
  !               vy = 0.0_dp
  !               wz = 0.0_dp

  !               Pamb = 101325.0_dp / adi % P_ref
  !               Tamb =    298.0_dp / adi % T_ref

  !               call Wmix_i_scalar ( thd , Ya , Ws_i )
  !               call ha_scalar ( thd , tamb , has )

  !               hm = 0.0_dp
  !               do l = 1 , nrv
  !                  hm = hm + has (l) * Ya (l)
  !               end do

  !               rho = Pamb / ( Tamb * Ws_i )


  !               T (i,j,k)   = Tamb ! initial guess for the temperature (mandatory)

  !               v (i,j,k,1) = rho
  !               v (i,j,k,2) = rho * ux
  !               v (i,j,k,3) = rho * vy
  !               v (i,j,k,4) = rho * wz
  !               v (i,j,k,5) = rho * hm - Pamb + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )

  !               do l = 1 , nrv+npv
  !                  v (i,j,k,niv+l) = rho * Ya (l)
  !               end do

  !            end if

  !         end do
  !      end do
  !   end do

        
  !   call comm_cons (v)
  !   call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
  !   call upd_extrapolation (v)
  !   do l = 1 , ndim+ndim
  !      call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
  !   end do
    
  ! end subroutine ini_detonation


  

!> \brief Klein random field for Digital filter IC.
!!
!! Create a random field for Synthetic Turbulent Generator (STG) involve
!! in equation 16 of Klein et al. article.
!!
!! References :
!! -# M. Klein, A. Sadiki and J. Janicka
!!    "A digital filter based generation of inflow data for spatially developing
!!     direct numerical or large eddy simulation", 
!!    Journal of computational physics, 2003
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
!   subroutine init_klein_random_field ( inp )


!     type (inp_type) , intent (inout)   :: inp  !< input derived type

!     logical          :: loop
!     integer (ip)     :: ok , i , j , k , l
!     real (dp)        :: epsv , factor_fw
!     real (dp)        :: nxls_i , nyls_i , nzls_i


!     ! Test if the neigh (i) == MPI_PROC_NULL and if BC for STG is needed
!     ! This to avoid allocating useless arrays (inner process)
!     !===================================================================
!     loop = .false.
!     call apply_stg ( loop )
!     if ( .not. loop ) return ! leave the subroutine


!     ! Length scale parameter between [2,..,100]
!     inp % ext % nxls = 10
!     inp % ext % nyls = inp % ext % nxls
!     inp % ext % nzls = inp % ext % nxls ! Homogeneous length scale

!     nxls_i = 1.0_dp / float (inp % ext % nxls)
!     nyls_i = 1.0_dp / float (inp % ext % nyls)
!     nzls_i = 1.0_dp / float (inp % ext % nzls)

!     ! Gaussian filter width such that nxfw >= 2*nxls
!     factor_fw = 2 ! factor filter width >=2
!     inp % ext % nxfw = factor_fw * inp % ext % nxls
!     inp % ext % nyfw = inp % ext % nxfw
!     inp % ext % nzfw = inp % ext % nxfw ! Homogeneous filter width

!     nxfw = inp % ext % nxfw
!     nyfw = inp % ext % nyfw
!     nzfw = inp % ext % nzfw


!     ! Initialization of ndim random fields

!     !===========================================================================
!     !==== Revoir peut etre la parallelisation du champ initial turbulent =======
!     ! avec call MPI_BCAST ( wrk , 2 , MPI_INTEGER , 0 , MPI_COMM_WORLD , mpicode ) ????
!     !===========================================================================

! !    allocate ( inp % ext % rand  ( -nxfw:nxfw , -nyfw+1:nyfw+nty , -nzfw+1:nzfw+ntz , 1:ndim ) , &

!     allocate ( inp % ext % rand  ( -nxfw:nxfw , -nyfw+sy:nyfw+ey , -nzfw+sz:nzfw+ez , 1:ndim ) , &
!                stat = ok)
!     if ( ok > 0 ) call abort_mpi ('error allocate ini inp % ext % rand')

!     do l = 1 , ndim
!        do k = -nzfw+sz , nzfw+ez
!           do j = -nyfw+sy , nyfw+ey
!              do i = -nxfw , nxfw
!                 call random_number (epsv)
!                 inp % ext % rand (i,j,k,l) = epsv + epsv - 1.0_dp
!              end do
!           end do
!        end do
!     end do


!     ! Initialization of 3D filter coefficients (eq.11)
!     allocate ( inp % ext % bijk  ( -nxfw:nxfw , -nyfw:nyfw , -nzfw:nzfw ) , &
!                stat = ok)
!     if ( ok > 0 ) call abort_mpi ('error allocate ini inp % ext % bijk')

!     call ini_3D_convolution_filter ( nxls_i , nyls_i , nzls_i , inp % ext % bijk )


!   end subroutine init_klein_random_field


!> \brief Klein 3D filter IC.
!!
!! Create 3D filter for Synthetic Turbulent Generator (STG), by multipling 
!! the 3 filter coefficients. Refered to eq.11 and 14 in the article.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine ini_3D_convolution_filter ( nxls_i , nyls_i , nzls_i , bijk )


  !   real (dp)                                                      :: nxls_i , nyls_i , nzls_i  !< inverse of the length scale parameter
  !   real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: bijk                      !< 3D digital filter


  !   real (dp) , allocatable , dimension (:)        :: bi , bj , bk

  !   integer (ip)  :: ok , i , j , k

  !   allocate ( bi ( -nxfw:nxfw ) , &
  !              bj ( -nyfw:nyfw ) , & 
  !              bk ( -nzfw:nzfw ) , &
  !              stat = ok )
  !   if ( ok > 0 ) call abort_mpi ('error allocate ini bi, bj, bk')


  !   ! 1D coefficient filters (eq.14)
  !   call filter_coef ( nxfw , nxls_i , bi )
  !   call filter_coef ( nyfw , nyls_i , bj )
  !   call filter_coef ( nzfw , nzls_i , bk )

  !   ! 3D filters (eq.11)
  !   do k = -nzfw , nzfw
  !      do j = -nyfw , nyfw
  !         do i = -nxfw , nxfw
  !            bijk (i,j,k) = bi (i) * bj (j) * bk (k)
  !         end do
  !      end do
  !   end do


  !   deallocate ( bi , bj , bk )


  ! end subroutine ini_3D_convolution_filter



!> \brief Calculate filter coefficient.
!!
!! Filter coefficient for synthetic turbulent generator (STG).
!! Defined by Klein et al. in eq.14 of their article.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  ! subroutine filter_coef ( nfw , nls_i , bi )


  !   integer (ip)                            , intent (in)    :: nfw   !< filter width such that nfw >= 2/nls_i
  !   real (dp)                               , intent (in)    :: nls_i !< inverse of the length scale
  !   real (dp) , allocatable , dimension (:) , intent (inout) :: bi    !< filter coefficient

  !   integer (ip)              :: i

  !   real (dp) , parameter     :: pi05      = 0.5_dp * pi
  !   real (dp)                 :: sum_b , wrk


  !   sum_b = 0.0_dp
  !   do i = -nfw , nfw
  !      wrk = i * nls_i
  !      wrk = exp ( -pi05 * wrk * wrk )
  !      sum_b = sum_b + wrk * wrk
  !   end do
  !   sum_b = sqrt (sum_b)

  !   do i = -nfw , nfw
  !      wrk = i * nls_i
  !      bi (i) = exp ( -pi05 * wrk * wrk ) / sum_b
  !   end do


  ! end subroutine filter_coef


!> \brief Initialized perturbation condition.
!!
!! Initialized variables to generate blowing and suction at specific boundaries
!! to help transition frome laminar to turbulent flow
!!
!! References :
!! -# M.F. Shahab, Etude numerique de l'influence de l'impact d'une onde 
!!    de choc et d'un transfert de chaleur sur une couche limite en developpement,
!!    PhD thesis, ENSMA, 2011
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
!   subroutine init_perturbation ( inp )


!     type (inp_type) , intent (inout)   :: inp  !< input derived type

!     logical          :: loop
!     integer (ip)     :: ok , i , j , k , l
!     real (dp)        :: factor , sumf
!     real (dp) , dimension (:) , allocatable   :: farray


!     ! Test if the neigh (i) == MPI_PROC_NULL and if BC for perturbation is needed
!     ! This to avoid allocating useless arrays (inner process)
!     loop = .false.
!     do i = 1 , nneigh
!        if ( neigh (i) == MPI_PROC_NULL .and. &
!             (    bc (S) == wallperturbation  &
! !            .or. bc (S) == wallinjection     &
!             .or. bc (W) == syntheticturb     ) ) loop = .true. ; exit ! Apply, on SOUTH face here.
!     end do
!     if ( .not. loop ) return ! leave the subroutine


!     ! Summation limit
!     inp % ptb % ylmax = 20
!     inp % ptb % zlmax = 10
!     inp % ptb % mmax = 5

!     ! Reduction factor
!     ! such as Z(l+1) = factor * Z(l) and sum(Z(1:zlmax)) = 1
!     factor = 1.0_dp / 1.25_dp


!     ! Initialized Yl array for y-direction function
!     allocate ( farray            ( inp % ptb % ylmax ) , &
!                inp % ptb % yl    ( inp % ptb % ylmax ) , &
!                stat = ok)
!     if ( ok > 0 ) call abort_mpi ('error allocate ini inp % ptb % yl')

!     sumf = 0.0_dp
!     do l = 1 , inp % ptb % ylmax
!        farray (l) = factor ** (l-1)
!        sumf       = sumf + farray (l)
!     end do
!     inp % ptb % yl (1) = 1.0_dp / sumf
!     do l = 2 , inp % ptb % ylmax
!        inp % ptb % yl (l) = farray (l) * inp % ptb % yl (1)
!     end do

!     deallocate ( farray )


!     ! Initialized Zl array for z-direction function
!     allocate ( farray            ( inp % ptb % zlmax ) , &
!                inp % ptb % zl    ( inp % ptb % zlmax ) , &
!                stat = ok)
!     if ( ok > 0 ) call abort_mpi ('error allocate ini inp % ptb % zl')

!     sumf = 0.0_dp
!     do l = 1 , inp % ptb % zlmax
!        farray (l) = factor ** (l-1)
!        sumf       = sumf + farray (l)
!     end do
!     inp % ptb % zl (1) = 1.0_dp / sumf
!     do l = 2 , inp % ptb % zlmax
!        inp % ptb % zl (l) = farray (l) * inp % ptb % zl (1)
!     end do

!     deallocate ( farray )


! !    ! Initialized Tm array for time function
! !    allocate ( farray            ( inp % ptb % mmax ) , &
! !               inp % ptb % tm    ( inp % ptb % mmax ) , &
! !               stat = ok)
! !    if ( ok > 0 ) call abort_mpi ('error allocate ini inp % ptb % tm')

! !    sumf = 0.0_dp
! !    do l = 1 , inp % ptb % mmax
! !       farray (l) = factor ** (l-1)
! !       sumf       = sumf + farray (l)
! !    end do
! !    inp % ptb % tm (1) = 1.0_dp / sumf
! !    do l = 2 , inp % ptb % mmax
! !       inp % ptb % tm (l) = farray (l) * inp % ptb % tm (1)
! !    end do

! !    deallocate ( farray )


!   end subroutine init_perturbation


end module ICs

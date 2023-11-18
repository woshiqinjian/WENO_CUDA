!------------------------------------------------------------------------------
! MODULE: parameters
!------------------------------------------------------------------------------
!> \brief Common parameters.
!!
!! This module contains all the commons parameters for the simulation
!! and the post-treatment.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
module parameters

  implicit none

  integer ( kind (0) ) , parameter   :: dp = kind ( 0.0d0 ) , & !< double precision
                                        sp = kind ( 0.0   ) , & !< simple precision
                                        ip = kind ( 0     )     !< integer precision

  real (dp) , parameter              :: pi = acos(-1.0_dp) !< pi definition

  integer (ip) , parameter :: ndimmax      = 3  !< number of physical dimensions
  integer (ip) , parameter :: niv          = 5  !< number of conservative variable
  integer (ip) , parameter :: ng           = 4  !< number of ghost points for WENO scheme
  integer (ip) , parameter :: len_default  = 50 !< default character length
  integer (ip) , parameter :: nplanesmax   = 20 !< maximum number of plane for stats
  integer (ip) , parameter :: rank_default = 0  !< default processor number

  integer (ip) :: nbit        !< nbits for the access direct files
  logical      :: ind_files   !< use individual binary files

  real (dp)    :: rho_0 , delta_U , delta_0  
  real (dp)    :: q0
  real (dp)    :: m_trg
  real (dp) , parameter :: f_ini = 0.1_dp

  real    (dp) :: Lx , Ly , Lz          !< length of the whole domain
  integer (ip) :: ntx , nty , ntz       !< number of points of the whole domain
  integer (ip) :: nxmax , nymax , nzmax !< maximum number of points in one direction among all the process
  integer (ip) :: sx , ex , &
                  sy , ey , &
                  sz , ez !< start/end points per subdomain
  integer (ip) :: ndim      !< CURRENT number of dimensions
  integer (ip) :: nv        !< CURRENT number of total variables = niv + nrv + npv + nvv
  integer (ip) :: nderiv    !< nderiv = ndim * ( ndim + 1 ) + ndim * ( nrv + npv + nvv ) number of derivative variables to communicate through MPI
                            !<        if LES : + ndim ( nvv + 1 ) + 1
  integer (ip) :: nderivSGS !< nderivSGS = ndim derivatives of the isotropic SGS tensor

  integer (ip) :: nxfw , nyfw , nzfw !< Klein filter width (such that nxfw >= 2*nxls)
  
  logical      :: weno_avg , vis , vis_dt , reaction , filter_reaction , eglib , LES , forcesum !< common logical variables
  real    (dp) :: CFL , CFLmin , Fo , CPR , max_rel_weight , percent_weight_SGS , zst !< common real variables

  integer (ip) :: idum !< initial value for Random Number Generators (Numerical Recipes)
  real    (dp) , dimension (:,:) , allocatable :: duidxj , kdelta

  !> variables related to I/O files
  integer , dimension (nplanesmax,ndimmax)  :: coord_vol_min , coord_vol_max
  logical , dimension (nplanesmax)          :: log_volume , log_planeXY , log_planeXZ , log_planeYZ
  integer (ip)                              :: nbins
  integer (ip)                              :: ntimes

  integer (ip) , parameter                  :: unit_adi       = 21  , &
                                               unit_grid      = 22  , &
                                               unit_inp       = 23  , &
                                               unit_ext       = 24  , &
                                               unit_HIT       = 25  , &
                                               unit_restart   = 26  , &
                                               unit_time_rest = 27  , &
                                               unit_time_stat = 28  , &
                                               unit_post      = 29  , &
                                               unit_plot      = 30  , &
                                               unit_stat      = 31  , &
                                               unit_spec      = 32  , &
                                               unit_stop      = 33  , &
                                               unit_tabchem   = 34  , &
                                               unit_KE        = 35  , &
                                               unit_moment    = 36  , &
                                               unit_skew      = 37  , &
                                               unit_flat      = 38  , &
                                               unit_rms       = 39  , &
                                               unit_mach      = 40  , &
                                               unit_mean      = 41  , &
                                               unit_theta     = 42  , &
                                               unit_ptheta    = 43  , &
                                               unit_taylor    = 44  , &
                                               unit_diss      = 45  , &
                                               unit_diss2     = 46  , &
                                               unit_eta       = 47  , &
                                               unit_LES       = 48  , &
                                               unit_statvol   = 100 , & ! to avoid collisions
                                               unit_statXY    = 200 , & ! to avoid collisions
                                               unit_statXZ    = 300 , & ! to avoid collisions
                                               unit_statYZ    = 400 , & ! to avoid collisions
                                               unit_probes    = 500 , & ! to avoid collisions
                                               unit_conv      = 600 , & ! to avoid collisions
                                               unit_pdf       = 700 , & ! to avoid collisions
                                               unit_cavg      = 800     ! to avoid collisions

  character (len_default) , parameter       :: file_adi       = 'todo.adi'        , &
                                               file_thd       = 'todo.thd'        , &
                                               file_tranfit   = 'tranfit.out'     , &
                                               file_grid      = 'grid.dat'        , &
                                               file_inp       = 'input.dat'       , &
                                               file_syntturb  = 'syntturb.dat'    , &
                                               file_tabchem   = 'tabchem.dat'     , &
!                                               file_saveascii = 'save.out'        , & ! old variable NOT USED anymore
                                               file_restart   = 'restart'         , &
                                               file_statvol   = 'volume'          , &
                                               file_statXY    = 'statXY'          , &
                                               file_statXZ    = 'statXZ'          , &
                                               file_statYZ    = 'statYZ'          , &
                                               file_time_rest = 'times_rest.out ' , &
                                               file_time_stat = 'times_stat.out ' , &
                                               file_plXY      = 'planesXY.out'    , &
                                               file_plXZ      = 'planesXZ.out'    , &
                                               file_plYZ      = 'planesYZ.out'    , &
                                               file_post      = 'post.dat'        , &
                                               file_plot      = 'plot'            , &
                                               file_inst      = 'inst'            , &
                                               file_reyavg    = 'reyavg'          , &
                                               file_favavg    = 'favavg'          , &
                                               file_stat      = 'stat'            , &
                                               file_conv      = 'convergence'     , &
                                               file_pdf       = 'pdfs'            , &
                                               file_spec      = 'spectra'         , &
                                               file_cavg      = 'cond_avgs'       , &
                                               file_KE        = 'KE_trace.dat'    , &
                                               file_moment    = 'moment.dat'      , &
                                               file_skew      = 'skewness.dat'    , &
                                               file_flat      = 'flatness.dat'    , &
                                               file_rms       = 'rms_values.dat'  , &
                                               file_mach      = 'mach_value.dat'  , &
                                               file_mean      = 'mean_values.dat' , &
                                               file_theta     = 'theta.dat'       , &
                                               file_ptheta    = 'ptheta.dat'      , &
                                               file_taylor    = 'taylor.dat'      , &
                                               file_diss      = 'diss.dat'        , &
                                               file_diss2     = 'diss2.dat'       , &
                                               file_eta       = 'kolmogorov.dat'  , &
                                               file_LES       = 'SGS_value.dat'   , &                                               
                                               file_stop      = 'STOP'

  character (len_default) , parameter       :: dir_parent     = './'              , &
                                               dir_restart    = 'restarts/'       , &
                                               dir_statvol    = 'statsvol/'       , &
                                               dir_statXY     = 'statsXY/'        , &
                                               dir_statXZ     = 'statsXZ/'        , &
                                               dir_statYZ     = 'statsYZ/'        , &
                                               dir_probes     = 'probes/'         , &
                                               dir_plot       = 'plots/'          , &
                                               dir_statvar    = 'statvar/'

  character (len_default) , parameter       :: format_restart = '(I5.5)'          , &
                                               format_nplane  = '(I2.2)'

  !> Type of initial conditions
  character (len_default) , parameter       :: Init_rank            = 'Init_rank'            , &
                                               Diffusion            = 'Diffusion'            , &
                                               Bogey                = 'Bogey'                , & ! shear layer
                                               Fu                   = 'Fu'                   , & ! shear layer
                                               Miller               = 'Miller'               , & ! shear layer
                                               MIXChengMiller       = 'MIXChengMiller'       , &
                                               Premix               = 'Premix'               , &
                                               Jet                  = 'Jet'                  , &
                                               DIFFMIXChengMiller   = 'DIFFMIXChengMiller'   , &
                                               Fedkiw               = 'Fedkiw'               , & ! shock tube
                                               DMR                  = 'DMR'                  , & ! double Mach reflection
                                               HIT                  = 'HIT'                  , & ! Decaying HIT
                                               Ambient              = 'Ambient'              , & ! Ambient conditions
                                               Poiseuille           = 'Poiseuille'           , & ! Poiseuille flow
                                               Vortex               = 'Vortex'               , & ! Isentropic vortex
                                               vortexring           = 'vortexring'           , & ! Isentropic vortex ring
                                               shear                = 'shearlayer'           , & ! Shear layer
                                               channel              = 'channel'              , & ! Channel flow
                                               Detonation           = 'Detonation'               ! Detonation wave

  !> Type of boundary conditions
  character (len_default) , parameter       :: inner                = 'inner'                , &
                                               symmetryplane        = 'symmetryplane'        , &
                                               periodic             = 'periodic'             , &
                                               extrapolation        = 'extrapolation'        , &
                                               inflow               = 'inflow'               , &
                                               noreflection         = 'noreflection'         , &
                                               adiabaticwall        = 'adiabaticwall'        , &
                                               slipwall             = 'slipwall'             , &
                                               movingwall           = 'movingwall'           , &
                                               syntheticturb        = 'syntheticturb'        , &
                                               shock                = 'shock'                , &
                                               postshock_slipwall   = 'postshock_slipwall'   , &
                                               inflowjet            = 'inflowjet'            , &
                                               wallinjection        = 'wallinjection'        , &
                                               wallperturbation     = 'wallperturbation'     , &
                                               inflowcheng          = 'inflowcheng'          , &
                                               inflowmiller         = 'inflowmiller'         , &
                                               inflowmixchengmiller = 'inflowmixchengmiller' , &
                                               inflowbogey          = 'inflowbogey'          , &
                                               premixfix            = 'premixfix'            , &
                                               nscbc_out            = 'nscbc_out'            , &
                                               supersonicflow       = 'supersonicflow'

  !> Combustion models
  character (len_default) , parameter       :: PSR                  = 'PSR'                  , &
                                               MIL                  = 'MIL'                  , &
                                               hybrid_PSR_MIL       = 'hybrid_PSR_MIL'

  !> SGS turbulent models
  character (len_default) , parameter       :: Smagorinsky          = 'Smagorinsky'          , &
                                               StructureFunction    = 'StructureFunction'    , &
                                               WALE                 = 'WALE'                 , &
                                               Yoshizawa            = 'Yoshizawa'

  !> Post-treatment data
  character (len_default) , parameter       :: rey_avg              = 'rey_avg' , &
                                               fav_avg              = 'fav_avg' , &
                                               decimal              = 'decimal' , &
                                               logarit              = 'logarit'

  !> Min/max bounds for physical quantities
  real (dp) , parameter                     :: min_Y    = 0.0_dp    , max_Y    = 1.0_dp   , &
                                               min_Z    = 0.0_dp    , max_Z    = 1.0_dp   , & ! bound for gas; can change for two phase flows
                                               min_FI   = 0.0_dp    , max_FI   = 1.0_dp   , &
                                               min_Da   = 1.0e-3_dp , max_Da   = 1.0e3_dp , &
                                               min_sdr  = 1.0e-2_dp , max_sdr  = 1.0e7_dp , &
                                               min_freq = 0.0_dp    , max_freq = 1.0e3_dp


end module parameters

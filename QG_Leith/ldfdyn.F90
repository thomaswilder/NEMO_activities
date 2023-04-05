MODULE ldfdyn
   !!======================================================================
   !!                       ***  MODULE  ldfdyn  ***
   !! Ocean physics:  lateral viscosity coefficient 
   !!=====================================================================
   !! History :  OPA  ! 1997-07  (G. Madec)  multi dimensional coefficients
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.7  ! 2014-01  (F. Lemarie, G. Madec)  restructuration/simplification of ahm specification,
   !!                 !                                  add velocity dependent coefficient and optional read in file
   !!           4.0.4 ! 2022-10  (T. Wilder)  implementing 2D and QG Leith viscosity
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_dyn_init  : initialization, namelist read, and parameters control
   !!   ldf_dyn       : update lateral eddy viscosity coefficients at each time step 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers   
   USE dom_oce         ! ocean space and time domain 
   USE phycst          ! physical constants
   USE ldfslp          ! lateral diffusion: slopes of mixing orientation
   USE ldfc1d_c2d      ! lateral diffusion: 1D and 2D cases
   USE eosbn2          ! equation of states: QG Leith
   USE zdfmxl          ! mixed layer depth: QG Leith
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module for ehanced bottom friction file
   USE timing          ! Timing
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_dyn_init   ! called by nemogcm.F90
   PUBLIC   ldf_dyn        ! called by step.F90

   !                                    !!* Namelist namdyn_ldf : lateral mixing on momentum *
   LOGICAL , PUBLIC ::   ln_dynldf_OFF   !: No operator (i.e. no explicit diffusion)
   LOGICAL , PUBLIC ::   ln_dynldf_lap   !: laplacian operator
   LOGICAL , PUBLIC ::   ln_dynldf_blp   !: bilaplacian operator
   LOGICAL , PUBLIC ::   ln_dynldf_lev   !: iso-level direction
   LOGICAL , PUBLIC ::   ln_dynldf_hor   !: horizontal (geopotential) direction
!  LOGICAL , PUBLIC ::   ln_dynldf_iso   !: iso-neutral direction                        (see ldfslp)
   INTEGER , PUBLIC ::   nn_ahm_ijk_t    !: choice of time & space variations of the lateral eddy viscosity coef.
   !                                        !  time invariant coefficients:  aht = 1/2  Ud*Ld   (lap case) 
   !                                           !                             bht = 1/12 Ud*Ld^3 (blp case)
   REAL(wp), PUBLIC ::   rn_Uv                 !: lateral viscous velocity  [m/s]
   REAL(wp), PUBLIC ::   rn_Lv                 !: lateral viscous length    [m]
   !                                        ! Smagorinsky viscosity  (nn_ahm_ijk_t = 32) 
   REAL(wp), PUBLIC ::   rn_csmc               !: Smagorinsky constant of proportionality 
   REAL(wp), PUBLIC ::   rn_minfac             !: Multiplicative factor of theorectical minimum Smagorinsky viscosity
   REAL(wp), PUBLIC ::   rn_maxfac             !: Multiplicative factor of theorectical maximum Smagorinsky viscosity
   !                                        ! 2D Leith viscosity  (nn_ahm_ijk_t = 33) 
   REAL(wp), PUBLIC ::   rn_c2dc               !: 2D Leith tuning parameter, typically set to 1
   !                                        ! QG Leith viscosity  (nn_ahm_ijk_t = 34) 
   REAL(wp), PUBLIC ::   rn_cqgc               !: QG Leith tuning parameter, typically set to 1
   !                                        ! iso-neutral laplacian (ln_dynldf_lap=ln_dynldf_iso=T)
   REAL(wp), PUBLIC ::   rn_ahm_b              !: lateral laplacian background eddy viscosity  [m2/s]

   !                                    !!* Parameter to control the type of lateral viscous operator
   INTEGER, PARAMETER, PUBLIC ::   np_ERROR  =-10                       !: error in setting the operator
   INTEGER, PARAMETER, PUBLIC ::   np_no_ldf = 00                       !: without operator (i.e. no lateral viscous trend)
   !                          !!      laplacian     !    bilaplacian    !
   INTEGER, PARAMETER, PUBLIC ::   np_lap    = 10   ,   np_blp    = 20  !: iso-level operator
   INTEGER, PARAMETER, PUBLIC ::   np_lap_i  = 11                       !: iso-neutral or geopotential operator
   !
   INTEGER           , PUBLIC ::   nldf_dyn         !: type of lateral diffusion used defined from ln_dynldf_... (namlist logicals)
   LOGICAL           , PUBLIC ::   l_ldfdyn_time    !: flag for time variation of the lateral eddy viscosity coef.

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ahmt, ahmf   !: eddy viscosity coef. at T- and F-points [m2/s or m4/s]
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dtensq       !: horizontal tension squared         (Smagorinsky only)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dshesq       !: horizontal shearing strain squared (Smagorinsky only)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   esqt, esqf   !: Square of the local gridscale (e1e2/(e1+e2))**2 (Smag) or ( e1e2 ) (2D/QG Leith)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dwzmagsq     !: square of magnitude of gradient of vorticity (2D/QG Leith)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ddivmagsq    !: square of magnitude of gradient of divergence (2D/QG Leith)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hdivnqg      !: Horizontal divergence on t-point (2D/QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zwz          !: Vorticity on f-point (2D/QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zstx, zsty   !: x and y components of stretching at T- points (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zwzdx, zwzdy !: x and y components of horizontal gradients of vertical vorticity at T- points (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zstlimx, zstlimy !: Limit of stretching term at T- points (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zbu          !: Buoyancy at T- point (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rbu, rro2    !: Burger number and square of Rossby number at T- points (QG Leith)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zbudxup, zbudyvp !: gradients of buoyancy - x and y components on U- point and V- points, resp. (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zbudx, zbudy !: x and y components of gradients in buoyancy at T- points (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tmpzstx      !: alternative stretching term computed in QG Leith for diagnostic purposes 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   mld_qg       !: QG Leith mixed layer depth

   REAL(wp) ::   r1_2    = 0.5_wp            ! =1/2
   REAL(wp) ::   r1_4    = 0.25_wp           ! =1/4
   REAL(wp) ::   r1_8    = 0.125_wp          ! =1/8
   REAL(wp) ::   r1_12   = 1._wp / 12._wp    ! =1/12
   REAL(wp) ::   r1_288  = 1._wp / 288._wp   ! =1/( 12^2 * 2 )

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_dyn_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn_init  ***
      !!                   
      !! ** Purpose :   set the horizontal ocean dynamics physics
      !!
      !! ** Method  :   the eddy viscosity coef. specification depends on:
      !!              - the operator:
      !!             ln_dynldf_lap = T     laplacian operator
      !!             ln_dynldf_blp = T   bilaplacian operator
      !!              - the parameter nn_ahm_ijk_t:
      !!    nn_ahm_ijk_t  =  0 => = constant
      !!                  = 10 => = F(z) :     = constant with a reduction of 1/4 with depth 
      !!                  =-20 => = F(i,j)     = shape read in 'eddy_viscosity.nc' file
      !!                  = 20    = F(i,j)     = F(e1,e2) or F(e1^3,e2^3) (lap or bilap case)
      !!                  =-30 => = F(i,j,k)   = shape read in 'eddy_viscosity.nc'  file
      !!                  = 30    = F(i,j,k)   = 2D (case 20) + decrease with depth (case 10)
      !!                  = 31    = F(i,j,k,t) = F(local velocity) (  |u|e  /12   laplacian operator
      !!                                                           or |u|e^3/12 bilaplacian operator )
      !!                  = 32    = F(i,j,k,t) = F(local deformation rate and gridscale) (D and L) (Smagorinsky)  
      !!                                                           (   L^2|D|      laplacian operator
      !!                                                           or  L^4|D|/8  bilaplacian operator )
      !!                  = 33    = F(i,j,k,t) = F( PV gradient, divergence, and gridscale) (laplacian operator) (2D Leith)
      !!                  = 34    = F(i,j,k,t) = F( QG PV gradient, divergence, and gridscale) (laplacian operator) (QG Leith)
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk                     ! dummy loop indices
      INTEGER  ::   ioptio, ierr, inum, ios, inn   ! local integer
      REAL(wp) ::   zah0, zah_max, zUfac           ! local scalar
      CHARACTER(len=5) ::   cl_Units               ! units (m2/s or m4/s)
      !!
      NAMELIST/namdyn_ldf/ ln_dynldf_OFF, ln_dynldf_lap, ln_dynldf_blp,   &   ! type of operator
         &                 ln_dynldf_lev, ln_dynldf_hor, ln_dynldf_iso,   &   ! acting direction of the operator
         &                 nn_ahm_ijk_t , rn_Uv    , rn_Lv,   rn_ahm_b,   &   ! lateral eddy coefficient
         &                 rn_csmc      , rn_minfac    , rn_maxfac,       &   ! Smagorinsky settings
         &                 rn_c2dc,                                       &   ! 2D Leith tuning parameter
         &                 rn_cqgc                                            ! QG Leith tuning parameter
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )
      READ  ( numnam_ref, namdyn_ldf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_ldf in reference namelist' )

      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namdyn_ldf, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn_ldf in configuration namelist' )
      IF(lwm) WRITE ( numond, namdyn_ldf )

      IF(lwp) THEN                      ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_dyn : lateral momentum physics'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_ldf : set lateral mixing parameters'
         !
         WRITE(numout,*) '      type :'
         WRITE(numout,*) '         no explicit diffusion                ln_dynldf_OFF = ', ln_dynldf_OFF
         WRITE(numout,*) '         laplacian operator                   ln_dynldf_lap = ', ln_dynldf_lap
         WRITE(numout,*) '         bilaplacian operator                 ln_dynldf_blp = ', ln_dynldf_blp
         !
         WRITE(numout,*) '      direction of action :'
         WRITE(numout,*) '         iso-level                            ln_dynldf_lev = ', ln_dynldf_lev
         WRITE(numout,*) '         horizontal (geopotential)            ln_dynldf_hor = ', ln_dynldf_hor
         WRITE(numout,*) '         iso-neutral                          ln_dynldf_iso = ', ln_dynldf_iso
         !
         WRITE(numout,*) '      coefficients :'
         WRITE(numout,*) '         type of time-space variation         nn_ahm_ijk_t  = ', nn_ahm_ijk_t
         WRITE(numout,*) '         lateral viscous velocity  (if cst)      rn_Uv      = ', rn_Uv, ' m/s'
         WRITE(numout,*) '         lateral viscous length    (if cst)      rn_Lv      = ', rn_Lv, ' m'
         WRITE(numout,*) '         background viscosity (iso-lap case)     rn_ahm_b   = ', rn_ahm_b, ' m2/s'
         !
         WRITE(numout,*) '      Smagorinsky settings (nn_ahm_ijk_t  = 32) :'
         WRITE(numout,*) '         Smagorinsky coefficient              rn_csmc       = ', rn_csmc
         WRITE(numout,*) '         factor multiplier for eddy visc.'
         WRITE(numout,*) '            lower limit (default 1.0)         rn_minfac     = ', rn_minfac
         WRITE(numout,*) '            upper limit (default 1.0)         rn_maxfac     = ', rn_maxfac
         !
         WRITE(numout,*) '      2D Leith settings (nn_ahm_ijk_t  = 33) :'
         WRITE(numout,*) '         2D Leith coefficient                 rn_c2dc       = ', rn_c2dc
         !
         WRITE(numout,*) '      QG Leith settings (nn_ahm_ijk_t  = 34) :'
         WRITE(numout,*) '         QG Leith coefficient                 rn_cqgc       = ', rn_cqgc
      ENDIF

      !
      !           !==  type of lateral operator used  ==!   (set nldf_dyn)
      !           !=====================================!
      !
      nldf_dyn = np_ERROR
      ioptio = 0
      IF( ln_dynldf_OFF ) THEN   ;   nldf_dyn = np_no_ldf   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_dynldf_lap ) THEN   ;                              ioptio = ioptio + 1   ;   ENDIF
      IF( ln_dynldf_blp ) THEN   ;                              ioptio = ioptio + 1   ;   ENDIF
      IF( ioptio /= 1   )   CALL ctl_stop( 'dyn_ldf_init: use ONE of the 3 operator options (NONE/lap/blp)' )
      !
      IF(.NOT.ln_dynldf_OFF ) THEN     !==  direction ==>> type of operator  ==!
         ioptio = 0
         IF( ln_dynldf_lev )   ioptio = ioptio + 1
         IF( ln_dynldf_hor )   ioptio = ioptio + 1
         IF( ln_dynldf_iso )   ioptio = ioptio + 1
         IF( ioptio /= 1   )   CALL ctl_stop( 'dyn_ldf_init: use ONE of the 3 direction options (level/hor/iso)' )
         !
         !                             ! Set nldf_dyn, the type of lateral diffusion, from ln_dynldf_... logicals
         ierr = 0
         IF( ln_dynldf_lap ) THEN         ! laplacian operator
            IF( ln_zco ) THEN                   ! z-coordinate
               IF ( ln_dynldf_lev )   nldf_dyn = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_hor )   nldf_dyn = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_iso )   nldf_dyn = np_lap_i   ! iso-neutral            (   rotation)
            ENDIF
            IF( ln_zps ) THEN                   ! z-coordinate with partial step
               IF ( ln_dynldf_lev )   nldf_dyn = np_lap     ! iso-level              (no rotation)
               IF ( ln_dynldf_hor )   nldf_dyn = np_lap     ! iso-level              (no rotation)
               IF ( ln_dynldf_iso )   nldf_dyn = np_lap_i   ! iso-neutral            (   rotation)
            ENDIF
            IF( ln_sco ) THEN                   ! s-coordinate
               IF ( ln_dynldf_lev )   nldf_dyn = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_hor )   nldf_dyn = np_lap_i   ! horizontal             (   rotation)
               IF ( ln_dynldf_iso )   nldf_dyn = np_lap_i   ! iso-neutral            (   rotation)
            ENDIF
         ENDIF
         !
         IF( ln_dynldf_blp ) THEN         ! bilaplacian operator
            IF( ln_zco ) THEN                   ! z-coordinate
               IF( ln_dynldf_lev )   nldf_dyn = np_blp   ! iso-level = horizontal (no rotation)
               IF( ln_dynldf_hor )   nldf_dyn = np_blp   ! iso-level = horizontal (no rotation)
               IF( ln_dynldf_iso )   ierr = 2            ! iso-neutral            (   rotation)
            ENDIF
            IF( ln_zps ) THEN                   ! z-coordinate with partial step
               IF( ln_dynldf_lev )   nldf_dyn = np_blp   ! iso-level              (no rotation)
               IF( ln_dynldf_hor )   nldf_dyn = np_blp   ! iso-level              (no rotation)
               IF( ln_dynldf_iso )   ierr = 2            ! iso-neutral            (   rotation)
            ENDIF
            IF( ln_sco ) THEN                   ! s-coordinate
               IF( ln_dynldf_lev )   nldf_dyn = np_blp   ! iso-level              (no rotation)
               IF( ln_dynldf_hor )   ierr = 2            ! horizontal             (   rotation)
               IF( ln_dynldf_iso )   ierr = 2            ! iso-neutral            (   rotation)
            ENDIF
         ENDIF
         !
         IF( ierr == 2 )   CALL ctl_stop( 'rotated bi-laplacian operator does not exist' )
         !
         IF( nldf_dyn == np_lap_i )   l_ldfslp = .TRUE.  ! rotation require the computation of the slopes
         !
      ENDIF
      !
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE( nldf_dyn )
         CASE( np_no_ldf )   ;   WRITE(numout,*) '   ==>>>   NO lateral viscosity'
         CASE( np_lap    )   ;   WRITE(numout,*) '   ==>>>   iso-level laplacian operator'
         CASE( np_lap_i  )   ;   WRITE(numout,*) '   ==>>>   rotated laplacian operator with iso-level background'
         CASE( np_blp    )   ;   WRITE(numout,*) '   ==>>>   iso-level bi-laplacian operator'
         END SELECT
         WRITE(numout,*)
      ENDIF
      
      !
      !           !==  Space/time variation of eddy coefficients  ==!
      !           !=================================================!
      !
      l_ldfdyn_time = .FALSE.                ! no time variation except in case defined below
      !
      IF( ln_dynldf_OFF ) THEN
         IF(lwp) WRITE(numout,*) '   ==>>>   No viscous operator selected. ahmt and ahmf are not allocated'
         RETURN
         !
      ELSE                                   !==  a lateral diffusion operator is used  ==!
         !
         !                                         ! allocate the ahm arrays
         ALLOCATE( ahmt(jpi,jpj,jpk) , ahmf(jpi,jpj,jpk) , STAT=ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ldf_dyn_init: failed to allocate arrays')
         !
         ahmt(:,:,:) = 0._wp                       ! init to 0 needed 
         ahmf(:,:,:) = 0._wp
         !
         !                                         ! value of lap/blp eddy mixing coef.
         IF(     ln_dynldf_lap ) THEN   ;   zUfac = r1_2 *rn_Uv   ;   inn = 1   ;   cl_Units = ' m2/s'   !   laplacian
         ELSEIF( ln_dynldf_blp ) THEN   ;   zUfac = r1_12*rn_Uv   ;   inn = 3   ;   cl_Units = ' m4/s'   ! bilaplacian
         ENDIF
         zah0    = zUfac *    rn_Lv**inn              ! mixing coefficient
         zah_max = zUfac * (ra*rad)**inn              ! maximum reachable coefficient (value at the Equator)
         !
         SELECT CASE(  nn_ahm_ijk_t  )             !* Specification of space-time variations of ahmt, ahmf
         !
         CASE(   0  )      !==  constant  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity. = constant = ', zah0, cl_Units
            ahmt(:,:,1:jpkm1) = zah0
            ahmf(:,:,1:jpkm1) = zah0
            !
         CASE(  10  )      !==  fixed profile  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( depth )'
            IF(lwp) WRITE(numout,*) '           surface viscous coef. = constant = ', zah0, cl_Units
            ahmt(:,:,1) = zah0                        ! constant surface value
            ahmf(:,:,1) = zah0
            CALL ldf_c1d( 'DYN', ahmt(:,:,1), ahmf(:,:,1), ahmt, ahmf )
            !
         CASE ( -20 )      !== fixed horizontal shape read in file  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F(i,j) read in eddy_viscosity.nc file'
            CALL iom_open( 'eddy_viscosity_2D.nc', inum )
            CALL iom_get ( inum, jpdom_data, 'ahmt_2d', ahmt(:,:,1) )
            CALL iom_get ( inum, jpdom_data, 'ahmf_2d', ahmf(:,:,1) )
            CALL iom_close( inum )
            DO jk = 2, jpkm1
               ahmt(:,:,jk) = ahmt(:,:,1)
               ahmf(:,:,jk) = ahmf(:,:,1)
            END DO
            !
         CASE(  20  )      !== fixed horizontal shape  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( e1, e2 ) or F( e1^3, e2^3 ) (lap. or blp. case)'
            IF(lwp) WRITE(numout,*) '           using a fixed viscous velocity = ', rn_Uv  ,' m/s   and   Lv = Max(e1,e2)'
            IF(lwp) WRITE(numout,*) '           maximum reachable coefficient (at the Equator) = ', zah_max, cl_Units, '  for e1=1°)'
            CALL ldf_c2d( 'DYN', zUfac      , inn        , ahmt, ahmf )         ! surface value proportional to scale factor^inn
            !
         CASE( -30  )      !== fixed 3D shape read in file  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F(i,j,k) read in eddy_viscosity_3D.nc file'
            CALL iom_open( 'eddy_viscosity_3D.nc', inum )
            CALL iom_get ( inum, jpdom_data, 'ahmt_3d', ahmt )
            CALL iom_get ( inum, jpdom_data, 'ahmf_3d', ahmf )
            CALL iom_close( inum )
            !
         CASE(  30  )       !==  fixed 3D shape  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( latitude, longitude, depth )'
            IF(lwp) WRITE(numout,*) '           using a fixed viscous velocity = ', rn_Uv  ,' m/s   and   Ld = Max(e1,e2)'
            IF(lwp) WRITE(numout,*) '           maximum reachable coefficient (at the Equator) = ', zah_max, cl_Units, '  for e1=1°)'
            CALL ldf_c2d( 'DYN', zUfac      , inn        , ahmt, ahmf )         ! surface value proportional to scale factor^inn
            CALL ldf_c1d( 'DYN', ahmt(:,:,1), ahmf(:,:,1), ahmt, ahmf )  ! reduction with depth
            !
         CASE(  31  )       !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( latitude, longitude, depth , time )'
            IF(lwp) WRITE(numout,*) '           proportional to the local velocity : 1/2 |u|e (lap) or 1/12 |u|e^3 (blp)'
            !
            l_ldfdyn_time = .TRUE.     ! will be calculated by call to ldf_dyn routine in step.F90
            !
         CASE(  32  )       !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( latitude, longitude, depth , time )'
            IF(lwp) WRITE(numout,*) '           proportional to the local deformation rate and gridscale (Smagorinsky)'
            !
            l_ldfdyn_time = .TRUE.     ! will be calculated by call to ldf_dyn routine in step.F90
            !
            !                          ! allocate arrays used in ldf_dyn. 
            ALLOCATE( dtensq(jpi,jpj,jpk) , dshesq(jpi,jpj,jpk) , esqt(jpi,jpj) , esqf(jpi,jpj) , STAT=ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ldf_dyn_init: failed to allocate Smagorinsky arrays')
            !
            DO jj = 1, jpj             ! Set local gridscale values
               DO ji = 1, jpi
                  esqt(ji,jj) = ( 2._wp * e1e2t(ji,jj) / ( e1t(ji,jj) + e2t(ji,jj) ) )**2 
                  esqf(ji,jj) = ( 2._wp * e1e2f(ji,jj) / ( e1f(ji,jj) + e2f(ji,jj) ) )**2 
               END DO
            END DO
            !
         CASE(  33  )       !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( latitude, longitude, depth , time )'
            IF(lwp) WRITE(numout,*) '           proportional to the PV gradient, divergence, and gridscale (2D Leith)'
            !
            l_ldfdyn_time = .TRUE.     ! will be calculated by call to ldf_dyn routine in step.F90
            !
            !                          ! allocate arrays used in ldf_dyn. 
            ALLOCATE( esqt(jpi,jpj) , esqf(jpi,jpj) , dwzmagsq(jpi,jpj,jpk) , ddivmagsq(jpi,jpj,jpk) ,      &
               &  zwz(jpi,jpj,jpk) , hdivnqg(jpi,jpj,jpk) , STAT=ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ldf_dyn_init: failed to allocate 2D Leith arrays')
            !
            DO jj = 1, jpj             ! Set local gridscale values
               DO ji = 1, jpi
                  esqt(ji,jj) = e1e2t(ji,jj)
                  esqf(ji,jj) = e1e2f(ji,jj)
               END DO
            END DO
            !
            !== initialise key variables with zeros ==!
            hdivnqg(:,:,:) = 0._wp
            zwz(:,:,:) = 0._wp
            dwzmagsq(:,:,:) = 0._wp
            ddivmagsq(:,:,:) = 0._wp
            !
         CASE(  34  )       !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( latitude, longitude, depth , time )'
            IF(lwp) WRITE(numout,*) '           proportional to the PV gradient, divergence, and gridscale (QG Leith)'
            !
            l_ldfdyn_time = .TRUE.     ! will be calculated by call to ldf_dyn routine in step.F90
            !
            !                          ! allocate arrays used in ldf_dyn. 
            ALLOCATE( esqt(jpi,jpj) , esqf(jpi,jpj) , dwzmagsq(jpi,jpj,jpk) , ddivmagsq(jpi,jpj,jpk) ,      &
               &  zbu(jpi,jpj,jpk) , zbudxup(jpi,jpj,jpk) , zbudyvp(jpi,jpj,jpk) , zwz(jpi,jpj,jpk) ,       &
               &  zstlimx(jpi,jpj,jpk) , zstlimy(jpi,jpj,jpk) , zstx(jpi,jpj,jpk) , zsty(jpi,jpj,jpk) ,     &
               &  rbu(jpi,jpj,jpk), rro2(jpi,jpj,jpk) , zwzdx(jpi,jpj,jpk) , zwzdy(jpi,jpj,jpk) ,           &
               &  zbudx(jpi,jpj,jpk) , zbudy(jpi,jpj,jpk) , hdivnqg(jpi,jpj,jpk) ,                          &
               &  tmpzstx(jpi,jpj,jpk) , mld_qg(jpi,jpj) , STAT=ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ldf_dyn_init: failed to allocate QG Leith arrays')
            !
            DO jj = 1, jpj             ! Set local gridscale values
               DO ji = 1, jpi
                  esqt(ji,jj) = e1e2t(ji,jj)
                  esqf(ji,jj) = e1e2f(ji,jj)
               END DO
            END DO
            !
            !== initialise key variables with zeros ==!
            dwzmagsq(:,:,:) = 0._wp
            ddivmagsq(:,:,:) = 0._wp
            zbu(:,:,:) = 0._wp
            zbudxup(:,:,:) = 0._wp
            zbudyvp(:,:,:) = 0._wp
            zbudx(:,:,:) = 0._wp
            zbudy(:,:,:) = 0._wp
            zstx(:,:,:) = 0._wp
            zsty(:,:,:) = 0._wp
            zstlimx(:,:,:) = 0._wp
            zstlimy(:,:,:) = 0._wp
            zwz(:,:,:) = 0._wp
            zwzdx(:,:,:) = 0._wp
            zwzdy(:,:,:) = 0._wp
            rbu(:,:,:) = 0._wp
            rro2(:,:,:) = 0._wp
            hdivnqg(:,:,:) = 0._wp
            tmpzstx(:,:,:) = 0._wp
            mld_qg(:,:) = 0._wp
            !
         CASE DEFAULT
            CALL ctl_stop('ldf_dyn_init: wrong choice for nn_ahm_ijk_t, the type of space-time variation of ahm')
         END SELECT
         !
         IF( .NOT.l_ldfdyn_time ) THEN             !* No time variation 
            IF(     ln_dynldf_lap ) THEN                 !   laplacian operator (mask only)
               ahmt(:,:,1:jpkm1) =       ahmt(:,:,1:jpkm1)   * tmask(:,:,1:jpkm1)
               ahmf(:,:,1:jpkm1) =       ahmf(:,:,1:jpkm1)   * fmask(:,:,1:jpkm1)
            ELSEIF( ln_dynldf_blp ) THEN                 ! bilaplacian operator (square root + mask)
               ahmt(:,:,1:jpkm1) = SQRT( ahmt(:,:,1:jpkm1) ) * tmask(:,:,1:jpkm1)
               ahmf(:,:,1:jpkm1) = SQRT( ahmf(:,:,1:jpkm1) ) * fmask(:,:,1:jpkm1)
            ENDIF
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE ldf_dyn_init


   SUBROUTINE ldf_dyn( kt, prd, pn2 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn  ***
      !! 
      !! ** Purpose :   update at kt the momentum lateral mixing coeff. (ahmt and ahmf)
      !!
      !! ** Method  :   time varying eddy viscosity coefficients:
      !!
      !!    nn_ahm_ijk_t = 31    ahmt, ahmf = F(i,j,k,t) = F(local velocity) 
      !!                         ( |u|e /12  or  |u|e^3/12 for laplacian or bilaplacian operator )
      !!
      !!    nn_ahm_ijk_t = 32    ahmt, ahmf = F(i,j,k,t) = F(local deformation rate and gridscale) (D and L) (Smagorinsky)  
      !!                         ( L^2|D|    or  L^4|D|/8  for laplacian or bilaplacian operator )
      !!
      !!    nn_ahm_ijk_t = 33    ahmt, ahmf = F(i,j,k,t) = F( grad( PV and divergence), and gridscale) (laplacian operator) (2D Leith)
      !!
      !!    nn_ahm_ijk_t = 34    ahmt, ahmf = F(i,j,k,t) = F( grad( QG PV and divergence), and gridscale) (laplacian operator) (QG Leith)
      !!
      !! ** note    :    in BLP cases the sqrt of the eddy coef is returned, since bilaplacian is en re-entrant laplacian
      !! ** action  :    ahmt, ahmf   updated at each time step
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step index
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   prd   ! in situ density
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   pn2   ! Brunt-Vaisala frequency (locally ref.)
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zu2pv2_ij_p1, zu2pv2_ij, zu2pv2_ij_m1, zemax   ! local scalar (option 31)
      REAL(wp) ::   zcmsmag, zstabf_lo, zstabf_up, zdelta, zdb     ! local scalar (option 32)
      REAL(wp) ::   zcm2dl, zsq2d , zztmpx, zztmpy                 ! local scalar (option 33)
      REAL(wp) ::   zcmqgl, zbuup, zbulw, zusq, znsq               ! local scalar (option 34)
      REAL(wp) ::   zker1, zker2, zqglep1, zqglep2, zsqqg          ! more local scalar (option 34)
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ldf_dyn')
      !
      SELECT CASE(  nn_ahm_ijk_t  )       !== Eddy vicosity coefficients ==!
      !
      CASE(  31  )       !==  time varying 3D field  ==!   = F( local velocity )
         !
         IF( ln_dynldf_lap   ) THEN        ! laplacian operator : |u| e /12 = |u/144| e
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zu2pv2_ij_m1 = ub(ji-1,jj  ,jk) * ub(ji-1,jj  ,jk) + vb(ji  ,jj-1,jk) * vb(ji  ,jj-1,jk)
                     zemax = MAX( e1t(ji,jj) , e2t(ji,jj) )
                     ahmt(ji,jj,jk) = SQRT( (zu2pv2_ij + zu2pv2_ij_m1) * r1_288 ) * zemax * tmask(ji,jj,jk)      ! 288= 12*12 * 2
                  END DO
               END DO
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     zu2pv2_ij_p1 = ub(ji  ,jj+1,jk) * ub(ji  ,jj+1,jk) + vb(ji+1,jj  ,jk) * vb(ji+1,jj  ,jk)
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zemax = MAX( e1f(ji,jj) , e2f(ji,jj) )
                     ahmf(ji,jj,jk) = SQRT( (zu2pv2_ij + zu2pv2_ij_p1) * r1_288 ) * zemax * fmask(ji,jj,jk)      ! 288= 12*12 * 2
                  END DO
               END DO
            END DO
         ELSEIF( ln_dynldf_blp ) THEN      ! bilaplacian operator : sqrt( |u| e^3 /12 ) = sqrt( |u/144| e ) * e
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zu2pv2_ij_m1 = ub(ji-1,jj  ,jk) * ub(ji-1,jj  ,jk) + vb(ji  ,jj-1,jk) * vb(ji  ,jj-1,jk)
                     zemax = MAX( e1t(ji,jj) , e2t(ji,jj) )
                     ahmt(ji,jj,jk) = SQRT(  SQRT( (zu2pv2_ij + zu2pv2_ij_m1) * r1_288 ) * zemax  ) * zemax * tmask(ji,jj,jk)
                  END DO
               END DO
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     zu2pv2_ij_p1 = ub(ji  ,jj+1,jk) * ub(ji  ,jj+1,jk) + vb(ji+1,jj  ,jk) * vb(ji+1,jj  ,jk)
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zemax = MAX( e1f(ji,jj) , e2f(ji,jj) )
                     ahmf(ji,jj,jk) = SQRT(  SQRT( (zu2pv2_ij + zu2pv2_ij_p1) * r1_288 ) * zemax  ) * zemax * fmask(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF
         !
         CALL lbc_lnk_multi( 'ldfdyn', ahmt, 'T', 1.,  ahmf, 'F', 1. )
         !
         !
      CASE(  32  )       !==  time varying 3D field  ==!   = F( local deformation rate and gridscale ) (Smagorinsky)
         !
         IF( ln_dynldf_lap .OR. ln_dynldf_blp  ) THEN        ! laplacian operator : (C_smag/pi)^2 L^2 |D|
            !
            zcmsmag   = (rn_csmc/rpi)**2                                            ! (C_smag/pi)^2
            zstabf_lo = rn_minfac * rn_minfac / ( 2._wp * 12._wp * 12._wp * zcmsmag ) ! lower limit stability factor scaling
            zstabf_up = rn_maxfac / ( 4._wp * zcmsmag * 2._wp * rdt )               ! upper limit stability factor scaling
            IF( ln_dynldf_blp ) zstabf_lo = ( 16._wp / 9._wp ) * zstabf_lo          ! provide |U|L^3/12 lower limit instead 
            !                                                                       ! of |U|L^3/16 in blp case
            DO jk = 1, jpkm1
               !
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     zdb =   ( ub(ji,jj,jk) * r1_e2u(ji,jj) - ub(ji-1,jj,jk) * r1_e2u(ji-1,jj) ) * r1_e1t(ji,jj) * e2t(ji,jj)  &
                        &  - ( vb(ji,jj,jk) * r1_e1v(ji,jj) - vb(ji,jj-1,jk) * r1_e1v(ji,jj-1) ) * r1_e2t(ji,jj) * e1t(ji,jj)
                     dtensq(ji,jj,jk) = zdb * zdb * tmask(ji,jj,jk)
                  END DO
               END DO
               !
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     zdb =   ( ub(ji,jj+1,jk) * r1_e1u(ji,jj+1) - ub(ji,jj,jk) * r1_e1u(ji,jj) ) * r1_e2f(ji,jj) * e1f(ji,jj)  &
                        &  + ( vb(ji+1,jj,jk) * r1_e2v(ji+1,jj) - vb(ji,jj,jk) * r1_e2v(ji,jj) ) * r1_e1f(ji,jj) * e2f(ji,jj)
                     dshesq(ji,jj,jk) = zdb * zdb * fmask(ji,jj,jk)
                     END DO
               END DO
               !
            END DO
            !
            CALL lbc_lnk_multi( 'ldfdyn', dtensq, 'T', 1. )  ! lbc_lnk on dshesq not needed
            !
            DO jk = 1, jpkm1
              !
               DO jj = 2, jpjm1                                ! T-point value
                  DO ji = fs_2, fs_jpim1
                     !
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zu2pv2_ij_m1 = ub(ji-1,jj  ,jk) * ub(ji-1,jj  ,jk) + vb(ji  ,jj-1,jk) * vb(ji  ,jj-1,jk)
                     !
                     zdelta         = zcmsmag * esqt(ji,jj)                                        ! L^2 * (C_smag/pi)^2
                     ahmt(ji,jj,jk) = zdelta * SQRT(          dtensq(ji  ,jj,jk) +                         &
                        &                            r1_4 * ( dshesq(ji  ,jj,jk) + dshesq(ji  ,jj-1,jk) +  &
                        &                                     dshesq(ji-1,jj,jk) + dshesq(ji-1,jj-1,jk) ) )
                     ahmt(ji,jj,jk) = MAX( ahmt(ji,jj,jk), SQRT( (zu2pv2_ij + zu2pv2_ij_m1) * zdelta * zstabf_lo ) ) ! Impose lower limit == minfac  * |U|L/2
                     ahmt(ji,jj,jk) = MIN( ahmt(ji,jj,jk),                                    zdelta * zstabf_up )   ! Impose upper limit == maxfac  * L^2/(4*2dt)
                     !
                  END DO
               END DO
               !
               DO jj = 1, jpjm1                                ! F-point value
                  DO ji = 1, fs_jpim1
                     !
                     zu2pv2_ij_p1 = ub(ji  ,jj+1,jk) * ub(ji  ,jj+1,jk) + vb(ji+1,jj  ,jk) * vb(ji+1,jj  ,jk)
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     !
                     zdelta         = zcmsmag * esqf(ji,jj)                                        ! L^2 * (C_smag/pi)^2
                     ahmf(ji,jj,jk) = zdelta * SQRT(          dshesq(ji  ,jj,jk) +                         &
                        &                            r1_4 * ( dtensq(ji  ,jj,jk) + dtensq(ji  ,jj+1,jk) +  &
                        &                                     dtensq(ji+1,jj,jk) + dtensq(ji+1,jj+1,jk) ) )
                     ahmf(ji,jj,jk) = MAX( ahmf(ji,jj,jk), SQRT( (zu2pv2_ij + zu2pv2_ij_p1) * zdelta * zstabf_lo ) ) ! Impose lower limit == minfac  * |U|L/2
                     ahmf(ji,jj,jk) = MIN( ahmf(ji,jj,jk),                                    zdelta * zstabf_up )   ! Impose upper limit == maxfac  * L^2/(4*2dt)
                     !
                  END DO
               END DO
               !
            END DO
            !
         ENDIF
         !
         IF( ln_dynldf_blp ) THEN      ! bilaplacian operator : sqrt( (C_smag/pi)^2 L^4 |D|/8)
            !                          !                      = sqrt( A_lap_smag L^2/8 )
            !                          ! stability limits already applied to laplacian values
            !                          ! effective default limits are 1/12 |U|L^3 < B_hm < 1//(32*2dt) L^4
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ahmt(ji,jj,jk) = SQRT( r1_8 * esqt(ji,jj) * ahmt(ji,jj,jk) )
                  END DO
               END DO
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     ahmf(ji,jj,jk) = SQRT( r1_8 * esqf(ji,jj) * ahmf(ji,jj,jk) )
                  END DO
               END DO
            END DO
            !
         ENDIF
         !
         CALL lbc_lnk_multi( 'ldfdyn', ahmt, 'T', 1. , ahmf, 'F', 1. )
         !
         !
      CASE(  33  )       !==  time varying 3D field  ==!   = F( PV gradient, divergence, and gridscale ) (2D Leith)
         !
         IF( ln_dynldf_lap ) THEN        ! laplacian operator
            !
            ! allocate local variables !
            zcm2dl = (rn_c2dc/rpi)**6        			! (C_2d/pi)^6
            !
            !== calculate vertical vorticity (f + zeta) on f-point ==!
            DO jk = 1, jpkm1                                 ! Horizontal slab
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     zwz(ji,jj,jk) = ff_f(ji,jj) + ( e2v(ji+1,jj  ) * vn(ji+1,jj  ,jk) - e2v(ji,jj) * vn(ji,jj,jk)            &
                        &          - e1u(ji  ,jj+1) * un(ji  ,jj+1,jk) + e1u(ji,jj) * un(ji,jj,jk)  ) * r1_e1e2f(ji,jj)
                  END DO
               END DO
            END DO
            !
            !== calculate gradients of vorticity, then square of magnitude (t-point) ==!
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     zztmpx = r1_2 * ( r1_e1f(ji,jj-1) * ( zwz(ji,jj-1,jk) - zwz(ji-1,jj-1,jk) )            &
                        &                 + r1_e1f(ji,jj) * ( zwz(ji,jj,jk) - zwz(ji-1,jj,jk) ) )
                     zztmpy = r1_2 * ( r1_e2f(ji-1,jj) * ( zwz(ji-1,jj,jk) - zwz(ji-1,jj-1,jk) )            &
                        &                 + r1_e2f(ji,jj) * ( zwz(ji,jj,jk) - zwz(ji,jj-1,jk) ) )
                     dwzmagsq(ji,jj,jk) = ( zztmpx * zztmpx + zztmpy * zztmpy ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
            CALL lbc_lnk_multi( 'ldfdyn', dwzmagsq, 'T', 1. )
!            !
!            CALL div_hor( kt )
            !
            DO jk = 1, jpkm1                                      !==  Horizontal divergence  ==!
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     hdivnqg(ji,jj,jk) = (  e2u(ji  ,jj) * e3u_n(ji  ,jj,jk) * un(ji  ,jj,jk)    &
                        &               - e2u(ji-1,jj) * e3u_n(ji-1,jj,jk) * un(ji-1,jj,jk)      &
                        &               + e1v(ji,jj  ) * e3v_n(ji,jj  ,jk) * vn(ji,jj  ,jk)      &
                        &               - e1v(ji,jj-1) * e3v_n(ji,jj-1,jk) * vn(ji,jj-1,jk)  )   &
                        &            * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
                  END DO  
               END DO  
            END DO
            !
            !== calculate gradients of divergence, then square of magnitude (f-point) ==!
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     zztmpx = r1_2 * ( r1_e1t(ji,jj+1) * ( hdivnqg(ji+1,jj+1,jk) - hdivnqg(ji,jj+1,jk) )               &
                        &  + r1_e1t(ji,jj) * ( hdivnqg(ji+1,jj,jk) - hdivnqg(ji,jj,jk) ) )
                     zztmpy = r1_2 * ( r1_e2t(ji+1,jj) * ( hdivnqg(ji+1,jj+1,jk) - hdivnqg(ji+1,jj,jk) )               &
                        &  + r1_e2t(ji,jj) * ( hdivnqg(ji,jj+1,jk) - hdivnqg(ji,jj,jk) ) )
                     ddivmagsq(ji,jj,jk) = ( zztmpx * zztmpx + zztmpy * zztmpy ) * fmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
            DO jk = 1, jpkm1	         !== 2D Leith viscosity coefficient on T-point ==!
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1 ! vector opt.
                     zsq2d = r1_4 * ( ddivmagsq(ji,jj,jk) + ddivmagsq(ji-1,jj,jk) + ddivmagsq(ji,jj-1,jk) +     &
                        &  ddivmagsq(ji-1,jj-1,jk) ) + dwzmagsq(ji,jj,jk)
                     ahmt(ji,jj,jk) = SQRT( zcm2dl * esqt(ji,jj)**3 * zsq2d )
                  END DO
               END DO
            END DO
            !
            DO jk = 1, jpkm1            !== 2D Leith viscosity coefficient on F-point ==!
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1 ! vector opt.
                     zsq2d = r1_4 * ( dwzmagsq(ji,jj,jk) + dwzmagsq(ji+1,jj,jk) + dwzmagsq(ji,jj+1,jk) +     &
                        &  dwzmagsq(ji+1,jj+1,jk) ) + ddivmagsq(ji,jj,jk)
                     ahmf(ji,jj,jk) = SQRT( zcm2dl * esqf(ji,jj)**3 * zsq2d )
                  END DO
               END DO
            END DO
            !
         ENDIF
         !
         CALL lbc_lnk_multi( 'ldfdyn', ahmt, 'T', 1.,  ahmf, 'F', 1. )
         !
         !
      CASE(  34  )       !==  time varying 3D field  ==!   = F(QG PV gradient, divergence, and gridscale) (QG Leith)
         !
         IF( ln_dynldf_lap ) THEN        ! laplacian operator
            !
            ! allocate local variables !
            zcmqgl = (rn_cqgc/rpi)**6         			! (C_qg/pi)^6
            zqglep1 = 1.e-12_wp
            zqglep2 = 1.e-24_wp
            !
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     IF( jk <= nmlnqg(ji,jj) ) mld_qg(ji,jj) = mld_qg(ji,jj) + e3w_n(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
            !== Assess the depth of the mixed layer. For jpk within mixed layer, choose 2D Leith routine, if below, choose QG Routine ==!
            !== Within mixed layer and at ocean bottom => 2D Leith scheme ==!
            !== begin calculation of stretching term d/dz[(f/(N**2))*b] ==!
            ! find buoyancy and interpolate onto w-grid !
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     IF( jk < 2 ) THEN
                        !== buoyancy at surface ==!
                        zbu(ji,jj,jk) = - grav * prd(ji,jj,jk)
                     ELSE
                        !== buoyancy below surface ==!
                        zbuup = - grav * prd(ji,jj,jk-1)
                        zbulw = - grav * prd(ji,jj,jk  )
                        zbu(ji,jj,jk) = 0.5_wp * ( zbuup + zbulw )
                     ENDIF
                  END DO
               END DO
            END DO
            !
            CALL lbc_lnk_multi( 'ldfdyn', zbu, 'T', 1. )
            !
            !== Calculate horizontal gradients of buoyancy and put on w-grid ==!
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     !== are we below the mixed layer? ==!
                     !!IF( jk > nmln(ji,jj) ) THEN
                     !== gradients of buoyancy on U and V grid at w point of cell ==!
                     zbudxup(ji,jj,jk) = r1_e1t(ji,jj) * ( zbu(ji+1,jj,jk) - zbu(ji,jj,jk) )
                     zbudyvp(ji,jj,jk) = r1_e2t(ji,jj) * ( zbu(ji,jj+1,jk) - zbu(ji,jj,jk) )
                     !!ENDIF
                  END DO
               END DO
            END DO
            !
            CALL lbc_lnk_multi( 'ldfdyn', zbudxup, 'U', 1. )
            CALL lbc_lnk_multi( 'ldfdyn', zbudyvp, 'V', 1. )
            !
            !== gradients of buoyancy on T- points ==!
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     !== are we below the mixed layer? ==!
                     !!IF( jk > nmln(ji,jj) ) THEN
                     zbudx(ji,jj,jk) = r1_2 * ( zbudxup(ji-1,jj,jk) + zbudxup(ji,jj,jk) )
                     zbudy(ji,jj,jk) = r1_2 * ( zbudyvp(ji,jj-1,jk) + zbudyvp(ji,jj,jk) )
                     !!ENDIF
                  END DO
               END DO
            END DO
            !
            CALL lbc_lnk_multi( 'ldfdyn', zbudx, 'T', 1. )
            CALL lbc_lnk_multi( 'ldfdyn', zbudy, 'T', 1. )
            !
            !== take vertical gradient and find stretching d/dz[(f * grad(b))/N^2] (t-point) ==!
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     !== are we below the mixed layer and above the sea floor? ==!
                     IF( jk > nmlnqg(ji,jj) .AND. jk < jpkm1 ) THEN
                        !== vertical gradient of x component ==!
                        zker1 = ( ff_t(ji,jj) * zbudx(ji,jj,jk  ) ) / MAX( pn2(ji,jj,jk  ), zqglep1 ) 
                        zker2 = ( ff_t(ji,jj) * zbudx(ji,jj,jk+1) ) / MAX( pn2(ji,jj,jk+1), zqglep1 ) 
                        zstx(ji,jj,jk) = ( ( zker1 - zker2 )  / e3w_n(ji,jj,jk) ) * tmask(ji,jj,jk)
                        !== vertical gradient of y component ==!
                        zker1 = ( ff_t(ji,jj) * zbudy(ji,jj,jk  ) ) / MAX( pn2(ji,jj,jk  ), zqglep1 )
                        zker2 = ( ff_t(ji,jj) * zbudy(ji,jj,jk+1) ) / MAX( pn2(ji,jj,jk+1), zqglep1 )
                        zsty(ji,jj,jk) = ( ( zker1 - zker2 )  / e3w_n(ji,jj,jk) ) * tmask(ji,jj,jk)
                     ENDIF
                  END DO
               END DO
            END DO
            !
            CALL lbc_lnk_multi( 'ldfdyn', zstx, 'T', 1. )
            CALL lbc_lnk_multi( 'ldfdyn', zsty, 'T', 1. )
            !
            !== calculate vertical vorticity (f+zeta) on f-point ==!
            !== calculate meridional gradient of coriolis parameter f and add to zwzdy ==!
            DO jk = 1, jpkm1                                 ! Horizontal slab
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     !== are we in the mixed layer or at the ocean bottom? ==!
!                     IF( jk <= nmln(ji,jj) .OR. jk == jpkm1 ) THEN
                     zwz(ji,jj,jk) = ff_f(ji,jj) + (  e2v(ji+1,jj  ) * vn(ji+1,jj  ,jk) - e2v(ji,jj) * vn(ji,jj,jk)            &
                        &          - e1u(ji  ,jj+1) * un(ji  ,jj+1,jk) + e1u(ji,jj) * un(ji,jj,jk)  ) * r1_e1e2f(ji,jj)
!                     ELSE
!                        !== are we below the mixed layer? ff_f(ji,jj) + meridional gradient? ==!
!                        zwz(ji,jj,jk) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ,jk) - e2v(ji,jj) * vn(ji,jj,jk)            &
!                           &          - e1u(ji  ,jj+1) * un(ji  ,jj+1,jk) + e1u(ji,jj) * un(ji,jj,jk)  ) * r1_e1e2f(ji,jj)
!                     ENDIF
                  END DO
               END DO
            END DO
            !
            !== calculate horizontal gradient of vertical vorticity on t-point ==!
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     zwzdx(ji,jj,jk) = r1_2 * ( r1_e1f(ji,jj-1) * ( zwz(ji,jj-1,jk) - zwz(ji-1,jj-1,jk) )            &
                        &                 + r1_e1f(ji,jj) * ( zwz(ji,jj,jk) - zwz(ji-1,jj,jk) ) )
                     zwzdy(ji,jj,jk) = r1_2 * ( r1_e2f(ji-1,jj) * ( zwz(ji-1,jj,jk) - zwz(ji-1,jj-1,jk) )            &
                        &                 + r1_e2f(ji,jj) * ( zwz(ji,jj,jk) - zwz(ji,jj-1,jk) ) )
                  END DO
               END DO
            END DO
            !
            CALL lbc_lnk_multi( 'ldfdyn', zwzdx, 'T', 1. )
            CALL lbc_lnk_multi( 'ldfdyn', zwzdy, 'T', 1. )
            !
            !== calculate the Burger number and square of Rossby number on t-point ==!
            !== calculate over entire domain for diagnostics ==!
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     !== grid scale velocity squared ==!
                     zztmpx = 0.5 * ( un(ji-1,jj  ,jk) + un(ji,jj,jk) )
                     zztmpy = 0.5 * ( vn(ji  ,jj-1,jk) + vn(ji,jj,jk) )
                     zusq =  zztmpx**2 + zztmpy**2
                     !== square of Rossby number U^2/(f^2 * A) ==!
                     rro2(ji,jj,jk) = ( zusq / ( MAX( ff_t(ji,jj)**2, zqglep2 ) * esqt(ji,jj) ) ) * fmask(ji,jj,jk)
                     !== averaging square of buoyancy frequency onto t-grid ==!
                     IF( jk < jpkm1 ) THEN
                        !== accounting for negative N^2 ==!
                        znsq = r1_2 * ( pn2(ji,jj,jk) + pn2(ji,jj,jk+1) )
                     ELSE
                        !== stratification is continuous at bottom ==!
                        znsq = pn2(ji,jj,jk)
                     ENDIF
                     !== Burger number (N^2 * delta_z^2)/(f^2 * A) ==!
                     rbu(ji,jj,jk) = ( MAX( znsq, zqglep1 )  * e3t_n(ji,jj,jk)**2 ) /    &
                        &  ( MAX( ff_t(ji,jj)**2, zqglep2 ) * esqt(ji,jj) )
                  END DO
               END DO
            END DO
            !
            CALL lbc_lnk_multi( 'ldfdyn', rro2, 'T', 1. )
            CALL lbc_lnk_multi( 'ldfdyn', rbu, 'T', 1. )
            !
            !== are we in the QG limit? Find the stretching value in x and y components ==!
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     !== are we below the mixed layer and above the sea floor? ==!
                     IF( jk > nmlnqg(ji,jj) .AND. jk < jpkm1 ) THEN
               		   !== x component of stretching ==!
                        zztmpx = MIN( ABS( zstx(ji,jj,jk) ),                                       &
                           &  ABS( ( zwzdx(ji,jj,jk) ) /                                           &
                           &  ( MAX( rbu(ji,jj,jk), rro2(ji,jj,jk) ) + zqglep2 ) ) )
                        tmpzstx(ji,jj,jk) = zwzdx(ji,jj,jk)/                                       &
                           &  ( MAX( rbu(ji,jj,jk), rro2(ji,jj,jk) ) + zqglep2 )
                        zstlimx(ji,jj,jk) = SIGN( zztmpx, zstx(ji,jj,jk) )
                        !== y component of stretching ==!
                        zztmpy = MIN( ABS( zsty(ji,jj,jk) ),                                       &
                           &  ABS( ( zwzdy(ji,jj,jk) ) /                                           &
                           &  ( MAX( rbu(ji,jj,jk), rro2(ji,jj,jk) ) + zqglep2 ) ) )
                        zstlimy(ji,jj,jk) = SIGN( zztmpy, zsty(ji,jj,jk) )
                     ENDIF
                  END DO
               END DO
            END DO
            !
            DO jk = 1, jpkm1                                      !==  Horizontal divergence  ==!
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     hdivnqg(ji,jj,jk) = (  e2u(ji  ,jj) * e3u_n(ji  ,jj,jk) * un(ji  ,jj,jk)    &
                        &               - e2u(ji-1,jj) * e3u_n(ji-1,jj,jk) * un(ji-1,jj,jk)      &
                        &               + e1v(ji,jj  ) * e3v_n(ji,jj  ,jk) * vn(ji,jj  ,jk)      &
                        &               - e1v(ji,jj-1) * e3v_n(ji,jj-1,jk) * vn(ji,jj-1,jk)  )   &
                        &            * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
                  END DO  
               END DO  
            END DO
            !
            !== calculate gradients of divergence, then square of magnitude (f-point) ==!
!            CALL div_hor( kt )
!            !
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     zztmpx = r1_2 * ( r1_e1t(ji,jj+1) * ( hdivnqg(ji+1,jj+1,jk) - hdivnqg(ji,jj+1,jk) )               &
                        &  + r1_e1t(ji,jj) * ( hdivnqg(ji+1,jj,jk) - hdivnqg(ji,jj,jk) ) )
                     zztmpy = r1_2 * ( r1_e2t(ji+1,jj) * ( hdivnqg(ji+1,jj+1,jk) - hdivnqg(ji+1,jj,jk) )               &
                        &  + r1_e2t(ji,jj) * ( hdivnqg(ji,jj+1,jk) - hdivnqg(ji,jj,jk) ) )
                     ddivmagsq(ji,jj,jk) = ( zztmpx * zztmpx + zztmpy * zztmpy ) * fmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
            !== square of magnitude of QG potential vorticity, see Pearson et al. (2017). On t-point ==!
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     zztmpx = zwzdx(ji,jj,jk) + zstlimx(ji,jj,jk)
                     zztmpy = zwzdy(ji,jj,jk) + zstlimy(ji,jj,jk)
                     dwzmagsq(ji,jj,jk) = ( zztmpx * zztmpx + zztmpy * zztmpy ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
            CALL lbc_lnk_multi( 'ldfdyn', dwzmagsq, 'T', 1. )
            !
            !== calculate visacosity coefficient ==!
            DO jk = 1, jpkm1	         !== QG Leith viscosity coefficient on T-point ==!
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1 ! vector opt.
                     zsqqg = r1_4 * ( ddivmagsq(ji,jj,jk) + ddivmagsq(ji-1,jj,jk) + ddivmagsq(ji,jj-1,jk) +     &
                        &  ddivmagsq(ji-1,jj-1,jk) ) + dwzmagsq(ji,jj,jk)
                     ahmt(ji,jj,jk) = SQRT( zcmqgl * esqt(ji,jj)**3 * zsqqg )
                  END DO
               END DO
            END DO
            !
            DO jk = 1, jpkm1            !== QG Leith viscosity coefficient on F-point ==!
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1 ! vector opt.
                     zsqqg = r1_4 * ( dwzmagsq(ji,jj,jk) + dwzmagsq(ji+1,jj,jk) + dwzmagsq(ji,jj+1,jk) +     &
                        &  dwzmagsq(ji+1,jj+1,jk) ) + ddivmagsq(ji,jj,jk)
                     ahmf(ji,jj,jk) = SQRT( zcmqgl * esqf(ji,jj)**3 * zsqqg )
                  END DO
               END DO
            END DO
            !
         ENDIF
         !
         CALL lbc_lnk_multi( 'ldfdyn', ahmt, 'T', 1.,  ahmf, 'F', 1. )
         !
         !== QG Leith diagnostics ==!
         CALL iom_put( "rro2"   , rro2(:,:,:) )      ! square of Rossby number T- point
         CALL iom_put( "rbu"    , rbu(:,:,:) )       ! Burger number T- point
         CALL iom_put( "zstx"   , zstx(:,:,:) )      ! x component of QG stretching T- point
         CALL iom_put( "zsty"   , zsty(:,:,:) )      ! y component of QG stretching T- point
         CALL iom_put( "zstlimx", zstlimx(:,:,:) )   ! x component of QG stretching T- point
         CALL iom_put( "zstlimy", zstlimy(:,:,:) )   ! y component of QG stretching T- point
         CALL iom_put( "zwzdx"  , zwzdx(:,:,:) )     ! x component of vorticity gradient T- point
         CALL iom_put( "zwzdy"  , zwzdy(:,:,:) )     ! y component of vorticity gradient T- point
         CALL iom_put( "zwz"    , zwz(:,:,:) )       ! QG vorticity at F- point
         CALL iom_put( "zbudx"  , zbudx(:,:,:) )     ! x component of buoyancy gradient T- point
         CALL iom_put( "zbudy"  , zbudy(:,:,:) )     ! y component of buoyancy gradient T- point
         CALL iom_put( "tmpzstx", tmpzstx(:,:,:) )   ! temp QG Leith stretching term
         CALL iom_put( "mld_qg" , mld_qg(:,:) )      ! QG Leith mixed layer depth
         !
      END SELECT
      !
      CALL iom_put( "ahmt_2d", ahmt(:,:,1) )   ! surface u-eddy diffusivity coeff.
      CALL iom_put( "ahmf_2d", ahmf(:,:,1) )   ! surface v-eddy diffusivity coeff.
      CALL iom_put( "ahmt_3d", ahmt(:,:,:) )   ! 3D      u-eddy diffusivity coeff.
      CALL iom_put( "ahmf_3d", ahmf(:,:,:) )   ! 3D      v-eddy diffusivity coeff.
      !
      IF( ln_timing )   CALL timing_stop('ldf_dyn')
      !
   END SUBROUTINE ldf_dyn

   !!======================================================================
END MODULE ldfdyn


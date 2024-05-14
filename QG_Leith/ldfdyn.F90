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
   !!   ldf_dyn_str   : update stretching term daily for use in QG Leith
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
   PUBLIC   ldf_dyn_str    ! called by ldf_dyn

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
   REAL(wp), PUBLIC ::   rn_c2dc_vor               !: 2D Leith tuning parameter for vorticity part, typically set to 1
   REAL(wp), PUBLIC ::   rn_c2dc_div               !: 2D Leith tuning parameter for divergence part, typically set to 1
   !                                        ! QG Leith viscosity  (nn_ahm_ijk_t = 34) 
   REAL(wp), PUBLIC ::   rn_cqgc_vor               !: QG Leith tuning parameter for vorticity part, typically set to 1
   REAL(wp), PUBLIC ::   rn_cqgc_div               !: QG Leith tuning parameter for divergence part, typically set to 1
!!   !                                        ! Leith viscosity parameter
!!   REAL(wp), PUBLIC ::   rn_minleith           !: Minimum value used in biharmonic Leith viscosity
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
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dzwzmagsq     !: square of magnitude of gradient of vorticity (2D/QG Leith)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ddivmagsq    !: square of magnitude of gradient of divergence (2D/QG Leith)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hdivnqg      !: Horizontal divergence on t-point (2D/QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zwz          !: Vorticity on f-point (2D/QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zstx, zsty   !: x and y components of stretching at T- points (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zwzdx, zwzdy !: x and y components of horizontal gradients of vertical vorticity at T- points (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hdivdx, hdivdy !: x and y components of horizontal gradients of divergence  at T- points (QG Leith)
!!   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zstlimx, zstlimy !: Limit of stretching term at T- points (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zbu          !: Buoyancy at T- point (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rbu, rro2, rfr2   !: Burger number, square of Rossby, Froude, and grid Reynolds number at T- points
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zbudxup, zbudyvp !: gradients of buoyancy - x and y components on U- point and V- points, resp. (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zbudx, zbudy !: x and y components of gradients in buoyancy at T- points (QG Leith)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tmpzstx      !: alternative stretching term computed in QG Leith for diagnostic purposes 
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   zrho10_3     !: mixed layer depth
   INTEGER,          ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   nmlnqg       !: number of levels in mixed layer
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ahmt_qg, ahmt_div   !: PV and Div contributions to QG Leith T-points [m2/s]

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
!      REAL(wp) ::   zztmp, zztmpx, zztmpy          ! local scalar
      CHARACTER(len=5) ::   cl_Units               ! units (m2/s or m4/s)
      !!
      NAMELIST/namdyn_ldf/ ln_dynldf_OFF, ln_dynldf_lap, ln_dynldf_blp,   &   ! type of operator
         &                 ln_dynldf_lev, ln_dynldf_hor, ln_dynldf_iso,   &   ! acting direction of the operator
         &                 nn_ahm_ijk_t , rn_Uv    , rn_Lv,   rn_ahm_b,   &   ! lateral eddy coefficient
         &                 rn_csmc      , rn_minfac    , rn_maxfac,       &   ! Smagorinsky settings
         &                 rn_c2dc_vor, rn_c2dc_div,                      &   ! 2D Leith tuning parameters
         &                 rn_cqgc_vor, rn_cqgc_div                           ! QG Leith tuning parameters
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
         WRITE(numout,*) '         2D Leith coefficient vorticity       rn_c2dc_vor   = ', rn_c2dc_vor
         WRITE(numout,*) '         2D Leith coefficient divergence      rn_c2dc_div   = ', rn_c2dc_div
         !
         WRITE(numout,*) '      QG Leith settings (nn_ahm_ijk_t  = 34) :'
         WRITE(numout,*) '         QG Leith coefficient vorticity       rn_cqgc_vor   = ', rn_cqgc_vor
         WRITE(numout,*) '         QG Leith coefficient divergence      rn_cqgc_div   = ', rn_cqgc_div
!!         !
!!         WRITE(numout,*) '      Biharmonic Leith settings (nn_ahm_ijk_t  = 33/34) :'
!!         WRITE(numout,*) '         minimum leith viscosity              rn_minleith   = ', rn_minleith
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
            ALLOCATE( esqt(jpi,jpj) , esqf(jpi,jpj) , dzwzmagsq(jpi,jpj,jpk) , ddivmagsq(jpi,jpj,jpk) ,     &
               &  zwz(jpi,jpj,jpk) , hdivnqg(jpi,jpj,jpk) ,                                                 &
               &  hdivdx(jpi,jpj,jpk) , hdivdy(jpi,jpj,jpk) , zwzdx(jpi,jpj,jpk) , zwzdy(jpi,jpj,jpk) ,     &
               &  STAT=ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ldf_dyn_init: failed to allocate 2D Leith arrays')
            !
            DO jj = 1, jpj             ! Set local gridscale values
               DO ji = 1, jpi
                  esqt(ji,jj) = e1e2t(ji,jj)
                  esqf(ji,jj) = e1e2f(ji,jj)
               END DO
            END DO
            !
            !== initialise key variables with zeros (probably don't need to do this??) ==!
            dzwzmagsq(:,:,:) = 0._wp
            ddivmagsq(:,:,:) = 0._wp
            hdivnqg(:,:,:) = 0._wp
            zwz(:,:,:) = 0._wp
!!            rre(:,:,:) = 0._wp
            zwzdx(:,:,:) = 0._wp
            zwzdy(:,:,:) = 0._wp
            hdivdx(:,:,:) = 0._wp
            hdivdy(:,:,:) = 0._wp
            !
         CASE(  34  )       !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( latitude, longitude, depth , time )'
            IF(lwp) WRITE(numout,*) '           proportional to the PV gradient, divergence, and gridscale (QG Leith)'
            !
            l_ldfdyn_time = .TRUE.     ! will be calculated by call to ldf_dyn routine in step.F90
            !
            !                          ! allocate arrays used in ldf_dyn. 
            ALLOCATE( esqt(jpi,jpj) , esqf(jpi,jpj) , dzwzmagsq(jpi,jpj,jpk) , ddivmagsq(jpi,jpj,jpk) ,           &
               &  zbu(jpi,jpj,jpk) , zbudxup(jpi,jpj,jpk) , zbudyvp(jpi,jpj,jpk) , zwz(jpi,jpj,jpk) ,             &
               &  zstx(jpi,jpj,jpk) , zsty(jpi,jpj,jpk) ,                                                         &
               &  rbu(jpi,jpj,jpk), rro2(jpi,jpj,jpk) , zwzdx(jpi,jpj,jpk) , zwzdy(jpi,jpj,jpk) ,                 &
               &  zbudx(jpi,jpj,jpk) , zbudy(jpi,jpj,jpk) , hdivnqg(jpi,jpj,jpk) , rfr2(jpi,jpj,jpk) ,            &
               &  tmpzstx(jpi,jpj,jpk) , zrho10_3(jpi, jpj) , nmlnqg(jpi, jpj) ,                                  &
               &  hdivdx(jpi,jpj,jpk) , hdivdy(jpi,jpj,jpk) , ahmt_qg(jpi,jpj,jpk) , ahmt_div(jpi,jpj,jpk) ,      &
               &  STAT=ierr )
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
            dzwzmagsq(:,:,:) = 0._wp
            ddivmagsq(:,:,:) = 0._wp
            zbu(:,:,:) = 0._wp
            zbudxup(:,:,:) = 0._wp
            zbudyvp(:,:,:) = 0._wp
            zbudx(:,:,:) = 0._wp
            zbudy(:,:,:) = 0._wp
            zstx(:,:,:) = 0._wp
            zsty(:,:,:) = 0._wp
!!            zstlimx(:,:,:) = 0._wp
!!            zstlimy(:,:,:) = 0._wp
            zwz(:,:,:) = 0._wp
            zwzdx(:,:,:) = 0._wp
            zwzdy(:,:,:) = 0._wp
            hdivdx(:,:,:) = 0._wp
            hdivdy(:,:,:) = 0._wp
!!            rre(:,:,:) = 0._wp
            rbu(:,:,:) = 0._wp
            rro2(:,:,:) = 0._wp
            rfr2(:,:,:) = 0._wp
            hdivnqg(:,:,:) = 0._wp
            tmpzstx(:,:,:) = 0._wp
            ahmt_qg(:,:,:) = 0._wp
            ahmt_div(:,:,:) = 0._wp
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
            !
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE ldf_dyn_init


   SUBROUTINE ldf_dyn( kt, kit000, prd, pn2, ahm_leith )
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
      INTEGER, INTENT(in) ::   kt                                             ! time step index
      INTEGER, INTENT(in) ::   kit000                                         ! first time step index
      REAL(wp), INTENT(in),    DIMENSION(:,:,:) ::   prd                      ! before in situ density
      REAL(wp), INTENT(in),    DIMENSION(:,:,:) ::   pn2                      ! before Brunt-Vaisala frequency
      REAL(wp), INTENT(inout), DIMENSION(:,:,:) ::   ahm_leith                ! ahmt for use in ldftra
      REAL(wp)                                  ::   zrho3   = 0.03_wp        ! density     criterion for mixed layer depth
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikt          ! local integer (option 34)
      REAL(wp) ::   zu2pv2_ij_p1, zu2pv2_ij, zu2pv2_ij_m1, zemax   ! local scalar (option 31)
      REAL(wp) ::   zcmsmag, zstabf_lo, zstabf_up, zdelta, zdb     ! local scalar (option 32)
      REAL(wp) ::   zcm2dl, zsq2d                                  ! local scalar (option 33)
      REAL(wp) ::   ahmt_max, ahmf_max, zztmpx, zztmpy             ! local scalar (option 33/34)
      REAL(wp) ::   zcmqgl, zsqqg, zztmp, zzdep, zu                ! local scalar (option 34)
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
         IF( ln_dynldf_lap .OR. ln_dynldf_blp  ) THEN
            !
		      ! allocate local variables !
		      zcm2dl = (1/rpi)**6        			! (1/pi)^6
            zztmp = (rn_c2dc_vor/rpi)**2        ! based on vorticity parameter
            zstabf_lo = rn_minfac * rn_minfac / ( 2._wp * 12._wp * 12._wp * zztmp ) ! lower limit stability factor scaling
            IF( ln_dynldf_blp ) zstabf_lo = ( 16._wp / 9._wp ) * zstabf_lo          ! lower limit biharmonic scaling factor
		      !
		      !== calculate vertical vorticity (f + zeta) on f-point ==!
		      DO jk = 1, jpkm1                                 ! Horizontal slab
		         DO jj = 1, jpjm1
		            DO ji = 1, fs_jpim1   ! vector opt.
		               zwz(ji,jj,jk) = ff_f(ji,jj) + ( e2v(ji+1,jj  ) * vb(ji+1,jj  ,jk) - e2v(ji,jj) * vb(ji,jj,jk)            &
		                  &          - e1u(ji  ,jj+1) * ub(ji  ,jj+1,jk) + e1u(ji,jj) * ub(ji,jj,jk)  ) * r1_e1e2f(ji,jj)
		            END DO
		         END DO
		      END DO
		      !
		      !== calculate gradients of vorticity, then square of magnitude (t-point) ==!
		      DO jk = 1, jpkm1
		         DO jj = 2, jpjm1
		            DO ji = 2, jpim1
		               zztmpx = r1_2 * ( ( r1_e1v(ji,jj-1) * ( zwz(ji,jj-1,jk) - zwz(ji-1,jj-1,jk) ) * vmask(ji  ,jj-1,jk) )            &
		                  &            + ( r1_e1v(ji,jj  ) * ( zwz(ji,jj  ,jk) - zwz(ji-1,jj  ,jk) ) * vmask(ji  ,jj  ,jk) ) )
		               zwzdx(ji,jj,jk) = zztmpx
		               zztmpy = r1_2 * ( ( r1_e2u(ji-1,jj) * ( zwz(ji-1,jj,jk) - zwz(ji-1,jj-1,jk) ) * umask(ji-1,jj  ,jk) )            &
		                  &            + ( r1_e2u(ji  ,jj) * ( zwz(ji  ,jj,jk) - zwz(ji,jj-1  ,jk) ) * umask(jj  ,jj  ,jk) ) )
		               zwzdy(ji,jj,jk) = zztmpy
		               dzwzmagsq(ji,jj,jk) = ( zztmpx * zztmpx + zztmpy * zztmpy )
		            END DO
		         END DO
		      END DO
		      !
		      CALL lbc_lnk_multi( 'ldfdyn', dzwzmagsq, 'T', 1., zwzdx, 'T', 1., zwzdy, 'T', 1. )
		      !
		      DO jk = 1, jpkm1                                      !==  Horizontal divergence  ==!
		         DO jj = 2, jpjm1
		            DO ji = fs_2, fs_jpim1   ! vector opt.
		               hdivnqg(ji,jj,jk) = (  e2u(ji  ,jj) * e3u_b(ji  ,jj,jk) * ub(ji  ,jj,jk)      &
		                  &                 - e2u(ji-1,jj) * e3u_b(ji-1,jj,jk) * ub(ji-1,jj,jk)      &
		                  &                 + e1v(ji,jj  ) * e3v_b(ji,jj  ,jk) * vb(ji,jj  ,jk)      &
		                  &                 - e1v(ji,jj-1) * e3v_b(ji,jj-1,jk) * vb(ji,jj-1,jk)  )   &
		                  &                 * r1_e1e2t(ji,jj) / e3t_b(ji,jj,jk)
		            END DO  
		         END DO  
		      END DO
		      !
		      CALL lbc_lnk_multi( 'ldfdyn', hdivnqg, 'T', 1. )
		      !
		      !== calculate gradients of divergence, then square of magnitude (f-point) ==!
		      DO jk = 1, jpkm1
		         DO jj = 1, jpjm1
		            DO ji = 1, jpim1
		               zztmpx = r1_2 * ( ( r1_e1u(ji,jj+1) * ( hdivnqg(ji+1,jj+1,jk) - hdivnqg(ji,jj+1,jk) ) * umask(ji,jj+1,jk) )               &
		                  &            + ( r1_e1u(ji,jj  ) * ( hdivnqg(ji+1,jj  ,jk) - hdivnqg(ji,jj  ,jk) ) * umask(ji,jj  ,jk) ) )
		               hdivdx(ji,jj,jk) = zztmpx
		               zztmpy = r1_2 * ( ( r1_e2v(ji+1,jj) * ( hdivnqg(ji+1,jj+1,jk) - hdivnqg(ji+1,jj,jk) ) * vmask(ji+1,jj,jk) )               &
		                  &            + ( r1_e2v(ji  ,jj) * ( hdivnqg(ji  ,jj+1,jk) - hdivnqg(ji  ,jj,jk) ) * vmask(jj  ,ji,jk) ) )
		               hdivdy(ji,jj,jk) = zztmpy
		               ddivmagsq(ji,jj,jk) = ( zztmpx * zztmpx + zztmpy * zztmpy )
		            END DO
		         END DO
		      END DO
		      !
		      CALL lbc_lnk_multi( 'ldfdyn', ddivmagsq , 'F', 1., hdivdx, 'F', 1., hdivdy, 'F', 1. )
		      !
		      DO jk = 1, jpkm1	         !== 2D Leith viscosity coefficient on T-point ==!
		         DO jj = 2, jpjm1
		            DO ji = fs_2, fs_jpim1 ! vector opt.
                     !
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zu2pv2_ij_m1 = ub(ji-1,jj  ,jk) * ub(ji-1,jj  ,jk) + vb(ji  ,jj-1,jk) * vb(ji  ,jj-1,jk)
                     !
		               zsq2d = ( rn_c2dc_vor**6 * dzwzmagsq(ji,jj,jk) ) +                                                            &
		                  &    ( rn_c2dc_div**6 * r1_4 * ( ddivmagsq(ji,jj,jk) + ddivmagsq(ji-1,jj,jk) + ddivmagsq(ji,jj-1,jk) +     &
		                  &      ddivmagsq(ji-1,jj-1,jk) ) )
                     !
                     zdelta = (rn_c2dc_vor/rpi)**2 * esqt(ji,jj)
                     !
                     ahmt(ji,jj,jk) = MAX( SQRT( zcm2dl * esqt(ji,jj)**3 * zsq2d, &
                        &                  SQRT( (zu2pv2_ij + zu2pv2_ij_m1) * zdelta * zstabf_lo ) ) ! Impose lower limit
		               ahmt_max = ( MIN( e1t(ji,jj), e2t(ji,jj) )**2 ) / ( 8.0_wp * rn_rdt )  
		               ahmt(ji,jj,jk) = MIN( ahmt(ji,jj,jk) , ahmt_max ) ! impose upper limit
		            END DO
		         END DO
		      END DO
		      !
		      DO jk = 1, jpkm1            !== 2D Leith viscosity coefficient on F-point ==!
		         DO jj = 1, jpjm1
		            DO ji = 1, fs_jpim1 ! vector opt.
                     !
                     zu2pv2_ij_p1 = ub(ji  ,jj+1,jk) * ub(ji  ,jj+1,jk) + vb(ji+1,jj  ,jk) * vb(ji+1,jj  ,jk)
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     !
		               zsq2d = ( rn_c2dc_vor**6 * r1_4 * ( dzwzmagsq(ji,jj,jk) + dzwzmagsq(ji+1,jj,jk) + dzwzmagsq(ji,jj+1,jk) +     &
		                  &  dzwzmagsq(ji+1,jj+1,jk) ) ) + ( rn_c2dc_div**6 * ddivmagsq(ji,jj,jk) )
                     !
                     zdelta = (rn_c2dc_vor/rpi)**2 * esqf(ji,jj)
                     !
                     ahmf(ji,jj,jk) = MAX( SQRT( zcm2dl * esqf(ji,jj)**3 * zsq2d, &
                        &                  SQRT( (zu2pv2_ij_p1 + zu2pv2_ij) * zdelta * zstabf_lo ) ) ! Impose lower limit
		               ahmf_max = ( MIN( e1f(ji,jj), e2f(ji,jj) )**2 ) / ( 8.0_wp * rn_rdt )  
		               ahmf(ji,jj,jk) = MIN( ahmf(ji,jj,jk) , ahmf_max ) ! impose upper limit
		            END DO
		         END DO
		      END DO
		      !
         ENDIF
         !
         IF( ln_dynldf_blp ) THEN ! bilaplacian operator, ahm_lap * delta^2 / 8 (Griffies and Hallberg, 2000)
            !
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ahmt(ji,jj,jk) = SQRT( r1_8 * ahmt(ji,jj,jk) * MIN( e1t(ji,jj), e2t(ji,jj) )**2 ) 
                  END DO
               END DO
               !
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     ahmf(ji,jj,jk) = SQRT( r1_8 * ahmf(ji,jj,jk) * MIN( e1f(ji,jj), e2f(ji,jj) )**2 )
                  END DO
               END DO
            END DO
            !
         ENDIF
         !
         CALL lbc_lnk_multi( 'ldfdyn', ahmt, 'T', 1.,  ahmf, 'F', 1. )
         !
         !== assigning for output and use in step.f90 ==!
         ahm_leith(:,:,:) = ahmt(:,:,:) ! can this be assigned only when GM/Redi turned on?
         !
         !== 2D Leith diagnostics ==!
         CALL iom_put( "zwzdx"   , zwzdx(:,:,:) )     ! x component of vorticity gradient T- point
         CALL iom_put( "zwzdy"   , zwzdy(:,:,:) )     ! y component of vorticity gradient T- point
         CALL iom_put( "hdivdx"  , hdivdx(:,:,:) )    ! x component of divergence gradient T- point
         CALL iom_put( "hdivdy"  , hdivdy(:,:,:) )    ! y component of divergence gradient T- point
         !
         !
      CASE(  34  )       !==  time varying 3D field  ==!   = F(QG PV gradient, divergence, and gridscale) (QG Leith)
         !
         IF( ln_dynldf_lap .OR. ln_dynldf_blp  ) THEN
            !
		      ! allocate local variables !
		      zcmqgl = (1/rpi)**6         			! (1/pi)^6
		      !
		      !== Compute the mixed layer depth based on a density criteria of zrho = 0.03 (see diahth.F90) ==!
		      ! initialization
		      zrho3 = 0.03_wp
		      DO jj = 1, jpj
		         DO ji = 1, jpi
		            nmlnqg(ji,jj) = mbkt(ji,jj)           ! Initialization to the number of T ocean points
		            zztmp = gdepw_b(ji,jj,mbkt(ji,jj)+1)
		            zrho10_3(ji,jj) = zztmp
		         END DO
		      END DO 
		      !
		      ! ------------------------- !
		      ! MLD: rho = rho10m + zrho3 !
		      ! ------------------------- !
		      DO jk = jpkm1, nlb10, -1    ! loop from bottom to nlb10
		         DO jj = 1, jpj
		            DO ji = 1, jpi
		               ikt = mbkt(ji,jj)
		               zzdep = gdepw_b(ji,jj,jk) * tmask(ji,jj,1)
		               zztmp = rhop(ji,jj,jk) - rhop(ji,jj,nla10)              ! delta rho(10m)
		               IF( zztmp > zrho3 ) THEN
		                  zrho10_3(ji,jj) = zzdep                              ! > 0.03
		                  nmlnqg(ji,jj) = MIN(jk, ikt) + 1                     ! Mixed layer level
		               ENDIF
		            END DO
		         END DO
		      END DO
		      !
		      !== calculate vertical vorticity (f+zeta) on f-point ==!
		      DO jk = 1, jpkm1                                 ! Horizontal slab
		         DO jj = 1, jpjm1
		            DO ji = 1, fs_jpim1   ! vector opt.
		               zwz(ji,jj,jk) = ff_f(ji,jj) + (  e2v(ji+1,jj  ) * vb(ji+1,jj  ,jk) - e2v(ji,jj) * vb(ji,jj,jk)            &
		                  &          - e1u(ji  ,jj+1) * ub(ji  ,jj+1,jk) + e1u(ji,jj) * ub(ji,jj,jk)  ) * r1_e1e2f(ji,jj)
		            END DO
		         END DO
		      END DO
		      !
		      !== calculate gradients of vorticity, then square of magnitude (t-point) ==!
		      DO jk = 1, jpkm1
		         DO jj = 2, jpjm1
		            DO ji = 2, jpim1
		               zwzdx(ji,jj,jk) = r1_2 * ( ( r1_e1v(ji,jj-1) * ( zwz(ji,jj-1,jk) - zwz(ji-1,jj-1,jk) ) * vmask(ji  ,jj-1,jk) )            &
		                  &                     + ( r1_e1v(ji,jj  ) * ( zwz(ji,jj  ,jk) - zwz(ji-1,jj  ,jk) ) * vmask(ji  ,jj  ,jk) ) )
		               zwzdy(ji,jj,jk) = r1_2 * ( ( r1_e2u(ji-1,jj) * ( zwz(ji-1,jj,jk) - zwz(ji-1,jj-1,jk) ) * umask(ji-1,jj  ,jk) )            &
		                  &                     + ( r1_e2u(ji  ,jj) * ( zwz(ji  ,jj,jk) - zwz(ji,jj-1  ,jk) ) * umask(jj  ,jj  ,jk) ) )
		            END DO
		         END DO
		      END DO
		      !
		      CALL lbc_lnk_multi( 'ldfdyn', zwzdx, 'T', 1., zwzdy, 'T', 1. )
		      !
		      !== Compute stretching term at first time step index and then at specified intervals ** runs at every timestep ** ==!
		      IF( kt == kit000 ) THEN       !! compute stretching subroutine
		         !
		         zstlimx(:,:,:) = 0._wp
		         zstlimy(:,:,:) = 0._wp
		         !
		         CALL ldf_dyn_str( kt, prd, pn2, zwzdx, zwzdy, nmlnqg, zstlimx, zstlimy )
		         !
		      ELSEIF( MOD(kt-1,1) == 0 ) THEN !! need to adjust this to account for user defined timesteps per day. See if it works first.
		         !
		         CALL ldf_dyn_str( kt, prd, pn2, zwzdx, zwzdy, nmlnqg, zstlimx, zstlimy )
		         !
		      ENDIF
		      !
		      DO jk = 1, jpkm1                                      !==  Horizontal divergence  ==!
		         DO jj = 2, jpjm1
		            DO ji = fs_2, fs_jpim1   ! vector opt.
		               hdivnqg(ji,jj,jk) = (  e2u(ji  ,jj) * e3u_b(ji  ,jj,jk) * ub(ji  ,jj,jk)      &
		                  &                 - e2u(ji-1,jj) * e3u_b(ji-1,jj,jk) * ub(ji-1,jj,jk)      &
		                  &                 + e1v(ji,jj  ) * e3v_b(ji,jj  ,jk) * vb(ji,jj  ,jk)      &
		                  &                 - e1v(ji,jj-1) * e3v_b(ji,jj-1,jk) * vb(ji,jj-1,jk)  )   &
		                  &                 * r1_e1e2t(ji,jj) / e3t_b(ji,jj,jk)
		            END DO  
		         END DO  
		      END DO
		      !
		      CALL lbc_lnk_multi( 'ldfdyn', hdivnqg, 'T', 1. )
		      !
		      !== calculate gradients of divergence, then square of magnitude (f-point) ==!
		      DO jk = 1, jpkm1
		         DO jj = 1, jpjm1
		            DO ji = 1, jpim1
		               zztmpx = r1_2 * ( ( r1_e1u(ji,jj+1) * ( hdivnqg(ji+1,jj+1,jk) - hdivnqg(ji,jj+1,jk) ) * umask(ji,jj+1,jk) )               &
		                  &            + ( r1_e1u(ji,jj  ) * ( hdivnqg(ji+1,jj  ,jk) - hdivnqg(ji,jj  ,jk) ) * umask(ji,jj  ,jk) ) ) * fmask(ji,jj,jk)
		               hdivdx(ji,jj,jk) = zztmpx
		               zztmpy = r1_2 * ( ( r1_e2v(ji+1,jj) * ( hdivnqg(ji+1,jj+1,jk) - hdivnqg(ji+1,jj,jk) ) * vmask(ji+1,jj,jk) )               &
		                  &            + ( r1_e2v(ji  ,jj) * ( hdivnqg(ji  ,jj+1,jk) - hdivnqg(ji  ,jj,jk) ) * vmask(jj  ,ji,jk) ) ) * fmask(ji,jj,jk)
		               hdivdy(ji,jj,jk) = zztmpy
		               ddivmagsq(ji,jj,jk) = ( zztmpx * zztmpx + zztmpy * zztmpy )
		            END DO
		         END DO
		      END DO
		      !
		      CALL lbc_lnk_multi( 'ldfdyn', ddivmagsq , 'F', 1. )
		      !
		      !== square of magnitude of QG potential vorticity, see Pearson et al. (2017). On t-point ==!
		      DO jk = 1, jpkm1
		         DO jj = 1, jpj
		            DO ji = 1, jpi
		               zztmpx = zwzdx(ji,jj,jk) + zstlimx(ji,jj,jk)
		               zztmpy = zwzdy(ji,jj,jk) + zstlimy(ji,jj,jk)
		               dzwzmagsq(ji,jj,jk) = ( zztmpx * zztmpx + zztmpy * zztmpy )
		            END DO
		         END DO
		      END DO
		      !
		      !== calculate viscosity coefficient ==!
		      DO jk = 1, jpkm1	         !== QG Leith viscosity coefficient on T-point ==!
		         DO jj = 2, jpjm1
		            DO ji = fs_2, fs_jpim1 ! vector opt.
		               ahmt_qg(ji,jj,jk) = dzwzmagsq(ji,jj,jk)
		               ahmt_div(ji,jj,jk) = r1_4 * ( ddivmagsq(ji,jj,jk) + ddivmagsq(ji-1,jj,jk) + ddivmagsq(ji,jj-1,jk) +     &
		                  &  ddivmagsq(ji-1,jj-1,jk) )
		               !== Set max value on viscosity coefficient ==!
		               zsqqg = ( rn_cqgc_vor**6 * dzwzmagsq(ji,jj,jk) ) +                                                            &
		                  &    ( rn_cqgc_div**6 * r1_4 * ( ddivmagsq(ji,jj,jk) + ddivmagsq(ji-1,jj,jk) + ddivmagsq(ji,jj-1,jk) +     &
		                  &      ddivmagsq(ji-1,jj-1,jk) ) )
		               ahmt_max = ( MIN( e1t(ji,jj), e2t(ji,jj) )**2 ) / ( 8.0_wp * rn_rdt )  
		               ahmt(ji,jj,jk) = MIN( SQRT( zcmqgl * esqt(ji,jj)**3 * zsqqg ), ahmt_max )
		            END DO
		         END DO
		      END DO
		      !
		      CALL lbc_lnk_multi( 'ldfdyn', ahmt_qg, 'T', 1.,  ahmt_div, 'T', 1. )
		      !
		      DO jk = 1, jpkm1            !== QG Leith viscosity coefficient on F-point ==!
		         DO jj = 1, jpjm1
		            DO ji = 1, fs_jpim1 ! vector opt.
		               !== Set max value of viscosity coefficient depending on stability criterion (Stevens, 1995) ==!
		               zsqqg = ( rn_cqgc_vor**6 * r1_4 * ( dzwzmagsq(ji,jj,jk) + dzwzmagsq(ji+1,jj,jk) + dzwzmagsq(ji,jj+1,jk) +     &
		                  &  dzwzmagsq(ji+1,jj+1,jk) ) ) + ( rn_cqgc_div**6 * ddivmagsq(ji,jj,jk) )
		               ahmf_max = ( MIN( e1f(ji,jj), e2f(ji,jj) )**2 ) / ( 8.0_wp * rn_rdt )  
		               ahmf(ji,jj,jk) = MIN( SQRT( zcmqgl * esqf(ji,jj)**3 * zsqqg ), ahmf_max )
		            END DO
		         END DO
		      END DO
		      !
         ENDIF
         !
         IF( ln_dynldf_blp ) THEN ! bilaplacian operator, ahm_lap * delta^2 / 8 (Griffies and Hallberg, 2000)
            !
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     !== Ensuring the viscosity never gets too small, needs to be made grid aware though ==!
!!                     ahmt(ji,jj,jk) = SQRT( MAX( r1_8 * ahmt(ji,jj,jk) * MIN( e1t(ji,jj), e2t(ji,jj) )**2, rn_minleith ) )
                     ahmt(ji,jj,jk) = SQRT( r1_8 * ahmt(ji,jj,jk) * MIN( e1t(ji,jj), e2t(ji,jj) )**2 ) 
                  END DO
               END DO
               !
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     !== Ensuring the viscosity never gets too small, needs to be made grid aware though ==!
!!                     ahmf(ji,jj,jk) = SQRT( MAX( r1_8 * ahmf(ji,jj,jk) * MIN( e1f(ji,jj), e2f(ji,jj) )**2, rn_minleith ) )
                     ahmf(ji,jj,jk) = SQRT( r1_8 * ahmf(ji,jj,jk) * MIN( e1f(ji,jj), e2f(ji,jj) )**2 )
                  END DO
               END DO
            END DO
            !
         ENDIF
         !
         CALL lbc_lnk_multi( 'ldfdyn', ahmt, 'T', 1.,  ahmf, 'F', 1. )
         !
         !== assigning for output and use in step.f90 ==!
         ahm_leith(:,:,:) = ahmt(:,:,:)
         !
         !== QG Leith diagnostics ==!
         CALL iom_put( "rro2"    , rro2(:,:,:) )      ! square of Rossby number T- point
         CALL iom_put( "rbu"     , rbu(:,:,:) )       ! Burger number T- point
         CALL iom_put( "rfr2"    , rfr2(:,:,:) )      ! square of Froude number T- point
         CALL iom_put( "zstx"    , zstx(:,:,:) )      ! x component of QG stretching T- point
         CALL iom_put( "zsty"    , zsty(:,:,:) )      ! y component of QG stretching T- point
         CALL iom_put( "zstlimx" , zstlimx(:,:,:) )   ! x component of QG stretching T- point
         CALL iom_put( "zstlimy" , zstlimy(:,:,:) )   ! y component of QG stretching T- point
         CALL iom_put( "zwzdx"   , zwzdx(:,:,:) )     ! x component of vorticity gradient T- point
         CALL iom_put( "zwzdy"   , zwzdy(:,:,:) )     ! y component of vorticity gradient T- point
         CALL iom_put( "hdivdx"  , hdivdx(:,:,:) )    ! x component of divergence gradient T- point
         CALL iom_put( "hdivdy"  , hdivdy(:,:,:) )    ! y component of divergence gradient T- point
         CALL iom_put( "zwz"     , zwz(:,:,:) )       ! QG vorticity at F- point
         CALL iom_put( "zbudx"   , zbudx(:,:,:) )     ! x component of buoyancy gradient T- point
         CALL iom_put( "zbudy"   , zbudy(:,:,:) )     ! y component of buoyancy gradient T- point
!!         CALL iom_put( "tmpzstx" , tmpzstx(:,:,:) )   ! temp QG Leith stretching term
         CALL iom_put( "ahmt_qg" , ahmt_qg(:,:,:) )   ! PV contribution to QG Leith T-point
         CALL iom_put( "ahmt_div", ahmt_div(:,:,:) )  ! Div contribution to QG Leith T-point
         !
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
   
   
   SUBROUTINE ldf_dyn_str( kt, prd, pn2, zwzdx, zwzdy, nmlnqg, zstlimx, zstlimy )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn_str  ***
      !! 
      !! ** Purpose :   compute stretching term for QG Leith viscosity
      !!
      !! ** Method  :   Input instantaneous density (prd), square of buoyancy frequency (pn2), and 
      !!                gradients of vorticity (zwz + f), then compute stretching term. 
      !!                The stretching term is then stored and remains constant for one desired timesteps.
      !! ** note    :    
      !! ** action  :   zstlimx, zstlimy updated per user requirements.
      !!----------------------------------------------------------------------
      INTEGER,  INTENT(in) ::   kt   ! time step index
      REAL(wp), INTENT(in),  DIMENSION(:,:,:) ::   prd                         ! before in situ density
      REAL(wp), INTENT(in),  DIMENSION(:,:,:) ::   pn2                         ! before Brunt-Vaisala frequency
      REAL(wp), INTENT(in),  DIMENSION(:,:,:) ::   zwzdx                       ! Horizontal vorticity gradient in x-direction
      REAL(wp), INTENT(in),  DIMENSION(:,:,:) ::   zwzdy                       ! Horizontal vorticity gradient in y-direction
      INTEGER,  INTENT(in),  DIMENSION(:,:)   ::   nmlnqg                      ! Level of mixed layer depth
      REAL(wp), INTENT(out), DIMENSION(:,:,:) ::   zstlimx                     ! Stretching in x direction
      REAL(wp), INTENT(out), DIMENSION(:,:,:) ::   zstlimy                     ! Stretching in y direction
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zbuup, zbulw, zusq, znsq, zztmpx, zztmpy      ! local scalar (option 34)
      REAL(wp) ::   zker1, zker2, zqglep1, zqglep2, zztmp, zbeta  ! more local scalar (option 34)
      REAL(wp), DIMENSION(jpi,jpj) ::   zwrk_2d                   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      zqglep1 = 1.e-12_wp
      zqglep2 = 1.e-24_wp
      !
      !== begin calculation of stretching term d/dz[(f/(N**2))*grad(b)] ==!
      !== find buoyancy and interpolate onto w-grid ==!
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( jk < 2 ) THEN
                  !== buoyancy at surface ==!
                  zbu(ji,jj,jk) = - grav * prd(ji,jj,jk)
               ELSE
                  !== buoyancy below surface ==!
                  zbuup = - grav * prd(ji,jj,jk-1)
                  zbulw = - grav * prd(ji,jj,jk  )
                  zbu(ji,jj,jk) = 0.5_wp * ( zbuup + zbulw ) * wmask(ji,jj,jk)
               ENDIF!

            END DO
         END DO
      END DO
      !
      !== Calculate horizontal gradients of buoyancy and put on w-grid ==!
      DO jk = 1, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               !== gradients of buoyancy on U and V grid at w point of cell ==!
               zbudxup(ji,jj,jk) = r1_e1u(ji,jj) * ( zbu(ji+1,jj,jk) - zbu(ji,jj,jk) ) * umask(ji,jj,jk)
               zbudyvp(ji,jj,jk) = r1_e2v(ji,jj) * ( zbu(ji,jj+1,jk) - zbu(ji,jj,jk) ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      CALL lbc_lnk_multi( 'ldfdyn', zbudxup, 'U', 1., zbudyvp, 'V', 1. )
      !
      !== gradients of buoyancy on W- points ==!
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zbudx(ji,jj,jk) = r1_2 * ( zbudxup(ji-1,jj,jk) + zbudxup(ji,jj,jk) ) * wmask(ji,jj,jk)
               zbudy(ji,jj,jk) = r1_2 * ( zbudyvp(ji,jj-1,jk) + zbudyvp(ji,jj,jk) ) * wmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      CALL lbc_lnk_multi( 'ldfdyn', zbudx, 'T', 1., zbudy, 'T', 1.  )
      !
      !== take vertical gradient and find stretching d/dz[(f * grad(b))/N^2] (t-point) ==!
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               !== are we below the mixed layer and above the sea floor? ==!
               IF( jk > nmlnqg(ji,jj) .AND. jk < mbkt(ji,jj)  ) THEN
                  !== vertical gradient of x component ==!
                  zker1 = ( ff_t(ji,jj) * zbudx(ji,jj,jk  ) ) / MAX( pn2(ji,jj,jk  ), zqglep1 ) 
                  zker2 = ( ff_t(ji,jj) * zbudx(ji,jj,jk+1) ) / MAX( pn2(ji,jj,jk+1), zqglep1 ) 
                  zstx(ji,jj,jk) = ( ( zker1 - zker2 ) / e3t_b(ji,jj,jk) ) * tmask(ji,jj,jk)
                  !== vertical gradient of y component ==!
                  zker1 = ( ff_t(ji,jj) * zbudy(ji,jj,jk  ) ) / MAX( pn2(ji,jj,jk  ), zqglep1 )
                  zker2 = ( ff_t(ji,jj) * zbudy(ji,jj,jk+1) ) / MAX( pn2(ji,jj,jk+1), zqglep1 )
                  zsty(ji,jj,jk) = ( ( zker1 - zker2 ) / e3t_b(ji,jj,jk) ) * tmask(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO
      !
      !== calculate the first baroclinic deformation radius Ld = (1/abs(f)*pi)*sum(Nsq)*delta_z ==!
      ! testing in idealised configuration which does not pass the equator
      zwrk_2d(:,:) = 0._wp
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               !== averaging square of buoyancy frequency onto t-grid ==!
               IF( jk < mbkt(ji,jj) ) THEN
                  !== accounting for negative N^2 ==!
                  zztmp = r1_2 * ( MAX( pn2(ji,jj,jk), zqglep1 ) + MAX( pn2(ji,jj,jk+1), zqglep1 ) ) * tmask(ji,jj,jk)
               ELSE
                  !== stratification is continuous at bottom ==!
                  zztmp = MAX( pn2(ji,jj,jk ), zqglep1 ) * tmask(ji,jj,jk)
               ENDIF
               !== first baroclinic deformation radius, Ld (Chelton et al., 1998) ==!
               IF( gphit(ji,jj) < -5 .OR. gphit(ji,jj) > 5 ) THEN ! outside equator
                  zwrk_2d(ji,jj) = zwrk_2d(ji,jj) +    &
                    &         ( SQRT( zztmp ) * e3t_b(ji,jj,jk) ) / ( ABS(ff_t(ji,jj)) * rpi )
               ELSE  ! near the equator, see Gill (1982) - equatorial radius of deformation.
                  zbeta = 2. * omega * COS( rad * gphit(ji,jj) ) / ra
                  zwrk_2d(ji,jj) = zwrk_2d(ji,jj) +    &
                    &         SQRT( ( SQRT( zztmp ) * e3t_b(ji,jj,jk) ) / ( 2 *  zbeta * rpi ) )
               ENDIF
            END DO
         END DO
      END DO
      !
      !== calculate the Burger number and square of Rossby number on t-point ==!
      !== calculate over entire domain for diagnostics ==!
      DO jk = 1, jpkm1
         DO jj = 2, jpj
            DO ji = 2, jpi
               !== grid scale velocity squared ==!
               zusq = 0.5_wp * ( ( ub(ji-1,jj  ,jk) * ub(ji-1,jj  ,jk) + ub(ji,jj,jk) * ub(ji,jj,jk) ) +                &
                  &              ( vb(ji  ,jj-1,jk) * vb(ji  ,jj-1,jk) + vb(ji,jj,jk) * vb(ji,jj,jk) ) )
               !== square of Rossby number U^2/(f^2 * A) ==!
               rro2(ji,jj,jk) = ( zusq / ( MAX( ff_t(ji,jj)**2, zqglep2 ) * esqt(ji,jj) ) ) * tmask(ji,jj,jk)
               !== averaging square of buoyancy frequency onto t-grid ==!
               !IF( jk < mbkt(ji,jj) ) THEN
               !   !== accounting for negative N^2 ==!
               !   znsq = r1_2 * ( MAX( pn2(ji,jj,jk), zqglep1 ) + MAX( pn2(ji,jj,jk+1), zqglep1 ) ) * tmask(ji,jj,jk)
               !ELSE
               !   !== stratification is continuous at bottom ==!
               !   znsq = MAX( pn2(ji,jj,jk ), zqglep1 )
               !ENDIF
               !== Burger number (N^2 * delta_z^2)/(f^2 * A) ==!
               !rbu(ji,jj,jk) = ( MAX( znsq          , zqglep1 ) * e3t_b(ji,jj,jk)**2 ) /    &
               !   &            ( MAX( ff_t(ji,jj)**2, zqglep2 ) * esqt(ji,jj) )
               !== computing Burger number with deformation radius ==!
               rbu(ji,jj,jk) = ( zwrk_2d(ji,jj) * zwrk_2d(ji,jj) ) / esqt(ji,jj)
               !== Froude number squared (Fr^2 = Ro^2/Bu) ==!
               rfr2(ji,jj,jk) = rro2(ji,jj,jk)/rbu(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      CALL lbc_lnk_multi( 'ldfdyn', rro2, 'T', 1., rbu, 'T', 1., rfr2, 'T', 1. )
      !
      !== are we in the QG limit? Find the stretching value in x and y components ==!
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               !== are we below the mixed layer and above the sea floor? ==!
               IF( jk > nmlnqg(ji,jj) .AND. jk < mbkt(ji,jj) ) THEN
                  !== x component of stretching ==!
                  zztmpx = MIN( ABS( zstx(ji,jj,jk) ),                                       &
                     &  ABS( ( zwzdx(ji,jj,jk) * rfr2(ji,jj,jk) ) /                          &
                     &  ( rro2(ji,jj,jk) + rfr2(ji,jj,jk)**2 + zqglep2 ) ) )
                  zstlimx(ji,jj,jk) = SIGN( zztmpx, zstx(ji,jj,jk) )
                  !== y component of stretching ==!
                  zztmpy = MIN( ABS( zsty(ji,jj,jk) ),                                       &
                     &  ABS( ( zwzdy(ji,jj,jk) * rfr2(ji,jj,jk) ) /                          &
                     &  ( rro2(ji,jj,jk) + rfr2(ji,jj,jk)**2 + zqglep2 ) ) )
                  zstlimy(ji,jj,jk) = SIGN( zztmpy, zsty(ji,jj,jk) )
               ENDIF
            END DO
         END DO
      END DO
      !
   END SUBROUTINE ldf_dyn_str

   !!======================================================================
END MODULE ldfdyn


MODULE north_bdy

   USE utils
   ! USE netcdf

   IMPLICIT NONE
   PUBLIC

   CONTAINS

   SUBROUTINE nbdy_resto( presto, gphit )
      !!---------------------------------
      !!     **ROUTINE: nbdy_resto
      !!
      !!  ** Purpose: To calculate a restoring boundary coefficient at the 
      !!					 northern boundary of a Southern Ocean channel model, IDEAL.
      !!
      !!	 ** Method: Take in basic value of resto from `make_dmp_file.F90` and pass it in here.
      !!					Then multiply it by a function that decreases from the northern boundary,
      !!					giving an approximate sponge layer of a specified distance, likely 100 km.
      !!					See pseudo code in notebook.
      !!
      !!-------------------------------------
      IMPLICIT NONE
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: presto
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: gphit
      REAL(wp) :: zlat0 = 100._wp									! width of sponge layer in km
      REAL(wp) :: zlat, zlatmax, zlatstart
      
      INTEGER :: jj, ji           ! dummy loop indices
      
      ! CALL check_nf90( nf90_get_var( ncin, gphit_id, gphit, (/ 1,1 /), (/ jpi, jpj /) ) ) 
      
      zlatmax = maxval(gphit(:,:)) 			! max latitude in the domain in km				
      zlatstart = zlatmax - zlat0		! start of the sponge layer in km
      DO jj = 1, jpj
      	DO ji = 1, jpi
      		zlat = ABS(gphit(ji,jj))
      		IF ( zlat >= zlatstart ) THEN
      			presto(ji,jj) = presto(ji,jj) * (1._wp - COS( 0.5_wp * rpi * ( zlat-zlatstart )/zlat0 ) )
      		ELSEIF ( zlat < zlatstart ) THEN
      			presto(ji,jj) = 0._wp
      		ENDIF
   		END DO
		END DO
      		

   END SUBROUTINE nbdy_resto

END MODULE north_bdy

! simulation of the ground fire by the prescribed heat flux (W/m2) and BC, OC, 
! SO2 emissions (kg) from the surface
! adapted from SFIRE model (see /WRF/phys/module_fr_fire*.F)
! ukhova 110221
module module_kuwait_fire_driver

implicit none
contains

subroutine fire_driver_kuwait_em_step (grid , config_flags          &
            ,fire_year,fire_month,fire_day,fire_hour                &
            ,fire_dur,fire_hflx,fire_bc_emis,fire_oc_emis           &
            ,fire_so2_emis                                          &
            ,curr_secs_in,gmt,julday_in                             &
            ,ids,ide, kds,kde, jds,jde                              &
            ,ims,ime, kms,kme, jms,jme                              &
            ,ips,ipe, kps,kpe, jps,jpe                              &
            ,rho,z_at_w,dz8w) 

    USE module_domain, only: domain
    USE module_configure , only : grid_config_rec_type
    implicit none

    TYPE(domain) , TARGET          :: grid   ! data
    TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags
    integer, intent(in):: &
             ids,ide, kds,kde, jds,jde                              &
            ,ims,ime, kms,kme, jms,jme                              &
            ,ips,ipe, kps,kpe, jps,jpe
    real,dimension(ims:ime, kms:kme, jms:jme),intent(in)::rho,dz8w,z_at_w

    REAL, INTENT(IN) :: curr_secs_in  !current forecast time (seconds)
    REAL, INTENT(IN) :: gmt
    REAL, INTENT(IN) :: julday_in


    REAL, DIMENSION( ims:ime, jms:jme),INTENT(IN ) :: fire_year,fire_month &
              ,fire_day,fire_hour,fire_dur,fire_hflx,fire_bc_emis,fire_oc_emis &
              ,fire_so2_emis

    ! local
    integer :: its,ite,jts,jte,kts,kte            ! atm tile
    integer:: ij

    INTEGER :: dbg_lvl = 1

    CALL wrf_debug ( dbg_lvl , 'fire_driver_kuwait_em_step: FIRE step start' )

    call fire_driver_kuwait_phys (&
          ids,ide, kds,kde, jds,jde,                  &
          ims,ime, kms,kme, jms,jme,                  &
          ips,ipe, kps,kpe, jps,jpe,                  &
          fire_year,fire_month,fire_day,fire_hour,    &
          fire_dur,fire_hflx,fire_bc_emis,fire_so2_emis,&
          fire_oc_emis,                               &
          curr_secs_in,gmt,julday_in,                 &
          rho,dz8w,                                   &
          grid%chem,                                  &            
          grid%num_tiles,                             & ! atm grid tiling
          grid%i_start,min(grid%i_end,ide-1),                    &
          grid%j_start,min(grid%j_end,jde-1),                    &
          grid%dt,grid%dx,grid%dy,                    &
          grid%msftx,grid%msfty,              &
          grid%grnhfx,                    & ! out redundant arrays, atm grid
          grid, config_flags)

    
    ! --- add heat flux to tendency variables by postulated decay
    do ij=1,grid%num_tiles
       ! FIRE works on domain by 1 smaller, in last row&col winds are not set properly
       its = grid%i_start(ij)             ! start atmospheric tile in i
       ite = min(grid%i_end(ij),ide-1)    ! end atmospheric tile in i
       jts = grid%j_start(ij)             ! start atmospheric tile in j
       jte = min(grid%j_end(ij),jde-1)    ! end atmospheric tile in j
       kts=kds
       kte=kde

       call kuwait_fire_tendency(                 &
            ids,ide-1, kds,kde, jds,jde-1, & ! domain dimensions
            ims,ime, kms,kme, jms,jme,      &
            its,ite, kts,kte, jts,jte,      & ! 
            grid%grnhfx,                    & ! fluxes on atm grid 
            grid%ht,z_at_w,dz8w,grid%mut,grid%c1h,grid%c2h,rho,&
            grid%rthfrten)                    ! out
    enddo


    CALL wrf_debug ( dbg_lvl , 'fire_driver_kuwait_em_step: FIRE step complete' )
            
end subroutine fire_driver_kuwait_em_step

subroutine fire_driver_kuwait_phys (&
    ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
    ims,ime, kms,kme, jms,jme,                    &
    ips,ipe, kps,kpe, jps,jpe,                    &
    fire_year,fire_month,fire_day,fire_hour,      &
    fire_dur,fire_hflx,fire_bc_emis,fire_so2_emis,&
    fire_oc_emis,                                 &
    curr_secs,gmt,julday,                         &
    rho,dz8w,                                     &
    chem,                                         &
    num_tiles,i_start,i_end,j_start,j_end,        & ! atm grid tiling
    dt,dx,dy,                                     & ! in scalars
    mapfac_mx, mapfac_my,                         &
    grnhfx,                                       & ! out
    grid, config_flags)

    USE module_domain, only: domain
    USE module_configure, only: grid_config_rec_type

    USE module_state_description
    USE module_data_mosaic_asect

    implicit none

        type(domain) , target :: grid                   ! state 

    ! --- incoming variables
    integer, intent(in)::                        &
        ids,ide, kds,kde, jds,jde,                    & ! atm domain bounds
        ims,ime, kms,kme, jms,jme,                    & ! atm memory bounds 
        ips,ipe, kps,kpe, jps,jpe,                    & ! atm patch bounds
        num_tiles                                       ! number of tiles

    REAL, DIMENSION( ims:ime, jms:jme),INTENT(IN ) :: fire_year,fire_month &
              ,fire_day,fire_hour,fire_dur,fire_hflx,fire_bc_emis,fire_oc_emis &
              ,fire_so2_emis

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem

    real,dimension(ims:ime, kms:kme, jms:jme),intent(in)::rho,dz8w

    type (grid_config_rec_type), intent(in)          :: config_flags 

    integer,dimension(num_tiles),intent(in) :: i_start,i_end,j_start,j_end  ! atm grid tiling

    real, intent(IN):: dt,dx,dy              ! start of this time step from the beginning of the simulation
    REAL, INTENT(IN) :: curr_secs
    REAL, INTENT(IN) :: gmt
    REAL, INTENT(IN) :: julday


    real,intent(in),dimension(ims:ime, jms:jme)::   mapfac_mx, mapfac_my    ! map factors
        
    ! --- outgoing variables
    real, intent(out), dimension(ims:ime, jms:jme):: grnhfx ! heat flux from surface (W/m^2) 

    ! --- local variables
    integer :: its,ite,jts,jte,ij,i,j
    INTEGER :: endday
    real*8, dimension (1,1) ::airmas,area

    integer :: beg_jul

    REAL :: gmtp,gmtm
    integer :: curr_hours,ifire
    Integer :: endhr,beghr,K,end_hour
    CHARACTER (LEN=200) :: message
    LOGICAL :: LEAP

    integer iphase, itype, p1st,l_bc,ko,n,l
    real dens_bc,bc_mass_per_sec,so2_mass_per_sec,mass1part

    INTEGER :: dbg_lvl = 1

    real, save :: size_dist(8) = (/ 0.0494,  0.3795,  0.4714,  0.0967,  0.003,  0.0,  0.0, 0.0 /)   ! fractions

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Setup the timing flags...
      ! (methodology is adapated from module_radiation_driver.F)

      CALL wrf_debug ( dbg_lvl , 'fire_driver_kuwait_phys: adding heat flux start' ) 

      do ij=1,num_tiles
        ! set up tile bounds    
        its = i_start(ij)  ! start atmospheric tile in i
        ite = i_end(ij)    ! end atmospheric tile in i
        jts = j_start(ij)  ! start atmospheric tile in j
        jte = j_end(ij)    ! end atmospheric tile in j

        do j=jts,jte
        do i=its,ite
           !ifire = 0

           if(fire_hflx(i,j).le.1) cycle

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! CALCULATING A TIME OF ERUPTION
            LEAP = .FALSE.
            IF((MOD(fire_year(i,j),4.0) .EQ. 0 .AND. MOD(fire_year(i,j),100.0).NE.0 ).OR. MOD(fire_year(i,j),400.0).EQ.0 ) THEN
                      LEAP = .TRUE.        
                  ENDIF
            IF (LEAP) THEN
              K = 1
            ELSE
                K = 2
            END IF
            beg_jul = ((275*fire_month(i,j))/9) - K*((fire_month(i,j)+9)/12) + fire_day(i,j) - 30

            ! The code aboive is correct, but WRF shows julday-1  
            beg_jul=beg_jul-1

            end_hour=fire_hour(i,j)+int(fire_dur(i,j)/3600.)
            end_hour=end_hour-24*int(end_hour/24)
            
            ! ukhova: 28 July 2019
            ! added 1 for Raikoke
            !endday=beg_jul+int(fire_dur(i,j)/(3600*24))+1

            ! FOR Pinatubo
            endday=beg_jul+int(fire_dur(i,j)/(3600*24))

            write(0,*)' julday = ',julday
            write(0,*)' beg_jul, fire_hour,beg_jul+fire_hour(i,j)/24.0 = ',beg_jul, fire_hour(i,j),beg_jul+fire_hour(i,j)/24.0
            write(0,*)' end_jul, end_hour,endday+end_hour/24.0 = ',endday,end_hour,endday+end_hour/24.0
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(julday.ge.beg_jul+fire_hour(i,j)/24.0 .and. julday.le.endday+end_hour/24.0)then
              CALL wrf_debug ( dbg_lvl , ' ***FIRE!!!!')

              ! heat flux W/m2
              grnhfx(i,j)=fire_hflx(i,j)

              ! emissions of BC
              ! for now just do itype=1
              itype   = 1
              iphase  = ai_phase
              
              dens_bc = 1.8 ! g/cm3
              ko      = 1

              !fire_bc_emis - total emissions of BC in kg during period of fire fire_dur (seconds)
              bc_mass_per_sec=1.0e09*fire_bc_emis(i,j)/fire_dur(i,j)    ! ug per second per cell

              so2_mass_per_sec=fire_so2_emis(i,j)/fire_dur(i,j)         ! kg per second per cell

              area(1,1)=(dx/mapfac_mx(i,j))*(dy/mapfac_my(i,j))         ! surface area in m2
              airmas(1,1)=rho(i,ko,j)*area(1,1)*dz8w(i,ko,j)            ! air mass inside cell (kg)

              do 1810 n = 1, nsize_aer(itype)
                  l_bc = lptr_bc_aer(n,itype,iphase)
                  if (l_bc >= p1st) then
                    chem(i,ko,j,l_bc) = chem(i,ko,j,l_bc)+(bc_mass_per_sec*size_dist(n)*dt)/airmas(1,1)

                    mass1part=0.523598*(dcen_sect(n,itype)**3)*dens_bc*1.0e06 !mass (ug) of a BC particle with D=dcen_sect (cm)

                    ! #/kg
                    l = numptr_aer(n,itype,iphase)
                    if (l >= p1st) chem(i,ko,j,l) = chem(i,ko,j,l) + (bc_mass_per_sec*size_dist(n)*dt)/(mass1part*airmas(1,1))
                  end if
              1810    continue

              !so2 emission
              chem(i,ko,j,p_so2) = chem(i,ko,j,p_so2) + so2_mass_per_sec*453125*dt/airmas(1,1)

            else
              grnhfx(i,j)=0.0
              CALL wrf_debug ( dbg_lvl , ' ***NO FIRE!!!!')
            endif

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
        enddo !its
        enddo !jts

      enddo  !num_tiles      
      CALL wrf_debug ( dbg_lvl , 'fire_driver_kuwait_phys: adding heat flux stop' ) 

end subroutine fire_driver_kuwait_phys


! --- this routine takes fire generated heat flux and
!     calculates their influence on the theta
! --- note that this tendency is valid at the Arakawa-A location
SUBROUTINE kuwait_fire_tendency( &
    ids,ide, kds,kde, jds,jde,   & ! dimensions
    ims,ime, kms,kme, jms,jme,   &
    its,ite, kts,kte, jts,jte,   &
    grnhfx,                      & ! heat fluxe
    zs,z_at_w,dz8w,mu,c1h,c2h,rho, &
    rthfrten)             ! theta and Qv tendencies

   use module_model_constants, only: cp

   IMPLICIT NONE

! --- incoming variables

   INTEGER , INTENT(in) :: ids,ide, kds,kde, jds,jde, &
                           ims,ime, kms,kme, jms,jme, &
                           its,ite, kts,kte, jts,jte

   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: grnhfx  ! flux from the surface (W/m^2)
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: zs  ! topography (m abv sealvl)
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: mu  ! dry air mass (Pa)
   REAL, INTENT(in), DIMENSION( kms:kme         ) :: c1h, c2h ! Hybrid coordinate weights

   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: z_at_w ! m abv sealvl
   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: dz8w   ! dz across w-lvl
   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: rho    ! density

! --- outgoing variables

   REAL, INTENT(out), DIMENSION( ims:ime,kms:kme,jms:jme ) :: rthfrten  ! theta tendency from fire (in mass units)

! --- local variables
   INTEGER :: i,j,k
   INTEGER :: i_st,i_en, j_st,j_en, k_st,k_en

   REAL :: cp_i,alfg_i,rho_i
   REAL :: z_w

   REAL :: alfg ! extinction depth of sfc fire heat"   "m"

   REAL, DIMENSION( its:ite,kts:kte,jts:jte ) :: hfx
   
    do j=jts,jte
        do k=kts,min(kte+1,kde)
           do i=its,ite
               rthfrten(i,k,j)=0.
           enddo
        enddo
    enddo


! --- set some local constants
   alfg = 50        ! meters
   cp_i = 1./cp     ! inverse of specific heat
   alfg_i = 1./alfg

! --- set loop indicies : note that

   i_st = MAX(its,ids+1)
   i_en = MIN(ite,ide-1)
   k_st = kts
   k_en = MIN(kte,kde-1)
   j_st = MAX(jts,jds+1)
   j_en = MIN(jte,jde-1)

! --- distribute flux over vertical

   DO j = j_st,j_en
      DO k = k_st,k_en
         DO i = i_st,i_en
            ! --- set z (in meters above ground)
            z_w = z_at_w(i,k,j) - zs(i,j) ! should be zero when k=k_st
            ! --- heat flux
            hfx(i,k,j) =  grnhfx(i,j) * cp_i * EXP( - alfg_i * z_w )
         END DO
      END DO
   END DO

! --- add flux divergence to tendencies
!
!   multiply by dry air mass (mu) to eliminate the need to
!   call sr. calculate_phy_tend (in dyn_em/module_em.F)

!  calculate_phy_tend couples the physics tendencies to the column mass (mu),
!  because prognostic equations are in flux form, but physics tendencies are
!  computed for uncoupled variables.
!
   DO j = j_st,j_en
      DO k = k_st,k_en-1
         DO i = i_st,i_en
            rho_i = 1./rho(i,k,j)
            rthfrten(i,k,j) = - (c1h(k)*mu(i,j)+c2h(k)) * rho_i * (hfx(i,k+1,j)-hfx(i,k,j)) / dz8w(i,k,j)
         END DO
      END DO
   END DO

   RETURN

END SUBROUTINE kuwait_fire_tendency


end module module_kuwait_fire_driver


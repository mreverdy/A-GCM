!****************************************************************************!
!                                                                            !
!                            GOECP 1.00 			             !
!                            						     !
!									     !
!                    Mathieu Reverdy / LMD / ESA                             !
!                      last update **/**/****                                !
!                                                                            !
!                                                                            !
!  1) Purpose : compute instantaneous and daily mean profiles of Scattering  !
!               Ration (SR), Color Ratio (CR) & Depolarization Ratio (DR)    !
!               over a model grid.                                           !
!                                                                            !
!               Compute various daily/monthly cloudiness over a model grid : !
!                  - Map of Low Mid High cloud cover                         !
!                  - 3D Cloud Fraction                                       !
!                  - 3D Cloud Phase                                          !
!                  - SR Histograms                                           !
!                                                                            !
!  2) Input   : SDS & META variables from CALIPSO hdf level1 Data files.     !
!                                                                            !
!  3) Output  : - 3D_CloudFraction :                                         !
!                    clcalipso(lon,lat,alt,time)                             !
!                    clrcalipso(lon,lat,alt,time)                            !
!                    uncalipso(lon,lat,alt,time)                             !
!                                                                            !
!               - Map of Low Mid High :                                      !
!                    cllcalipso(lon,lat,time)                                !
!                    clmcalipso(lon,lat,time)                                !
!                    clhcalipso(lon,lat,time)                                !
!                    cltcalipso(lon,lat,time)                                !
!                    clccalipso(lon,lat,time)                                !
!                                                                            !
!               - SR_Histo :                                                 !
!                    cfad_lidarsr532_Occ(lon,lat,alt,box,time)               !
!                    cfad_lidarsr532_Occ2(lon,lat,alt,box2,time)             !
!                                                                            !
!               - 3D_CloudPhase :                                            !
!                    ice_cloud(lon,lat,alt,time)                             !
!                    water_cloud(lon,lat,alt,time)                           !
!                                                                            !
!               - instant_SR_CR_DR                                           !
!                    longitude(it)                                           !
!                    latitude(it)                                            !
!                    altitude(it)                                            !
!                    time(it)                                                !
!                    SE(it)                                                  !
!                    instant_SR(it,alt)                                      !
!                    instant_CR(it,alt)                                      !
!                    instant_DR(it,alt)                                      !
!                                                                            !
!  4) Grid :                                                                 !
!               -CFMIP2   : 2° x 2° x 40levels from 0 to 19.2km              !
!                           (180,90,40)                                      !
!               -CFMIP1   : 1° x 1° x 40levels from 0 to 19.2km              !
!                           (360,180,40)                                     !
!               -CFMIP2.5 : 2.5° x 2.5° x 40levels from 0 to 19.2km          !
!                           (144,72,40)                                      !
!               -CFMIP    : 3.75° x 2.5° x 40levels from 0 to 19.2km         !
!                           (96,72,40)                                       !
!               -LMDZ     : 3.75° x 2.5° x 19levels from 0 to 40.4823km      !
!                           (96,72,19)                                       !
!               -LMDZ40   : 3.75° x 2.5° x 40levels from 0 to 40.4823km      !
!                           (96,72,40)                                       !
!               -NASA     : 5° x 5° x 40levels each 480m from 0 to 19.2km    !
!                           (73,37,41)                                       !
!                                                                            !
!  5) Compilation : use the makefile "makefiles.sh"                          !
!      makefiles.sh $1                                                       !
!      ifort $1.f90 -I/usr/include/hdf -L/usr/lib64/hdf -lmfhdf -ldf -ljpeg  !
!                   -lz -I/opt/netcdf-3.6.0-p1-ifort-64/include              !
!                   -L/opt/netcdf-3.6.0-p1-ifort-64/lib  -lnetcdf -o $1.e    !
!                                                                            !
!  6) Last updates & bug fix :                                               !
!   - 15/11/13 Draft 		                                             !
!									     !
!                                                                            !
!****************************************************************************!

!************************* SUBROUTINE SCHEME ********************************!
!                                                                            !
! sdsread8            : x2 for time                                          !
! sdsread             : x8 for atb,atb1064,mol,perp,temp,lat,lon,SE          !
! metaread            : x2 for altl & altm                                   !
!                                                                            !
!----------------------------------------------------------------------------!
!                                                                            !
! interp              : x4 for Pressure and molecular                        !
! atb_mol             : x4 for normalized ratio & molecular calculation      !
! atb_mol_interp      : x2 for molecular extrapolation                       !
! SE_alt_chim         : x4 for adding the SE to atb & mol                    !
!                                                                            !
!----------------------------------------------------------------------------!
!                                                                            !
! vertical_mean       : x6 for the average of atb, mol, perp, temp           !
! zero_detect         : x8 for the detection of empty boxes                  !
! Surf_detect2        : x8 for adding the Surface Elevation                  !
! SE_km_2_pres2       : x2 for adding the Surface Elevation                  !
! SR_CR_DEPOL_mean    : x16 for the SR CR & DR calculation                   !
! filtre_2lvl         : x4 for the SR calculation with a filter up to 21km   !
!                          (over 21kilometers grid such as LMDZ/LMDZ40 grid) !
!                                                                            !
!----------------------------------------------------------------------------!
!                                                                            !
! SR_CR_DR_2nc        : x1 for the recording of instant_SR instant_CR and    !
!                       instant_DR variables in netcdf files                 !
!                                                                            !
!----------------------------------------------------------------------------!
!                                                                            !
! fraction_subgrid2_8km   : x1 for cloud,clear,uncertain flags calculation   !
!                                                                            !
!----------------------------------------------------------------------------!
!                                                                            !
! create_mapnc        : x1 for the creation of the MapLowMidHigh file        !
! create_mapnc2       : x1 for the creation of the GEWEX MapLowMidHigh file  !
! map_recvar2nc2      : x1 for the recording of map variables                !
! map_recvar2nc7      : x1 for the recording of GEWEX map variables          !
! create_diagnc       : x1 for the creation of the SR_histograms file        !
! diag_recvar2nc3     : x1 for the recording of the SR_histograms variables  !
! create_profnc       : x1 for the creation of the 3D_CloudFraction files    !
! prof_recvar2nc      : x1 for the recording of the 3D_CloudFraction var     !
! create_depolnc3D    : x1 for the creation of the 3D_CloudPhase files       !
! depol_recvar2nc     : x1 for the recording of the 3D_CloudPhase var        !
!                                                                            !
!----------------------------------------------------------------------------!

program A-GCM

  use netcdf
  implicit none


!****************************************************************************!
!*!!!!!!!!!!!! DEFINITIONS & DECLARATIONS OF VARIABLES !!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


!----------------------------------------------------------------------------!
! Date and time variable                                                     !
!----------------------------------------------------------------------------!
!               variables used in subroutine DATE_AND_TIME()                 !
!****************************************************************************!
!   value     : 8 dimension value which contains the year,month,day,hour,    !
!               minute, seconds and milliseconds of the real-time.           !
!   datenc    : date from the real-system clock, format: yyyymmdd.           !
!   timenc    : time from the real-system clock, format: hhmmss.sss.         !
!   zone      : represente the difference with respect to Coordinated        !
!               Universal Time (UTC), format: (+-)hhmm.                      !
!                                                                            !
!----------------------------------------------------------------------------!
!               variables used all along the program                         !
!****************************************************************************!
!   date      : first value of the variable time in the CALIPSO hdf file.    !
!   year      : year of the CALIPSO hdf file (derived from date var).        !
!   month     : month of the CALIPSO hdf file (derived from date var).       !
!   day       : day of the CALIPSO hd file (derived from date & month var).  !
!   jour      : day of the CALIPSO hd file (derived from date & month var).  !
!   resd      : number of days since 2000/01/01 for the trimonth period.     !
!   resh      : number of hour since 2000/01/01 for the trimonth period.     !
!                                                                            !
!----------------------------------------------------------------------------!
!               variables used in the function interdat()                    !
!****************************************************************************!
!   interdat  : routine which allow to perform the calculation of the number !
!               of days since the 2000/01/01 to the middle of the period     !
!               looked by the routine.                                       !
!   date1     : input of the routine interdat. It is the date of the midle   !
!               on the period looked by the routine, format: yyyymmddhh      !
!   date2     : input of the routine interdat. It is the date of the         !
!               first day since the year 2000, format: yyyymmddhh.           !
!                                                                            !
!----------------------------------------------------------------------------!
      integer,dimension(8)  ::  value
      character  ::  datenc*8, timenc*10, zone*5 
      character  ::  blanks*15,forme*15
      character(len=15),parameter ::  version="Num_version"
      integer  ::  date, year, month, day, jour          
      integer  ::  date1, date2, resh                  
      integer  ::  interdat, ret1
      real  ::  resd, test    

!----------------------------------------------------------------------------!
! Calculation of time elpased                                                !
!----------------------------------------------------------------------------!
      integer  ::  t1, t2, ir
      real     ::  t_cpu, t_cpu_0, t_cpu_1, tempstot
 
!----------------------------------------------------------------------------!
! Name of keyboarding data and some char variables                           !
!----------------------------------------------------------------------------!
!                            input files                                     !
!****************************************************************************!
!   file2     : path and name of hdf CALIPSO input file                      !
!   file3     : path and name of the file which lists all the hdf CALIPSO    !
!               file name and path (=file2)                                  !
!                                                                            !
!----------------------------------------------------------------------------!
!                            output files                                    !
!****************************************************************************!
!   file4     : instant_SR_CR_DR output file name                            !
!   file5     : 3D_CloudFraction output file name                            !
!   file6     : MapLowMidHigh output file name                               !
!   file7     : SR_histo output file name                                    !
!   file10    : 3D_CloudFraction_Phase output file name                               !
!   file11    : MapLowMidHigh_Phase output file name                               !
!                                                                            !
!----------------------------------------------------------------------------!
!                             model grid                                     !
!****************************************************************************!
!   model     : name of model selected (at that time only lmdz is available) !
!   gcm       : name of the grid selected description in section 4)          !
!                                                                            !
!----------------------------------------------------------------------------!
!                           switch variables                                 !
!****************************************************************************!
!   sauve     : = "chim" "wrf" or "off" to record output files in ASCII      !
!   switch    : = "night" or "day" to select day or night mode               !
!   switch2   : = "sat" or "cloudy" to select sat or cloudy mode             !
!   alt_pres  : = "altitude" or "pressure" to select altitude or pressure    !
!               mode.                                                        !
!                                                                            !
!----------------------------------------------------------------------------!
!                          other char variables                              !
!****************************************************************************!
!   metal     : name of lidar altitude variable metadata over 583lvl         !
!   metam     : name of meteo altitude variable metadata over 33lvl          !
!   sds_varname : name of SDS variables from hdf CALIPSO files               !
!   datec     : date in characters converted from date                       !
!   numfichc  : number of hdf files read in characters, converted from       !
!               numfich                                                      !
!                                                                            !
!----------------------------------------------------------------------------!
      character  ::  file2*1024,  file3*132, file4*132 , file5*132   
      character  ::  file6*132, file7*132, file8*132, file9*132
      character  ::  file10*132, file11*132, file12*132,file13*132 
      character  ::  file66*132   
      character  ::  model*30, filetmp*1024, filetmp2*1024,instantname*1024
      character  ::  metal*30,metam*30
      character  ::  command*1024, command2*1024, command3*1024,command4*1024
      character  ::  sauve*3,gcm*8, idiagc*2, idepc*2
      character  ::  sds_varname*100
      character  ::  switch*5, alt_pres*8, switch2*6,instant_switch*8
      character  ::  datec*6, numfichc*4, yearc*4
      character  ::  grid*8, altfile*10, lonfile*18, latfile*17
  100 format(A200)
                 
                                       
!----------------------------------------------------------------------------!
! id used in netcdf output files with subroutine create_*nc()                !
!----------------------------------------------------------------------------!
!   dimidsp   : dimension id of the MaP3D variables recorded in the ncdf     !
!               files                                                        !
!   dimidsm   : dimension id of the MaPLowMidHigh variables recorded in the  !
!               ncdf files                                                   !
!   dimidsd   : dimension id of diagSR variable recorded in the ncdf file    !
!   nan       : nan is the NaN value in the ncdf files                       !
!                                                                            !
!----------------------------------------------------------------------------!
      real,parameter  ::   nan=-9999.
      real*4,parameter  ::  SeuilMol1km = 0.00015 , SeuilTemp1km = 6.5
      real*4  ::  seuilatb, deltatb  
      integer  ::  dimidsp(4),dimidsp2(5), dimidsm(3), dimidsm2(4),dimidphase(4),dimidsd(5),dimidsdb(5)
      integer  ::  dimidhist3(4), dimidsd2(4),dimidsd3(5),dimidsd4(6),dimidpha(6)
      integer  ::  dimidhist(4),dimidhist2(4)

!----------------------------------------------------------------------------!
! Output variable dimensions                                                 !
!----------------------------------------------------------------------------!
!   latmax    : number of the latitude boxes, according to the model chosen  !
!               CFMIP=181, LMDZ=72, NASA=37                                  !
!   lonmax    : number of the longitude boxes, CFMIP=361, LMDZ=96, NASA=73   !
!   altmax    : number of the altitude boxes, CFMIP=41, LMDZ=19, NASA=41     !
!   diagmax   : number of the diagSR boxes, always equal to 15               !
!   daymax    : maximum number of day in a month (for the daily var)         !
!                                                                            !
!----------------------------------------------------------------------------!
      integer  ::  latmax, lonmax, altmax, altmax2
      integer, parameter  ::  diagmax = 19 , diagmax2 = 11, depolmax = 21
      integer, parameter  ::  daymax = 31, pr2max=301, permax=241, tempmax=39
      integer  ::  toplowl,topmidl,tophighl,toplvlsat1,toplvlsat2
      integer  ::  altstart,altend,nol
                  
                     
!----------------------------------------------------------------------------!
! Loop variables                                                             !
!----------------------------------------------------------------------------!
!                        obs & model variables                               !
!****************************************************************************!
!   i         : loop index of the number of profil in the hdf file read      !
!   iz        : loop index of the altitude boxes                             !
!   altitude  : number of level of the lidar variables                       !
!                                                                            !
!----------------------------------------------------------------------------!
!                            obs variables                                   !
!****************************************************************************!
!   ilid      : loop index of the lidar altitude (583)                       !
!   altitude2 : number of level of the meteo variables                       !
!                                                                            !
!----------------------------------------------------------------------------!
!                        CHIMERE & WRF variables                             !
!****************************************************************************!
!   j         : loop index of the latitude boxes                             !
!   k         : loop index of the longitude boxes                            !
!                                                                            !
!----------------------------------------------------------------------------!
!                    model variables (daily/monthly)                         !
!****************************************************************************!
!   ilon      : loop index of the longitude boxes in the daily/monthly var   !
!   ilat      : loop index of the latitude boxes in the daily/monthly var    !
!   idiag     : loop index of the diagSR boxes in the daily/monthly var      !
!                                                                            !
!----------------------------------------------------------------------------!
!                         Unspecific variables                               !
!****************************************************************************!
!   numfich   : number of hdf file read, = number of line in file3           !
!   comptpf   : number of boxes crossed by the satellite                     !
!   box       : number total of boxes in the WRF or CHIMERE grid             !
!                                                                            !
!----------------------------------------------------------------------------!
      integer  ::  iz, i, j, k, ilat, ilon, ialt, ilid, idiag, idep, ipr2,iperp,itemp  
      integer(kind=2),parameter :: altitude2 = 33  
      integer(kind=2),parameter :: altitude = 583    
      integer  ::  numfich, comptpf,  box , nphase       

   
!----------------------------------------------------------------------------!
! Status variables                                                           !
!----------------------------------------------------------------------------!
!   err       : return 0 if the file could be opened                         !
!   OK_buffer : return 0 if the allocation is OK                             !
!                                                                            !
!----------------------------------------------------------------------------!
      integer  ::  err , OK_buffer                     
                     

!----------------------------------------------------------------------------!
! Input variable dimensions                                                  !
!----------------------------------------------------------------------------!
!   it        : number of profil (=nprof) in the hdf file read (about 60000) !
!                                                                            !
!----------------------------------------------------------------------------!
      integer :: it
      integer*4,dimension(:,:),allocatable :: indretmean, indtotmean, indtot
      integer*4,dimension(:,:),allocatable :: indretlowmean,indretmidmean
      integer*4,dimension(:,:),allocatable :: indrethighmean 
      integer*4,dimension(:),allocatable :: indret, indretlow, indretmid
      integer*4,dimension(:),allocatable :: indrethigh

                     
!----------------------------------------------------------------------------!
! SDS variables from CALIPSO hdf files, extracted by routine sdsread (1bit) &!
! sdsread8 (8bit)                                                            !
!----------------------------------------------------------------------------!
!   atb       : Total_Attenuated_Backscatter_532 per meter per steradian !
!               has dimension (altitude,nprof)                               !
!   mol       : Molecular_Number_Density in count numerical CN has dimension !
!               (altitude2,nprof)                                            !
!   pres      : Pressure in hPa has dimension (altitude2)                    !
!   temps     : Profil_UTC_time has format yymmdd.ffffffff, has dimension    !
!               (nprof)                                                      !
!   vartmp8   : Temporal var for the time reading, has dimension (nprof)     !
!   vartmp    : Temporal var for the lat,lon,SE reading, has dimension       !
!               (nprof)                                                      !
!   temps2    : Time converted in UTC fractionned hour of a day (ex : 1.5 =  !
!               1h30 am), has dimension (nprof)                              !
!   lat       : Latitude in degrees, has dimension (nprof)                   !
!   lon       : Longitude in degrees, has dimension (nprof)                  !
!   SE        : Surface_Elevation in kilometers, has dimension (nprof)       !
!                                                                            !
!----------------------------------------------------------------------------!
      real,dimension(:,:),allocatable  ::  atb, mol,  pres, atb2, perp, temp
      real*8,dimension(:),allocatable  ::  temps, temps2,temps2wrf
      real*8,dimension(:,:),allocatable  :: vartmp8
      real,dimension(:,:),allocatable  :: vartmp
      real,dimension(:),allocatable  ::  lat, lon
      real,dimension(:),allocatable  ::  SE,SEwrf
      real,dimension(:),allocatable  ::  latwrf,lonwrf
                    
!----------------------------------------------------------------------------!
! SDS variables interpolated verticaly from CALIPSO variables                !
!----------------------------------------------------------------------------!
!   mol2      : Molecular interpolated from 33lvl to 583lvl in CN            !
!               dim=(altitude,nprof)                                         !
!   pres2     : Pressure interpolated from 33lvl to 583lvl in hPa            !
!               dim=(altitude,nprof)                                         !
!   mol3      : Molecular converted from mol2(CN) in km-1 sr-1               !
!               dim=(altitude,nprof)                                         !
!   SEp       : Surface_Elevation in hPa derived from SE(km), dim=(nprof)    !
!                                                                            !
!----------------------------------------------------------------------------!
      real,dimension(:,:),allocatable  ::  mol2, pres2, mol3, mol4, temp2 
      real,dimension(:),allocatable  ::  SEp
                  
                   
!----------------------------------------------------------------------------!
! Output & Grid variables                                                    !
!----------------------------------------------------------------------------!
!   latmod    : values of model latitude from 90 to -90                      !
!               dimension : CFMIP=181, LMDZ=72, NASA=37                      !
!   lonmod    : values of model longitude from 180 to -180                   !
!               dimension : CFMIP=361, LMDZ=96, NASA=73                      !
!   altmod    : values of model altitude from 0 to 19.2km (CFMIP) or 40.5km  !
!               (LMDZ)                                                       !
!               dimension : CFMIP=41, LMDZ=19, NASA=41                       !
!               model, has dimension 41 or 19                                !
!   prestop   : same var as altmod in pressure mod                           !
!   srmod     : values of diagSR, dimension always equal to 15               !
!               values : -1,0,0.01,1.2,2,3,5,10,20,30,40,50,60,80,100        !
!   pr2moy    : observed atb averaged on altmod lvl with the routine         !
!               vertical_mean,following the satellite,                       !
!               dim=(altmax,nprof)                                           !
!   molmoy    : mol3 averaged on altmod lvl with the routine vertical_mean,  !
!               following the satellite, dim=(altmax,nprof)                  !
!   srmoy     : scattering ratio calculated with pr2moy and molmoy,          !
!               dim=(altmax,nprof)                                           !
!   indice    : number of atb iteration in an altitude box                   !
!   indicem   : number of mol iteration in an altitude box                   !
!   indiceh   : number of time iteration in a lat/lon/alt box with CHIM/WRF  !
!               model                                                        !
!   mheure    : time averaged on lon/lat/alt grid with CHIM/WRF model        !
!                                                                            !
!----------------------------------------------------------------------------!
      real*4,dimension(:),allocatable  ::  latmod, lonmod, srmod, prestop, depolmod, pr2mod, srdepmod,atbrmod, tempmod!, crmod
      real*4,dimension(:),allocatable  ::  latmid, lonmid      
      real,dimension(:),allocatable  ::  altmod, altmid,tempmid
      real,dimension(:,:),allocatable  :: altmod_bound, tempmod_bound
      real,dimension(:,:),allocatable  ::  pr2moy,  molmoy, srmoy,depolmoy, pr2moy2,crmoy,perpmoy,parmoy, tempmoy
      real,dimension(:,:),allocatable  ::  indice, indicem, indiceh, indice2,indicep,indicep2,indicetemp 
      real*4,dimension(:,:),allocatable  ::  mheure
      real*4,dimension(:,:),allocatable  :: SRwrf,CRwrf,DEPOLwrf            

!----------------------------------------------------------------------------!
! LMDZ output variables                                                      !
!----------------------------------------------------------------------------!
!               instantaneous fraction : dim=(altmax,nprof)                  !
!****************************************************************************!
!   uncertfraction   : uncertain fraction flag 0/1 for each profil and each  !
!                      altitude, calculated with routine fraction_subgrid2   !
!   satfraction      : fully attenuated flag 0/1                             !
!   cloudfraction    : cloudy fraction flag 0/1                              !
!   clearfraction    : clear fraction flag 0/1                               !
!   nanfraction      : NaN fraction flag 0/1                                 !
!   sefraction       : Surface_Elevation fraction flag 0/1                   !
!   fractot          : all 6 fraction unified in a single var, by flag 1,2,3,!
!                      4,5,6 used in the routine fraction_subgrid3 :         !
!                      1=sat, 2=clear, 3=uncert, 4=nan, 5=SE, 6=cloud        !
!                                                                            !
!----------------------------------------------------------------------------!
!               daily fraction : dim=(lonmax,latmax,altmax,day)              !
!****************************************************************************!
!   uncertfractday   : uncertfraction averaged on LMDZ/CFMIP/NASA grid day by!
!                      day.                                                  !
!   satfractday      : satfraction averaged on LMDZ/CFMIP/NASA grid          !
!   cloudfractday    : cloudfraction averaged on LMDZ/CFMIP/NASA grid        !
!   clearfractday    : clearfraction averaged on LMDZ/CFMIP/NASA grid        !
!   nanfractday      : nanfraction averaged on LMDZ/CFMIP/NASA grid          !
!   sefractday       : sefraction averaged on LMDZ/CFMIP/NASA grid           !
!   indday           : number of sat/cloud/clear/uncert values on            !
!                      lon/lat/alt/day boxes                                 !
!   inddaytot        : number of sat/cloud/clear/uncert/nan/se values on     !
!                      lon/lat/alt/day boxes                                 !
!                                                                            !
!----------------------------------------------------------------------------!
!              monthly fraction : dim=(lonmax,latmax,altmax)                 !
!****************************************************************************!
!   monthuncertfract : uncertfractday averaged on all days of the month      !
!   monthsatfract    : satfractday averaged on all days of the month         !
!   monthcloudfract  : cloudfractday on all days of the month                !
!   monthclearfract  : clearfractday on all days of the month                !
!   monthnanfract    : nanfractday averaged on all days of the month         !
!   monthsefract     : sefractday averaged on all days of the month          !
!   indpermonth      : number of sat/cloud/clear/uncert fractday on          !
!                      lon/lat/alt boxes                                     !
!   indpermonthtot   : number of sat/cloud/clear/uncert/nan/se fractday on   !
!                      lon/lat/alt boxes                                     !
!                                                                            !
!----------------------------------------------------------------------------!
!                    instantaneous isccp : dim=(nprof)                       !
!****************************************************************************!
!   isccplow        : isccp low cloud flag 0/1 for each profil               !
!   isccpmid        : isccp mid cloud flag 0/1 for each profil               !
!   isccphigh       : isccp high cloud flag 0/1 for each profil              !
!   colcloud        : isccp column cloud flag 0/1 for each profil            !
!   colclear        : isccp column clear flag 0/1 for each profil            !
!                                                                            !
!----------------------------------------------------------------------------!
!                daily isccp : dim=(lonmax,latmax,daymax)                    !
!****************************************************************************!
!   isccplowday     : isccplow averaged on LMDZ/CFMIP/NASA grid              !
!   isccpmidday     : isccpmid averaged on LMDZ/CFMIP/NASA grid              !
!   isccphighday    : isccphigh averaged on LMDZ/CFMIP/NASA grid             !
!   colcloudday     : colcloud averaged on LMDZ/CFMIP/NASA grid              !
!   colclearday     : colclear averaged on LMDZ/CFMIP/NASA grid              !
!   isccpindday     : number of isccp low/mid/high on a lon/lat/day box      !
!                                                                            !
!----------------------------------------------------------------------------!
!                      monthly isccp : dim=(lonmax,latmax)                   !
!****************************************************************************!
!   monthisccplow   : isccplowday averaged on all days of the month          !
!   monthisccpmid   : isccpmidday averaged on all days of the month          !
!   monthisccphigh  : isccphighday averaged on all days of the month         !
!   monthcolcloud   : colcloudday averaged on all days of the month          !
!   monthcolclear   : colclearday averaged on all days of the month          !
!   isccpdaypermonth: number of monthisccp low/mid/high on a lon/lat box     !
!                                                                            !
!----------------------------------------------------------------------------!
!                            diagSR variables                                !
!****************************************************************************!
!   diagSR          : number of occurence of a SR value in one of the 15     !
!                     different interval defined with srmod                  !
!                     dim=(lonmax,latmax,altmax,daymax,diagmax)              !
!   monthdiagSR     : diagSR averaged on all days of the month               !
!                     dim=(lonmax,latmax,altmax,diagmax)                     !
!   sumdiag         : sum of diagSR at an altitude                           !
!                                                                            !
!----------------------------------------------------------------------------!
      real,dimension(:,:),allocatable :: uncertfraction,satfraction,       &
                                           cloudfraction,clearfraction,      &
                                           nanfraction,sefraction,           &
                                           rejfraction,fractot,cloudfraction2

      real,dimension(:,:,:,:),allocatable :: cloudfractday, clearfractday, &
                                           !    satfractday, sefractday ,     &
                                               uncertfractday, &! nanfractday,  &
                                               indday!,inddaytot  

      real,dimension(:,:,:),allocatable :: monthcloudfract,monthclearfract,&
                                             monthuncertfract, &!monthnanfract, &
                                            ! monthsatfract, monthsefract,    &
                                             indpermonth,indphasepermonth

      integer,dimension(:,:,:),allocatable  ::  indnan

      integer,dimension(:),allocatable :: isccplow, isccpmid, isccphigh,      &
                                         colcloud!, colclear  
      integer,dimension(:,:),allocatable :: isccpliq, isccpice
      integer,dimension(:,:,:),allocatable :: isccpun

 
       integer  ::  indbase
       real,dimension(:,:,:),allocatable  ::  heightday2
       real,dimension(:),allocatable  ::  height2
       real,dimension(:,:),allocatable  ::  monthheight2

       real,dimension(:,:,:),allocatable  ::  heightday,indheight
       real,dimension(:,:),allocatable  ::  monthheight,indmonthheight
       real,dimension(:),allocatable  ::  heightmod
       integer,dimension(:,:,:),allocatable  ::  hlow,hmid,hhigh,hcol
       real,dimension(:,:,:),allocatable  ::  hheight
       real,dimension(:),allocatable  ::  height
       integer  ::  iheight, ihist, ihisttemp, icat
       integer,parameter  ::  heightmax=41, catmax=5
       integer,parameter  ::  histmax=11, histmax2=41, histtempmax=29,histtempmax2=28
       integer,dimension(histmax2)  ::  histmod2
       real,dimension(histmax)  ::  histmod
 
       real,dimension(histtempmax)  ::  histtempmod
       real,dimension(histtempmax2)  ::  histtempmod2
 
       integer,dimension(:,:,:),allocatable  :: hcoltemp,hlowtemp,hmidtemp,hhightemp
       real,dimension(:,:,:),allocatable  :: coltemp,lowtemp,midtemp,hightemp
       real,dimension(:,:,:),allocatable  :: indcoltemp,indlowtemp,indmidtemp,indhightemp
 
       real,dimension(:,:),allocatable  ::  monthlowtemp,monthmidtemp,monthhightemp,monthcoltemp
       real,dimension(:,:),allocatable  ::  indmonthlowtemp,indmonthmidtemp,indmonthhightemp,indmonthcoltemp


      real,dimension(:,:,:),allocatable :: isccplowday, isccpmidday,       &
                                             isccphighday, isccpindday,      &
                                             isccpinddaylow, isccpinddaymid, &
                                             colcloudday





      real*8  ::  colclearres

      real  ::  isccptemp

      real,dimension(:,:,:,:),allocatable :: isccpiceday, isccpliqday
      real,dimension(:,:,:,:,:),allocatable :: isccpunday
      real,dimension(:,:,:,:),allocatable :: isccpphaseday
      real,dimension(:,:,:),allocatable :: isccpdustday

      real,dimension(:,:,:),allocatable :: monthisccpice, monthisccpliq
      real,dimension(:,:,:,:),allocatable :: monthisccpun
      real,dimension(:,:,:),allocatable :: monthisccpphase,indmonthphase
      real,dimension(:,:,:),allocatable :: indmonthphase2
      integer,dimension(:,:,:),allocatable  ::  indmonthphase3D

      real,dimension(:,:),allocatable :: monthisccplow, monthisccpmid,     &
                                           monthisccphigh, isccpdaypermonth, &
                                           isccpdaypermonthlow,              &
                                           isccpdaypermonthmid,              &
                                           monthcolcloud, monthcolclear

     
      real,dimension(:,:,:,:),allocatable :: diagSR
      real,dimension(:,:,:,:,:),allocatable :: diagSRpha

      real,dimension(:,:,:,:,:),allocatable :: diagPHA


      integer*4  ::  sumdiag


!----------------------------------------------------------------------------!
! META variables                                                             !
!----------------------------------------------------------------------------!
!   altl   : altitude of lidar lvl in kilometer, dim=(altitude)              !
!   altm   : altitude of meteo lvl in kilometer, dim=(altitude2)             !
!                                                                            !
!----------------------------------------------------------------------------!
      real*4,dimension(:),allocatable  ::  altl, altm  
      
      integer  ::  seuilsnrlow, seuilsnrhigh
      integer  ::  nanprof, nansat, nanmid, nanlow
      integer  ::  icewaterres
      real,parameter  ::  alpha=0.0028, beta=0.0458, alpha1=3., beta1=0.0576

      real,parameter  ::  alpha50=9.0322e+03, beta50=-2.1358e+03, gamma50=173.3963, delta50=-3.9514, epsilon50=0.2559, zeta50=-9.4776e-04

      real,parameter  ::  Ahoi=0.1667, Bhoi=-0.01
      real  ::  depoltmp, perptmp1, perptmp2
      real,dimension(:,:),allocatable :: icecloud, watercloud,phasecloud
      real,dimension(:,:,:),allocatable :: uncloud
      real,dimension(:,:,:,:),allocatable :: icecloudfractday, watercloudfractday,indphaseday
      real,dimension(:,:,:,:),allocatable :: phasefractday,inddayphase
      real,dimension(:,:,:,:,:),allocatable :: uncloudfractday
      real,dimension(:,:,:,:),allocatable :: indphasefractday
      real,dimension(:,:,:,:,:),allocatable :: indphaseunday

      real,dimension(:,:,:),allocatable :: monthicecloud, monthwatercloud
      real,dimension(:,:,:),allocatable :: indphasemonth
      real,dimension(:,:,:,:),allocatable :: monthuncloud
      real,dimension(:,:,:),allocatable :: monthphasecloud
      real,dimension(:,:,:),allocatable  ::  monthcftemp,monthcftempice
      real,dimension(:,:,:),allocatable  ::  monthcftempliq
      real,dimension(:,:,:),allocatable  :: monthcftempphase,indmonthphasetemp
      real,dimension(:,:,:),allocatable  :: indcftemppermonth
      real,dimension(:,:,:,:),allocatable  :: indcftempphase,indcftemp
      real,dimension(:,:,:,:),allocatable  :: cftempday,cftempphaseday
      real,dimension(:,:,:,:),allocatable  :: cftempiceday,cftempliqday
      real,dimension(:,:),allocatable  ::  cftemp,cftempliq,cftempice
      real  ::  cfsumtemp

! variable pour fichier phase occurrences
real*4,dimension(:,:,:),allocatable  :: tot_ind,cloud_ind,ice_ind,water_ind
real*4,dimension(:,:,:,:),allocatable  :: un_ind



      metal='Lidar_Data_Altitudes'   ! name of meta var     
      metam='Met_Data_Altitudes' 



!****************************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************!
 

!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                  SELECTION OF INPUT/OUPUT PARAMETER                      *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

!------------------- Open output grid file -----------------------!
!                                                                 !
!   The format you have to enter is :                             !
!   "Map3D330m_",period,day/night,grid,version                    !
!   ex : "Map3D330m_200701_night_CFMIP_1.0"                       !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'output Prof file name = '
     read *, file5
     if (err==0) exit 
  enddo

!------------------- Open input data file ------------------------!
!                                                                 !
!   You have to enter the path and the name of the file  which    !
!   lists all the hdf CALIPSO                                     !
!   ex : "/homedata/mreverdy/GOCCP/Liste/200701_night"                       !    
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'input file name = '
     read *, file3
     open(unit=1,file=file3,iostat=err,status='OLD')
     if (err==0) exit
     print *,'--- input file not found' 
  enddo

!--------------------- Select the model --------------------------!
!                                                                 !
!   model = "lmdz" "chim" or "wrf" to select the LMDZ CHIMERE or  !
!           WRF model                                             !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'Enter the model : '
     read *, model
     if (err==0) exit 
  enddo
  
!----------------- Select the day or night -----------------------!
!                                                                 !
!   Select day or night Data version, because the process is not  !
!   the same during night and day                                 !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'Enter night or day : '
     read *, switch
     if (err==0) exit 
  enddo
 
!-------------------- Select the grid ----------------------------!
!                                                                 !
!    - "CFMIP" : 1° x 1° x 41levels each 480m from 0 to 19.2km    !
!             (361,181,41)                                        !
!    - "LMDZ"  : 3.75° x 2.53° x 19levels from 0 to 40.4823km     !
!             (96,72,19)                                          !
!    - "NASA"  : 5° x 5° x 41levels each 480m from 0 to 19.2km    !
!             (73,37,41)                                          !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'Enter the grid : '
     read *, gcm
     if (err==0) exit 
  enddo

!----------------- Select pressure or altitude -------------------!
!                                                                 !
!   Select "pressure" or "altitude" version of ouput data         !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'Enter vertical unit : '
     read *, alt_pres
     if (err==0) exit 
  enddo


!----------------- Select sat or cloudy mode ---------------------!
!                                                                 !
!   Select "sat" or "cloudy" mode in order to count the first     !
!   fully attenuated point as a cloudy point in the cloudfraction !
!   if "cloudy" mode is selected, and as a fully attenuated point !
!   if "sat" mode is selected. This change appears in the routine !
!   subgrid_fraction2.                                            !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'Enter sat or cloudy mode : '
     read *, switch2
     if (err==0) exit 
  enddo

print *, 'input parameters entered'

!****************************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************!

!_____________________________________________________________________________
!
! Temps CPU de calcul initial.
  call cpu_time(t_cpu_0)

! Temps elapsed de reference.
  call system_clock(count=t1, count_rate=ir)
!_____________________________________________________________________________
!


!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*             LOADING LON-LAT-ALT-SR-DEPOL-PR2 GRID VECTORS                *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

! Read parameter file
open(5,file='/homedata/mreverdy/GOCCP/'//trim(gcm)//'.p')
   read(5,*)grid, altmax2, altmax, latmax, lonmax, toplowl, topmidl, tophighl, altfile, latfile, lonfile
close(5)

print *, 'Grid parameter file read' 
allocate(heightmod(heightmax))
allocate(latmod(latmax),lonmod(lonmax),prestop(altmax),altmod(altmax2),srmod(diagmax),pr2mod(pr2max),atbrmod(permax), tempmod(tempmax),              &
         altmid(altmax),latmid(latmax-1),lonmid(lonmax-1), depolmod(depolmax), srdepmod(pr2max),altmod_bound(altmax,2))
allocate(tempmod_bound(tempmax-1,2),tempmid(tempmax-1))
  heightmod(:)=0;
  prestop(:)=0;latmod(:)=0;lonmod(:)=0;
  altmod(:)=0;srmod(:)=0; depolmod(:)=0; pr2mod(:)=0;atbrmod(:)=0;tempmod(:)=0;
  altmid(:)=0; latmid(:)=0; lonmid(:)=0;
  srdepmod(:)=0;

altmod_bound(:,:)=0; tempmod_bound(:,:)=0;tempmid(:)=0;
  
! loading the grid of the diagSR boxes value
open(17,file='/homedata/mreverdy/GOCCP/grilles_lmdz/srmod10')
!print *, 'open the file'

  do idiag=1,diagmax
     read(17,*)srmod(idiag)
  enddo
close(17)

! loading the grid of the DepolSR boxes value
open(6,file='/homedata/mreverdy/GOCCP/grilles_lmdz/depolmod')

  do idep=1,depolmax
     read(6,*)depolmod(idep)
  enddo
close(6)

! loading the grid of the pr2 boxes value
open(18,file='/homedata/mreverdy/GOCCP/grilles_lmdz/atbmod301')
  do ipr2=1,pr2max
     read(18,*)pr2mod(ipr2)
  enddo
close(18)

open(23,file='/homedata/mreverdy/GOCCP/grilles_lmdz/atbrmod241')
  do ipr2=1,permax
     read(23,*)atbrmod(ipr2)
  enddo
close(23)

! loading the grid of the DepolSR boxes value
open(20,file='/homedata/mreverdy/GOCCP/grilles_lmdz/tempmod39')

  do itemp=1,tempmax
     read(20,*)tempmod(itemp)
   enddo
close(20)

        do itemp=1,tempmax-1
            tempmid(itemp) = (tempmod(itemp)+tempmod(itemp+1))/2
            tempmod_bound(itemp,1)=tempmod(itemp);
            tempmod_bound(itemp,2)=tempmod(itemp+1);
         enddo

! Computing the Height grid
 do iheight=1,heightmax-1
    heightmod(iheight+1)=heightmod(iheight)+0.5
 enddo



! loading the level grid (altitude or pressure)
    if(alt_pres=='altitude')then

open(15,file='/homedata/mreverdy/GOCCP/grilles_lmdz/'//altfile)
         do iz=1,altmax2
            read(15,*)altmod(iz)
         enddo
         do iz=1,altmax
            altmid(iz) = (altmod(iz)+altmod(iz+1))/2
            altmod_bound(iz,1)=altmod(iz);
            altmod_bound(iz,2)=altmod(iz+1);
         enddo

         do iz=1,altmax     
            if(altmod(iz+1).GE.8.64)then
               seuilsnrhigh=iz
               exit
            endif
         enddo

         do iz=altmax,1,-1     
            if(altmod(iz+1).LE.2.4)then
                seuilsnrlow=iz
                exit
            endif
         enddo

    elseif(alt_pres=='pressure')then
         if(trim(gcm).eq.'LMDZ')then
            open(15,file='/homedata/mreverdy/GOCCP/grilles_lmdz/pression_lmdz2.txt')
         else if(trim(gcm).eq.'CFMIP')then
            open(15,file='/homedata/mreverdy/GOCCP/grilles_lmdz/pres_cfmip')
         elseif(trim(gcm).eq.'NASA')then
            open(15,file='/homedata/mreverdy/GOCCP/grilles_lmdz/pres_cfmip')   
 
         endif
            do iz=1,altmax
            read(15,*)prestop(iz)           ! lmdz milieu de la couche 
            enddo       
   endif
close(15)


open(10,file='/homedata/mreverdy/GOCCP/grilles_lmdz/'//lonfile)
      do ilon=1,lonmax
         read(10,*)lonmod(ilon)           !lmdz
      enddo
      do ilon=1,lonmax-1
         lonmid(ilon)=(lonmod(ilon)+lonmod(ilon+1))/2
      enddo
close(10)

open(21,file='/homedata/mreverdy/GOCCP/grilles_lmdz/'//latfile)
      do ilat=1,latmax
         read(21,*)latmod(ilat)           ! lmdz
      enddo
      do ilat=1,latmax-1
         latmid(ilat)=(latmod(ilat)+latmod(ilat+1))/2
      enddo
close(21)
  print *, 'latmod lonmod altmod srmod depolmod pr2mod ok'

!****************************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************!



!____________________________________________________________________________
!
! Initialization of box and numfich

      box=0
      numfich=0

! Reading the list of hdf files

 888  read(1,100,end=999)file2

      print *, 'Processing with ',trim(file2)
      numfich=numfich+1
      print *, 'lecture du fichier numero ',numfich

! initialization of var
   it=0; comptpf=0; date=0; year=0; month=0; day=0; 
filetmp2=trim(file2)


!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*             BEGINING OF THE READING OF THE SDS/META VAR                  *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

print *, 'Read the Calipso file data ...'


!****************************** READING SDS VAR ****************************!

call nprof(filetmp2,20,it)       ! find the number of profil it

! empty file checking
if(it.lt.500)then
goto 887
endif

! Allocation of interpolated variables
allocate(lat(it),lon(it),SE(it),temps(it),temps2(it), &
          mol2(altitude,it),   & 
         mol3(altitude, it), temp2(altitude, it), stat = OK_buffer)

if(alt_pres=='pressure')then
   allocate(pres2(altitude, it))
   pres2(:,:)=0
endif

! check the allocation 
 if (OK_buffer/=0) print *,'--- buffer allocation error '   

! Initialization of interpolated variables
   temps2(:)=0; temps(:)=0; mol2(:,:)=0; mol3(:,:)=0; temp2(:,:)=0;!mol4(:,:)=0;
   lat(:)=0;lon(:)=0;SE(:)=0;! indtot(:)=0; indret(:)=0;
 



!---------------- Select the type of hdf file ( prov/launch )----------------!
!                                                                            !
! Depending on the version of the hdf files (prov or launch), the time       !
! variable is different. That's why, the process isn't the same in the 2cases!
!                                                                            !
!----------------------------------------------------------------------------!

if(file2(67:70)/="Laun")then

   !     Retrieve data for time variable.
   sds_varname='Profile_UTC_Time'
   call sdsread8(vartmp8,filetmp2,sds_varname,ret1)
   
   if(ret1.eq.-1)then
      deallocate(vartmp8)
      deallocate(lat,lon,SE,temps,temps2,mol2, mol3, mol4,indtot, indret, stat = OK_buffer)
      deallocate(indretlow,indretmid,indrethigh)
      command3='echo '//trim(file2)//' >> /homedata/mreverdy/GOCCP/cal_corrompu'
      call system(trim(command3))
      goto 887
   endif
   	
   temps(:)=vartmp8(1,:);
   deallocate(vartmp8)


!------------------ Calculation of the date : type Prov ---------------------!
!                                                                            !
!   temps(nprof) is in International Atomic Time in UTC : yymmdd.ffffffff    !
!                                                                            !
!   yy  =  Last two digits of year where 07 represents 2007                  !
!   mm  =  Month in two-character subfield with values 01-12                 !
!   dd  =  Day of month in two-character subfield with values 01-28, -29,..  !
!   "."  =  Period as a separator                                            !
!   ffffffff  =  Fractional part of day                                      !
!                                                                            !
!----------------------------------------------------------------------------!

            date=int(temps(1))           ! date of the read file
            year=2000+int(date/10000)                      ! year
            month=int((date-int(date/10000)*10000)/100)    ! month
            day=date-int(date/10000)*10000-month*100       ! day
            jour=day
            print *, 'Processing for the ',day,'/',month,'/',year
            print *, 'File type = Prov'

else 

!------------ Calculation of the date : type Launch (old data) --------------!
!                                                                            !
! If the file is a launch one then, the UTC time is calculated from the TAI  !
! time in seconds.                                                           !
! This calculation is necessary for the June July and August months period.  !
!                                                                            !
!----------------------------------------------------------------------------!

   sds_varname='Profile_Time'
   call sdsread8(vartmp8,filetmp2,sds_varname,ret1)
   temps(:)=vartmp8(1,:);
   deallocate(vartmp8)


   date=int(temps(1)) 
   year=2006

   if(date.lt.425865606)then   ! first value of the time in June
      month=6
      day=int((date-423273606)/86400)+1  
   elseif(date.lt.428544006)then   ! first value of time in July
         month=7
         day=int((date-425865606)/86400)+1
   else
         month=8
         day=int((date-428544006)/86400)+1   ! first value of time in August
   endif

   print *, 'Processing for the ',day,'/',month,'/',year
   print *, 'Fichier Launch'

endif ! end selection of file type


!---------------- Time in day for the output netcdf files -------------------!
!                                                                            !
! Select the number of days since 2000-01-01 00:00:00 to the mid of the run  !
! period.                                                                    !
!                                                                            !
! ex : for monthly period, date is set to the 15th of the month              !
!      for a trimonthly period date is set to the 15th of the second month   !
!                                                                            !
!----------------------------------------------------------------------------!

if(numfich.eq.1)then 
     print *, file3(25:30)
     date1=((year*100+month)*100+15)*100
     date2=2000010100
     print *, date1
     print *, date2

     resh=interdat(date2,date1) ! calculation of time between the 2 period 
     print *, resh              ! in hour
     resd=resh/24               ! converted in days
     print *, resd
endif


!     Retrieving data for atb variable.
sds_varname='Total_Attenuated_Backscatter_532'
call sdsread(atb,filetmp2,sds_varname)

!     Retrieving data for atb variable.
sds_varname='Attenuated_Backscatter_1064'
call sdsread(atb2,filetmp2,sds_varname)

!     Retrieving data for atb variable.
sds_varname='Perpendicular_Attenuated_Backscatter_532'
call sdsread(perp,filetmp2,sds_varname)

!     Retrieving data for lat variable.
sds_varname='Latitude'
call sdsread(vartmp,filetmp2,sds_varname)
lat(:)=vartmp(1,:);
deallocate(vartmp)

!     Retrieving data for lon variable.
sds_varname='Longitude'
call sdsread(vartmp,filetmp2,sds_varname)
lon(:)=vartmp(1,:);
deallocate(vartmp)

!     Retrieving data for mol variable.
sds_varname='Molecular_Number_Density'
call sdsread(mol,filetmp2,sds_varname) 

!     Retrieve data for temp variable.
sds_varname='Temperature'
call sdsread(temp,filetmp2,sds_varname)


!     Retrieving data for surf_elevation !!
sds_varname='Surface_Elevation'
call sdsread(vartmp,filetmp2,sds_varname)
SE(:)=vartmp(1,:);
deallocate(vartmp)


!***************************** READING META VAR *****************************!

!     Retrieving data for altitude variable.
call metaread(altl,metal,filetmp2)
!     Retrieving data for altitude variable.
call metaread(altm,metam,filetmp2)


print *, 'HDF Calipso File read'
print *, 'Input variables read'
print *, ''
!_____________________________________________________________________________
!



print *, 'Interpolation of data & molecular calculation start ...'

!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*           BEGINNING OF CALCULATION ON THE INPUT VARIABLES                *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!



!************** CONVERTION OF THE UTC FRACTION IN UT HOUR *******************!


!******************************* WRF  MODE *********************************!
if(model.eq.'wrf')then

do i=1,it          !!!! loop on each profil   
   if((file2(67:70)=="Prov").or.(file2(67:70)=="AL_L"))then
   ! Calculation of hour UT with UTC time fraction
      temps2(i)=(temps(i)-date)*24  ! temps = fraction of hour in prov file type
   else 
   ! Calculation of hour UT with TAI time
      if(file2(67:70)=="Laun")then
         if(date.lt.425865606)then
            temps2(i)=(((temps(i)-423273606)/86400+1)-day)*24
         elseif(date.lt.428544006)then
            temps2(i)=(((temps(i)-425865606)/86400+1)-day)*24    
         else
            temps2(i)=(((temps(i)-428544006)/86400+1)-day)*24  
         endif
      endif
   endif
enddo

!else
!   deallocate(temps)
endif
!_____________________________________________________________________________
!


!****************************** LAUNCH FILE *********************************!   
if(file2(67:70)=="Laun")then
   do i=1,it          !!!! loop on each profil
      jour=day
      date=60000+day+100*month
   enddo
endif
!_____________________________________________________________________________
!


!  INTERPOLATION OF PRESSURE TEMPERATURE AND MOLECULAR FROM 33LVL TO 583LVL  !
if(alt_pres=='pressure')then
 do i=1,it          !!!! loop on each profil
   call interp(pres,pres2,altm,altl,i,it)
   call interp(mol,mol2,altm,altl,i,it)     ! atb molecular without dimension

 enddo        !!!! end of profil loop
else
 do i=1,it          !!!! loop on each profil
   temps2(i)=(temps(i)-date)*24  ! temps = fraction of hour in prov file type
   call interp(mol,mol2,altm,altl,i,it)     ! atb molecular without dimension
   call interp(temp,temp2,altm,altl,i,it)
   call atb_temp_interp(temp2,altl,i,it,SeuilTemp1km,SE)         ! molecular extrapolation
   call SE_alt_mol(SE,altl,temp2,altitude,it,i)! add the SE to mol3

enddo        !!!! end of profil loop
endif


deallocate(temps)

!************************ CALCULATION OF ATB-MOLECULAR **********************!
 

!---------------------------- day or night mode -----------------------------!
!                                                                            !
!  The calculation of normalized ratio is calculated at different altitude   !
!  levels whether it is day or night time.                                   !
!  This difference appears in the "atb_mol" routine:                         !
!     - between 20km and 25km during day time                                !
!     - between 22km and 25km during night time                              !
!                                                                            !
!----------------------------------------------------------------------------!

!-------------------------------- PSC filter --------------------------------!
!                                                                            !
!  During the polar winter, the cloud could appear up to 20km altitude at    !
!  latitude higher than 60° south. That's why from June to October, the      !
!  altitude calculation of the normalized ratio must be higher at these      !
!  latitudes.                                                                !
!  The calculation altitude is ranged :                                      !
!     - between 28.5km and 35km during day time                              !
!     - between 28.5km and 33km during night time                            !
!                                                                            !
!----------------------------------------------------------------------------!
if(switch.eq.'day')then   !!! DAY mode  !!!

 do i=1,it          !!!! loop on each profil   
    ! PSC filter during the polar winter
   if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then 
      ! normalized ratio betwen 28.5 & 35km & molecr calculation in km-1 sr-1  
      call atb_mol(atb,mol2,mol3,i,it,17,42) ! Antarctic PSCs season 
   elseif((lat(i).ge.60).and.((month.eq.12).or.(month.eq.1).or.(month.eq.2)))then
      call atb_mol(atb,mol2,mol3,i,it,17,42) ! Arctic PSCs season   Pitts et al 2011   
   else
      ! normalized ratio betwen 20 & 25km & molecular calculation in km-1 sr-1
      call atb_mol(atb,mol2,mol3,i,it,62,92) 
   endif



      call atb_mol_interp(mol3,altl,i,it,SeuilMol1km,SE)         ! molecular extrapolation
      call SE_alt_mol(SE,altl,mol3,altitude,it,i)! add the SE to mol3
      call SE_alt_atb(SE,altl,atb,altitude,it,i) ! add the SE to atb

 enddo         !!!! end of profil loop

elseif(switch.eq.'night')then  !!! NIGHT mode !!!

 do i=1,it          !!!! loop on each profil  
  ! PSC filter during the polar winter
  if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then
     ! normalized ratio betwen 28.5 & 33km & molecr calculation in km-1 sr-1  
     call atb_mol(atb,mol2,mol3,i,it,24,42)  
  elseif((lat(i).ge.60).and.((month.eq.12).or.(month.eq.1).or.(month.eq.2)))then
      call atb_mol(atb,mol2,mol3,i,it,24,42) ! Arctic PSCs season   Pitts et al 2011
  else
     ! normalized ratio betwen 22& 25km & molecular calculation in km-1 sr-1
     call atb_mol(atb,mol2,mol3,i,it,62,78)  
  endif

!mol4(:,i)=mol3(:,i); ! keep the mol3 before extrapolation in order to count the retrieval
      call atb_mol_interp(mol3,altl,i,it,SeuilMol1km,SE)        ! molecular extrapolation
      call SE_alt_mol(SE,altl,mol3,altitude,it,i) ! add the SE to mol3
      call SE_alt_atb(SE,altl,atb,altitude,it,i)  ! add the SE to atb
 enddo         !!!! end of profil loop


endif  ! end of day/night loop

print *, 'Interpolation of data & molecular calculation done'
!_____________________________________________________________________________
!




!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*            SELECTION OF THE MODEL GRID CFMIP / LMDZ / NASA               *!
!*               AND ALLOCATION OF THE OUTPUTS VARIABLES                    *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

    
 file6='MapLowMidHigh'//trim(file5(17:))//trim(version)    ! Output Map file name
 file66='MapHigh'//trim(file5(17:))//trim(version)    ! Output Map file name
file7='SR_histo'//trim(file5(17:)) //trim(version)          ! Output diagSR file name
 file12='SR_histo_Phase'//trim(file5(17:)) //trim(version)          ! Output diagSR file name
 file10='3D_CloudFraction_Phase'//trim(file5(17:))//trim(version)    ! Output depolSR file name
file13='3D_CloudFraction_Temp'//trim(file5(17:))//trim(version)  
 file11='MapLowMidHigh_Phase'//trim(file5(17:))//trim(version)    ! Output Map file name
! file12='Phase_histo'//trim(file5(17:))//trim(version)    ! Output Map file name


print *, 'altmax=',altmax


! Allocation / initialization of variable 
print *, 'Allocation / Initialization of variables...'

  allocate(indicep(altmax,it),indicep2(altmax,it),pr2moy2(altmax,it),indice2(altmax,it),crmoy(altmax,it))
  allocate(pr2moy(altmax,it),molmoy(altmax,it),srmoy(altmax,it),depolmoy(altmax,it),perpmoy(altmax,it),parmoy(altmax,it),tempmoy(altmax,it))
  allocate(indice(altmax,it),indicem(altmax,it),indicetemp(altmax,it)) 


if(numfich.eq.1)then ! Allocation of monthly variables


   allocate(indday(latmax-1,lonmax-1,altmax,daymax))
   allocate(cloudfractday(latmax-1,lonmax-1,altmax,daymax),                      &
            clearfractday(latmax-1,lonmax-1,altmax,daymax))
   allocate(uncertfractday(latmax-1,lonmax-1,altmax,daymax))
   allocate(icecloudfractday(latmax-1,lonmax-1,altmax,daymax),               &
            phasefractday(latmax-1,lonmax-1,altmax,daymax),                      &
            watercloudfractday(latmax-1,lonmax-1,altmax,daymax))
   allocate(uncloudfractday(latmax-1,lonmax-1,altmax,daymax,catmax))

   allocate(indphaseday(latmax-1,lonmax-1,altmax,daymax))
   allocate(indphaseunday(latmax-1,lonmax-1,altmax,daymax,catmax))

   allocate(isccplowday(latmax-1,lonmax-1,daymax),                               &
            isccpmidday(latmax-1,lonmax-1,daymax),                               &
            isccphighday(latmax-1,lonmax-1,daymax))

   allocate(heightday(latmax-1,lonmax-1,daymax),                                &
            indheight(latmax-1,lonmax-1,daymax))
   allocate(heightday2(latmax-1,lonmax-1,daymax))

   allocate(colcloudday(latmax-1,lonmax-1,daymax),                               &
   
            isccpinddaylow(latmax-1,lonmax-1,daymax),                            &
            isccpinddaymid(latmax-1,lonmax-1,daymax),                            &
            isccpindday(latmax-1,lonmax-1,daymax))


allocate(cftempday(latmax-1,lonmax-1,tempmax-1,daymax), &
         cftempiceday(latmax-1,lonmax-1,tempmax-1,daymax), &
         cftempliqday(latmax-1,lonmax-1,tempmax-1,daymax), &
         indcftemp(latmax-1,lonmax-1,tempmax-1,daymax), &
         indcftempphase(latmax-1,lonmax-1,tempmax-1,daymax))

   allocate(lowtemp(latmax-1,lonmax-1,daymax),                         &
            midtemp(latmax-1,lonmax-1,daymax),                          &
            hightemp(latmax-1,lonmax-1,daymax),                          &
            coltemp(latmax-1,lonmax-1,daymax))
   allocate(indlowtemp(latmax-1,lonmax-1,daymax),                         &
            indmidtemp(latmax-1,lonmax-1,daymax),                          &
            indhightemp(latmax-1,lonmax-1,daymax),                          &
            indcoltemp(latmax-1,lonmax-1,daymax))



   allocate(diagSR(lonmax-1,latmax-1,altmax,diagmax-1))
   allocate(diagSRpha(lonmax-1,latmax-1,altmax,diagmax-8,3))

   allocate(diagPHA(lonmax-1,latmax-1,altmax,tempmax-1,3))

   allocate(indnan(latmax-1,lonmax-1,altmax))


allocate(indtotmean(lonmax-1,latmax-1),indtot(lonmax-1,latmax-1))
   allocate(isccpliqday(latmax-1,lonmax-1,daymax,4),isccpiceday(latmax-1,lonmax-1,daymax,4))
   allocate(isccpunday(latmax-1,lonmax-1,daymax,4,catmax),                               &
            isccpphaseday(latmax-1,lonmax-1,daymax,4))

endif

if(alt_pres=='pressure')then
   allocate(SEp(it))
   SEp(:)=0
endif

  pr2moy(:,:)=0; molmoy(:,:)=0;tempmoy(:,:)=0;  
  srmoy(:,:)=0; depolmoy(:,:)=0;
  indice(:,:)=0; indicem(:,:)=0;
  pr2moy2(:,:)=0;indice2(:,:)=0;crmoy(:,:)=0;
  indicep2(:,:)=0; indicep(:,:)=0; indicetemp(:,:)=0;

if(numfich.eq.1)then 
   indday(:,:,:,:)=0;
  cloudfractday(:,:,:,:)=0;clearfractday(:,:,:,:)=0;
  uncertfractday(:,:,:,:)=0;
  icecloudfractday(:,:,:,:)=0; watercloudfractday(:,:,:,:)=0
  uncloudfractday(:,:,:,:,:)=0;
  phasefractday(:,:,:,:)=0; 

indcftempphase(:,:,:,:)=0
indcftemp(:,:,:,:)=0
cftempday(:,:,:,:)=0
cftempiceday(:,:,:,:)=0
cftempliqday(:,:,:,:)=0

  indphaseday(:,:,:,:)=0; indphaseunday(:,:,:,:,:)=0; 
  heightday(:,:,:)=0; indheight(:,:,:)=0;
  heightday2(:,:,:)=0;
  isccplowday(:,:,:)=0;isccpmidday(:,:,:)=0;isccphighday(:,:,:)=0;
  colcloudday(:,:,:)=0;isccpindday(:,:,:)=0;
  isccpinddaylow(:,:,:)=0;isccpinddaymid(:,:,:)=0;
  diagSR(:,:,:,:)=0;
  diagSRpha(:,:,:,:,:)=0;

  diagPHA(:,:,:,:,:)=0;

  isccpliqday(:,:,:,:)=0;isccpiceday(:,:,:,:)=0;
  isccpunday(:,:,:,:,:)=0; isccpphaseday(:,:,:,:)=0;

lowtemp(:,:,:)=0;midtemp(:,:,:)=0;hightemp(:,:,:)=0;coltemp(:,:,:)=0;
indlowtemp(:,:,:)=0;indmidtemp(:,:,:)=0;indhightemp(:,:,:)=0;indcoltemp(:,:,:)=0;


indtotmean(:,:)=0; indtot(:,:)=0;

endif
print *, ''
!_____________________________________________________________________________
!


!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*       AVERAGING OF OBSERVATIONS DATA OVER THE VERTICAL MODEL BOXES       *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


print *, 'Begin of vertical average'

!**************** VERTICAL AVERAGE OF OBSERVATIONS DATA *********************!



if(alt_pres=='altitude')then


   do i=1, it         !!! BEGIN OF IT LOOP
      do iz=altmax,1,-1 
         do ilid=1,altitude
            if ( (altl(ilid).ge.altmod(iz)).and.(altl(ilid).lt.altmod(iz+1)))&
            then   

             call vertical_mean(temp2,tempmoy,atb,indicetemp,i,iz,ilid,it,altmax,&
                                 altitude,4)    
              call vertical_mean(atb,pr2moy,atb,indice,i,iz,ilid,it,altmax,  &
                                  altitude,1)      
              call vertical_mean(mol3,molmoy,atb,indicem,i,iz,ilid,it,altmax,&
                                  altitude,3)
                call vertical_mean(atb2,pr2moy2,atb,indice2,i,iz,ilid,it,altmax,&
                                  altitude,3)    
               call vertical_mean(perp,perpmoy,atb,indicep,i,iz,ilid,it,altmax,&
                                  altitude,3)    
               call vertical_mean(perp,parmoy,atb,indicep2,i,iz,ilid,it,altmax,&
                                  altitude,2)                   

            endif
         enddo
      enddo

      do iz=altmax,1,-1 
         call zero_detect(pr2moy,i,iz,it,altmax)
         call zero_detect(molmoy,i,iz,it,altmax)
         call zero_detect(tempmoy,i,iz,it,altmax)  
         call zero_detect(pr2moy2,i,iz,it,altmax)
         call zero_detect(parmoy,i,iz,it,altmax)
         call zero_detect(perpmoy,i,iz,it,altmax)
      enddo

          call Surf_detect2(SE,altmod,pr2moy,altmax,it,i,alt_pres)
          call Surf_detect2(SE,altmod,molmoy,altmax,it,i,alt_pres)  
          call Surf_detect2(SE,altmod,tempmoy,altmax,it,i,alt_pres)
          call Surf_detect2(SE,altmod,pr2moy2,altmax,it,i,alt_pres)
          call Surf_detect2(SE,altmod,parmoy,altmax,it,i,alt_pres)
          call Surf_detect2(SE,altmod,perpmoy,altmax,it,i,alt_pres)

     do iz=altmax,1,-1 
       if((trim(gcm).eq.'LMDZ').or.(trim(gcm).eq.'WRF'))then
         if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then
            call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it, &
                                 altmax)
            call SR_CR_DEPOL_mean(pr2moy2,pr2moy,crmoy,indice2,indice,i,iz,it,   &
                                 altmax)
            call SR_CR_DEPOL_mean(parmoy,perpmoy,depolmoy,indicep2,indicep,i,iz,it,   &
                                 altmax)
         else
            call filtre_2lvl(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,      &
                            altmax,gcm)
            call SR_CR_DEPOL_mean(pr2moy2,pr2moy,crmoy,indice2,indice,i,iz,it,   &
                                 altmax)
            call SR_CR_DEPOL_mean(perpmoy,parmoy,depolmoy,indicep,indicep2,i,iz,it,   &
                                 altmax)
         endif  
       
       elseif(trim(gcm).eq.'LMDZ40')then
         if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then
           call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it, &
                                 altmax)
         else
           call filtre_2lvl(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,      &
                            altmax,gcm)
         endif
       
       elseif((trim(gcm).eq.'CFMIP').or.(trim(gcm).eq.'CFMIP1').or.(trim(gcm)&
       .eq.'CFMIP2.5').or.(trim(gcm).eq.'CFMIP80').or.(trim(gcm).eq.'CFMIP2').or.(trim(gcm).eq.'CFMIP160').or.(trim(gcm).eq.'CFMIP320'))then


         call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   &
                               altmax)
         call SR_CR_DEPOL_mean(pr2moy2,pr2moy,crmoy,indice2,indice,i,iz,it,   &
                               altmax)
         call SR_CR_DEPOL_mean(perpmoy,parmoy,depolmoy,indicep2,indicep,i,iz,it,   &
                               altmax)


       elseif(trim(gcm).eq.'NASA')then 
         call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   &
                               altmax)
         !call filtre_2lvl(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,altmax) 
       endif
     enddo

   !    Detection of the Surface Elevation & set -888 for value under this 
   !    threshold

      call Surf_detect(SE,altmod,srmoy,altmax,it,i)
      call Surf_detect(SE,altmod,crmoy,altmax,it,i)
      call Surf_detect(SE,altmod,depolmoy,altmax,it,i)
      

   enddo           !!! END OF IT LOOP


elseif(alt_pres=='pressure')then
   
do i=1, it         !!! BEGIN OF IT LOOP
      do iz=altmax,1,-1
         do ilid=1,altitude
         if ( (pres2(ilid,i).gt.prestop(iz)).and.(pres2(ilid,i).lt.          &
            prestop(iz-1)) )then
              call vertical_mean(atb,pr2moy,atb,indice,i,iz,ilid,it,altmax,  &
                                  altitude,1)     
              call vertical_mean(mol3,molmoy,atb,indicem,i,iz,ilid,it,altmax,&
                                  altitude,3)    
            endif
         enddo
      enddo
      do iz=altmax,1,-1
         call zero_detect(pr2moy,i,iz,it,altmax)
         call zero_detect(molmoy,i,iz,it,altmax)
      enddo


         call SE_km_2_Pres2(SE,SEp,altl,pres2,prestop,pr2moy,altmax,it,i) 
         call SE_km_2_Pres2(SE,SEp,altl,pres2,prestop,molmoy,altmax,it,i) 

     do iz=altmax,1,-1 
       if((trim(gcm).eq.'LMDZ').or.(trim(gcm).eq.'WRF'))then
         if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then
            call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it, &
                                 altmax)
         else
            call filtre_2lvl(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,      &
                            altmax,gcm)
         endif  
       
       elseif(trim(gcm).eq.'LMDZ40')then
         if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then
           call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it, &
                                 altmax)
         else
           call filtre_2lvl(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,      &
                            altmax,gcm)
         endif
       
       elseif((trim(gcm).eq.'CFMIP').or.(trim(gcm).eq.'CFMIP1').or.(trim(gcm)&
       .eq.'CFMIP2.5').or.(trim(gcm).eq.'CFMIP80').or.(trim(gcm).eq.'CFMIP2').or.(trim(gcm).eq.'CFMIP160').or.(trim(gcm).eq.'CFMIP320'))then


    if ((pr2moy(iz,i).eq.(-9999)).or.(molmoy(iz,i).eq.(-9999))) &
       then
       continue
     else


         do ilid=1,altitude
             if ( (altl(ilid).gt.altmod(iz)).and.(altl(ilid).lt.altmod(iz+1)))&
                  then  
                     if(mol4(ilid,i).ne.-9999)then
                        indret(i)=indret(i)+1
                     endif
             endif
          enddo




        if (iz.lt.toplowl) then
          do ilid=1,altitude
             if ( (altl(ilid).gt.altmod(iz)).and.(altl(ilid).lt.altmod(iz+1)))&
                  then  
                     if(mol4(ilid,i).ne.-9999)then
                        indretlow(i)=indretlow(i)+1
                     endif
             endif
          enddo
       endif

       if ((iz.ge.toplowl).and.(iz.le.topmidl))then
          do ilid=1,altitude
             if ( (altl(ilid).gt.altmod(iz)).and.(altl(ilid).lt.altmod(iz+1)))&
                  then  
                     if(mol4(ilid,i).ne.-9999)then
                        indretmid(i)=indretmid(i)+1
                     endif
             endif
          enddo
       endif

       if(iz.gt.topmidl) then
          do ilid=1,altitude
             if ( (altl(ilid).gt.altmod(iz)).and.(altl(ilid).lt.altmod(iz+1)))&
                  then  
                     if(mol4(ilid,i).ne.-9999)then
                        indrethigh(i)=indrethigh(i)+1
                     endif
             endif
          enddo
       endif

       test=indretlow(i)+indretmid(i)+indrethigh(i)
       if(test.ne.indret(i))print *, 'error indice retrieval', indret

    endif

         call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   &
                               altmax)
         call SR_CR_DEPOL_mean(pr2moy2,pr2moy,crmoy,indice2,indice,i,iz,it,   &
                               altmax)
         call SR_CR_DEPOL_mean(parmoy,perpmoy,depolmoy,indicep2,indicep,i,iz,it,   &
                               altmax)

       elseif(trim(gcm).eq.'NASA')then 
         call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   &
                               altmax)
   
       endif
     enddo


      call SE_km_2_Pres(SE,SEp,altl,pres2,prestop,srmoy,altmax,it,i) 
 enddo
endif            !!! END OF IT LOOP


! Mode wrf model
if(model.eq.'wrf')goto 622

print *, 'deallocate input var'


! Deallocate SDS/META variables
deallocate(atb,stat = OK_buffer)
deallocate(atb2,stat = OK_buffer)
deallocate(perp,stat = OK_buffer)
deallocate(mol,stat = OK_buffer)

deallocate(mol2,stat = OK_buffer)
deallocate(mol3,stat = OK_buffer)
deallocate(temp,stat = OK_buffer)
deallocate(temp2,stat = OK_buffer)



if(alt_pres=='pressure')deallocate(pres2,stat = OK_buffer)
deallocate(altl,stat = OK_buffer)
deallocate(altm,stat = OK_buffer)



!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                 CLOUD DIAGNOSTICS FOR ONE PROFIL                         *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

print *, 'allocate flag var'
! Allocation / initialization of instantaneous fraction variables
   allocate(cloudfraction(altmax,it),clearfraction(altmax,it),satfraction    &
            (altmax,it), uncertfraction(altmax,it),rejfraction(altmax,it))
   allocate(nanfraction(altmax,it),sefraction(altmax,it))!,fractot(altmax,it))

cloudfraction(:,:)=0; clearfraction(:,:)=0; satfraction(:,:)=0; 
uncertfraction(:,:)=0; nanfraction(:,:)=0;sefraction(:,:)=0;! fractot(:,:)=0 
rejfraction(:,:)=0;

print *, 'plop'



!****************************************************************************!
!******************** INSTANTANEOUS FRACTION DIAGNOSTIC *********************! 
!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                   RECORD INSTANTANEOUS SR FILES                          *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

seuilatb=2.5e-03
deltatb=0

! select the instant file type: 
!  - on = instant classic file activated
!  - off = desactivated
!  - fraction = instant fraction file with GOCCP cloud mask
!instant_switch='on'

SELECTCASE('instant_switch')


!****************************************************************************!
!****************************************************************************!
!****************************************************************************!
!*                               RECORD ON                                  *!
!****************************************************************************!
!****************************************************************************!
!****************************************************************************!
! instant classic file with SR value
CASE('on')


allocate(cloudfraction2(altmax,it))
cloudfraction2(:,:)=0;

  write(numfichc,'(i4)')numfich
  write(datec,'(I6)')date
  write(yearc,'(I4)')year

command4='echo '//trim(file2)//'| cut -d/ -f8 | cut -d. -f2 > /homedata/mreverdy/GOCCP/src/instantname'//yearc//datec(3:6)//'_'//trim(switch)//'_'//trim(gcm)


CALL SYSTEM(trim(command4))
open(10,file='/homedata/mreverdy/GOCCP/src/instantname'//yearc//datec(3:6)//'_'//trim(switch)//'_'//trim(gcm))
read(10,*)instantname
close(10)
!print *, 'avant routine',srmoy(7,49261)
!print *, molmoy(7,49261),pr2moy(7,49261)

do i=1,it       !!!!! BEGIN OF IT LOOP 
! flag 0/1        

  call fraction_subgrid2_8km(seuilsnrlow,seuilsnrhigh,srmoy,pr2moy,indice,   &
                             molmoy,indicem,satfraction, &
                             cloudfraction,clearfraction,uncertfraction,     &
                             nanfraction,sefraction,rejfraction,i,altmax,it, &
                             toplowl,topmidl,switch,switch2) 
  call fraction_subgrid3_8km(seuilsnrlow,seuilsnrhigh,altmod,srmoy,pr2moy,     &
                             indice,molmoy,indicem,cloudfraction2,i,altmax,it,  &
                             switch,switch2)  



enddo      !!! END OF IT LOOP   



where(rejfraction.eq.1)
srmoy=-777.
crmoy=-777.
depolmoy=-777.
pr2moy=-777.
molmoy=-777.
parmoy=-777.
perpmoy=-777.
!tempmoy=-777.
endwhere

where(nanfraction.eq.1)
srmoy=-9999.
crmoy=-9999.
depolmoy=-9999.
pr2moy=-9999.
molmoy=-9999.
parmoy=-9999.
perpmoy=-9999.
!tempmoy=-9999.
endwhere

where(sefraction.eq.1)
srmoy=-888.
crmoy=-888.
depolmoy=-888.
pr2moy=-888.
molmoy=-888.
parmoy=-888.
perpmoy=-888.
tempmoy=-888.
endwhere

where((srmoy.ne.-9999.).and.(srmoy.ne.-777.).and.(srmoy.ne.-888.))
pr2moy=pr2moy/indice
molmoy= molmoy/indicem
parmoy= parmoy/indicep2
perpmoy= perpmoy/indicep
endwhere

where(tempmoy.ne.-888.)
tempmoy= tempmoy/indicetemp
endwhere

print *, 'recording instant SR CR DR file'


!!!!! RECORD INSTANT SR FILES WITH ATB ATBper ATBpar ATBmol
  file4='instant_SR_CR_DR_'//trim(instantname)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(version)//'.nc'

print *, file4
  call SR_CR_DR_ATB_nc(file4,altmid,altmod_bound,resd,altmax,switch,gcm,it,lat,lon,SE,temps2,&
                  srmoy,crmoy,depolmoy,pr2moy,molmoy,perpmoy,parmoy,tempmoy,cloudfraction2)


deallocate(cloudfraction2)


!****************************************************************************!
!****************************************************************************!
!****************************************************************************!
!*                               RECORD OFF                                 *!
!****************************************************************************!
!****************************************************************************!
!****************************************************************************!
! Do not record any instant file
CASE('off')
do i=1,it       !!!!! BEGIN OF IT LOOP 
! flag 0/1        

  call fraction_subgrid2_8km(seuilsnrlow,seuilsnrhigh,srmoy,pr2moy,indice,   &
                             molmoy,indicem,satfraction, &
                             cloudfraction,clearfraction,uncertfraction,     &
                             nanfraction,sefraction,rejfraction,i,altmax,it, &
                             toplowl,topmidl,switch,switch2) 


!*************** instant SR corrected by delta atb ****************!

enddo      !!! END OF IT LOOP   


where(rejfraction.eq.1)
srmoy=-777.
crmoy=-777.
depolmoy=-777.
pr2moy=-777.
molmoy=-777.
parmoy=-777.
perpmoy=-777.
!tempmoy=-777.
endwhere

where(nanfraction.eq.1)
srmoy=-9999.
crmoy=-9999.
depolmoy=-9999.
pr2moy=-9999.
molmoy=-9999.
parmoy=-9999.
perpmoy=-9999.
!tempmoy=-9999.
endwhere

where(sefraction.eq.1)
srmoy=-888.
crmoy=-888.
depolmoy=-888.
pr2moy=-888.
molmoy=-888.
parmoy=-888.
perpmoy=-888.
tempmoy=-888.
endwhere

where((srmoy.ne.-9999.).and.(srmoy.ne.-777.).and.(srmoy.ne.-888.))
pr2moy=pr2moy/indice
molmoy= molmoy/indicem
parmoy= parmoy/indicep2
perpmoy= perpmoy/indicep
endwhere

where(tempmoy.ne.-888.)
tempmoy= tempmoy/indicetemp
endwhere






ENDSELECT

if(allocated(temps2)) deallocate(temps2,stat = OK_buffer)
deallocate(SE,stat = OK_buffer)



!************** INSTANTANEOUS ISCCP LOW MID HIGH FRACTION ******************!
!!!!!!!!!!!!!!!!!! DEFINITION OF TOP AND BASE LAYER !!!!!!!!!!!!!!!!!!!!!!!!!
! low level < 3.2km
! 3.2 <= mid level < 6.5km
! high level >= 6.5km

!  680 hPa ===> 3.5km avec équilibre hydrostatique 
                        !  P=P0.exp(-z/H) , H=8.5
           !  440 hPa ===> 7.2km
          ! top lvl of isccp 

print *, 'diagnostic fraction nuage subgrid'


! Allocate / initialization of instantaneous isccp variables
 allocate(isccplow(it),isccpmid(it),isccphigh(it),colcloud(it))!, colclear(it))
 allocate(watercloud(altmax,it),icecloud(altmax,it),uncloud(altmax,it,catmax),phasecloud(altmax,it))
allocate(cftemp(tempmax-1,it),cftempliq(tempmax-1,it),cftempice(tempmax-1,it))
!cftempice

 allocate(height(it),height2(it),isccpliq(4,it),isccpice(4,it),isccpun(4,it,catmax))


   height(:)=0; height2(:)=0;
   icecloud(:,:)=0; watercloud(:,:)=0; uncloud(:,:,:)=0;phasecloud(:,:)=0;
   isccplow(:)=0; isccpmid(:)=0; isccphigh(:)=0;
   colcloud(:)=0; !colclear(:)=0;
   isccpice(:,:)=0; isccpliq(:,:)=0; isccpun(:,:,:)=0;

cftemp(:,:)=0; 
cftempliq(:,:)=0; 
cftempice(:,:)=0;


!! looking for limit before the 2 different SNR
do  iz=altmax,1,-1 
   if(altmod(iz).lt.8.16)then
      toplvlsat1=iz+1
      exit
   endif      
enddo

!print *, "toto1"

!nol=nol_switch

do i=1,it       !!!!! BEGIN OF IT LOOP 

   nanprof=0
   perptmp1 =0


! CLoud Phase diagnostic with the equation of the liquid relation between ATB
! and ATBr: ATBr = 1.3919 * ATB² + 0.0176 * ATB

 do  iz=altmax,toplvlsat1,-1  
   if(cloudfraction(iz,i).gt.0.)then
       if(perpmoy(iz,i).gt.parmoy(iz,i))then   ! noisy point = unphysical value

          uncloud(iz,i,5)=uncloud(iz,i,5)+1.   
          phasecloud(iz,i)=7.
          
       else

!! Curve for discrimination between ice & liquid particles
           perptmp1 = (pr2moy(iz,i)**5)*alpha50 + (pr2moy(iz,i)**4)*beta50 + (pr2moy(iz,i)**3)*gamma50 + (pr2moy(iz,i)**2)*delta50 + pr2moy(iz,i)*epsilon50 + zeta50

!! Curve for discrimination between ho particles & other
           perptmp2 = pr2moy(iz,i)*Ahoi + Bhoi

            if( (perpmoy(iz,i)-perptmp1).ge.0. )then
                 if (tempmoy(iz,i).gt.0 ) then
                    uncloud(iz,i,3)=uncloud(iz,i,3)+1.    !fake ice particles  
                    watercloud(iz,i)=watercloud(iz,i)+1.  !fake ice particles 
                    phasecloud(iz,i)=5.
                 else
                    icecloud(iz,i)=icecloud(iz,i)+1.      ! ice particles
                    phasecloud(iz,i)=2.
                 endif

            elseif( (perpmoy(iz,i)-perptmp2).ge.0. )then
                    if (tempmoy(iz,i).gt.-42 ) then
                       watercloud(iz,i)=watercloud(iz,i)+1.
                       phasecloud(iz,i)=1.
                    else
                       icecloud(iz,i)=icecloud(iz,i)+1.      ! fakeliq particles
                       uncloud(iz,i,2)=uncloud(iz,i,2)+1.    !fakeliq particles  
                       phasecloud(iz,i)=4.                   !fakeliq particles
                    endif
                 else
                   uncloud(iz,i,4)=uncloud(iz,i,4)+1. !! horizontally oriented particles
                   phasecloud(iz,i)=6.              !! horizontally oriented particles
            endif
        
         endif
    endif
  enddo 


 do  iz=toplvlsat1-1,1,-1  
   if(cloudfraction(iz,i).gt.0.)then
       if(perpmoy(iz,i).gt.parmoy(iz,i))then   ! noisy point = unphysical value

          uncloud(iz,i,5)=uncloud(iz,i,5)+1.  
          phasecloud(iz,i)=7.           
       else

!! Curve for discrimination between ice & liquid particles
           perptmp1 = (pr2moy(iz,i)**5)*alpha50 + (pr2moy(iz,i)**4)*beta50 + (pr2moy(iz,i)**3)*gamma50 + (pr2moy(iz,i)**2)*delta50 + pr2moy(iz,i)*epsilon50 + zeta50

!! Curve for discrimination between ho particles & other
           perptmp2 = pr2moy(iz,i)*Ahoi + Bhoi
   
            if( (perpmoy(iz,i)-perptmp1).ge.0. )then
                 if (tempmoy(iz,i).gt.0 ) then
                    uncloud(iz,i,3)=uncloud(iz,i,3)+1.    !fake ice particles
                    watercloud(iz,i)=watercloud(iz,i)+1.  !fake ice particles   
                    phasecloud(iz,i)=5.
                 else
                    icecloud(iz,i)=icecloud(iz,i)+1.
                    phasecloud(iz,i)=2.
                 endif

            elseif( (perpmoy(iz,i)-perptmp2).ge.0. )then
                    if (tempmoy(iz,i).gt.-42 ) then
                       watercloud(iz,i)=watercloud(iz,i)+1.
                       phasecloud(iz,i)=1.
                    else
                       icecloud(iz,i)=icecloud(iz,i)+1.      ! fakeliq particles
                       uncloud(iz,i,2)=uncloud(iz,i,2)+1.    !fake liq particles  
                       phasecloud(iz,i)=4.                   !fake liq particles
                    endif
                 else
                   uncloud(iz,i,4)=uncloud(iz,i,4)+1.  !! horizontally oriented particles
                   phasecloud(iz,i)=6.                 !! horizontally oriented particles
            endif

       toplvlsat2=0
    if( srmoy(iz,i).gt.30)then
       toplvlsat2=iz

       if(toplvlsat2.gt.altmax-2)then
          print *, toplvlsat2
       endif

       goto 99 
    endif

   endif
  endif
 enddo 

99 continue

if(toplvlsat2.ne.0)then

! Other level below cloud with SR>30 phase = undefined
do iz=toplvlsat2-1,1,-1
    if(cloudfraction(iz,i).gt.0.)then
      uncloud(iz,i,1)=uncloud(iz,i,1)+1.
      phasecloud(iz,i)=3.  !! undefined cloud
    endif
enddo
       toplvlsat2=0
endif




  do  iz=1,altmax  
    if(cloudfraction(iz,i).gt.0.)then
icewaterres=watercloud(iz,i)+icecloud(iz,i)+uncloud(iz,i,1)+uncloud(iz,i,4)+uncloud(iz,i,5)
if((icewaterres.gt.1.).or.(icewaterres.eq.0.))then 
print *, 'error sum phase=',icewaterres
print *, watercloud(iz,i),icecloud(iz,i),uncloud(iz,i,:)
stop
endif
     
endif
  enddo



altend=0
altstart=0

if(nol.eq.1)then
!print *, "Non Over Lap MODE"

!!! NON OVERLAP MODE
!!! select only clouds in the highest isccp layer
!!! if isccpmid=1 then isccphigh must be 0
!!! A isccp layer can not have a cloud above 
B34: do  iz=altmax,1,-1  
     if(cloudfraction(iz,i).gt.0.)then


! Search the isscp high cloud fraction
         if (iz.ge.topmidl) then
     altstart=altmax
     altend=topmidl
            exit B34
         endif

! Search the isscp mid cloud fraction
         if ((iz.ge.toplowl).and.(iz.lt.topmidl)) then
         altstart=topmidl-1
         altend=toplowl
            exit B34
         endif

! Search the isscp low cloud fraction
        if (iz.lt.toplowl) then
            altstart=toplowl-1
            altend=1
            exit B34
       endif

    endif

enddo B34

else

!print *, "Over Lap MODE"

altend=1
altstart=altmax

endif

if(altend.ne.0)then   ! to avoid clear sky profil error
indbase=0
 do  iz=altend,altstart  

   if(cloudfraction(iz,i).gt.0.)then

! Calculation of cloudfraction as a function of the temperature
  do itemp=1,tempmax-1
    if ( (tempmoy(iz,i).gt.tempmod(itemp)).and.     &
       (tempmoy(iz,i) .le.tempmod(itemp+1)) )then

       cftemp(itemp,i)=cftemp(itemp,i)+1   ! cf tot

       if( watercloud(iz,i).ne.0 )then   
          cftempliq(itemp,i)=cftempliq(itemp,i)+1  ! cf liquid
       elseif( icecloud(iz,i).ne.0 )then  
          cftempice(itemp,i)=cftempice(itemp,i)+1  ! cf ice
       endif

    endif
  enddo


! Indice de depolarisation (indices: 0=clear-sky, 1=liquid only; 2=mixte; 3=ice only)
       
       depoltmp=0
       icewaterres=0


! Search the isscp low cloud fraction
   if (iz.lt.toplowl) then
      isccplow(i)=1

      if( icecloud(iz,i).gt.0 )then
         isccpice(1,i)=1 
      endif

      if( watercloud(iz,i).gt.0 )then
         isccpliq(1,i)=1
      endif
      
      do icat=1,catmax
         if(uncloud(iz,i,icat).gt.0)then
            isccpun(1,i,icat)=1
         endif
      enddo

   endif

! Search the isscp mid cloud fraction
   if ((iz.ge.toplowl).and.(iz.lt.topmidl)) then
      isccpmid(i)=1
            
      if(icecloud(iz,i).gt.0)then
         isccpice(2,i)=1 
      endif

      if(watercloud(iz,i).gt.0)then
         isccpliq(2,i)=1
      endif

      do icat=1,catmax
         if(uncloud(iz,i,icat).gt.0)then
            isccpun(2,i,icat)=1
         endif
      enddo

   endif

! Search the isscp high cloud fraction
   if (iz.ge.topmidl) then
      isccphigh(i)=1

height(i)=altmid(iz)+0.24

indbase=indbase+1
if(indbase.eq.1)then
  height2(i)=altmid(iz)-0.24
endif

      if(icecloud(iz,i).gt.0)then
         isccpice(3,i)=1 
      endif
      if(watercloud(iz,i).gt.0)then
         isccpliq(3,i)=1
      endif

      do icat=1,catmax
         if(uncloud(iz,i,icat).gt.0)then
            isccpun(3,i,icat)=1
         endif
      enddo

   endif

! Search cloud fraction on the column
            colcloud(i)=1

      if(icecloud(iz,i).gt.0)then
         isccpice(4,i)=1 
      endif

      if(watercloud(iz,i).gt.0)then
         isccpliq(4,i)=1
      endif

      do icat=1,catmax
         if(uncloud(iz,i,icat).gt.0)then
            isccpun(4,i,icat)=1
         endif
      enddo

    endif  !! endif cloudfraction > 0

         if((nanfraction(iz,i).eq.1.).or.(rejfraction(iz,i).eq.1.).or.(sefraction(iz,i).eq.1.))then
            nanprof=nanprof+1
         endif



  enddo
endif 


!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                  DAILY AVERAGE OF CLOUD DIAGNOSTICS                      *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


       
!******************* READING THE LAT LON MODEL GRID FILES *******************!

!---------------------------- lat lon grid file -----------------------------!
!                                                                            !
!   There is 3 grid in this program :                                        !
!                                                                            !
!       -CFMIP : 1° x 1° x 41levels each 480m from 0 to 19.2km               !
!                (361,181,41)                                                !
!       -LMDZ  : 3.75° x 2.53° x 19levels from 0 to 40.4823km                !
!                (96,72,19)                                                  !
!       -NASA  : 5° x 5° x 41levels each 480m from 0 to 19.2km               !
!               (73,37,41)                                                   !
!                                                                            !
!----------------------------------------------------------------------------!


!************************** BEGIN THE DAILY AVERAGE *************************!

!---------------------------- lat lon grid file -----------------------------!
!                                                                            !
!   increase the daily variables over the lat/lon (& alt & diag) grid        !
!                                                                            !
!----------------------------------------------------------------------------!

! First level always equal to 0 because the average is perfomed between the 
! level i and i-1.
!print *, "toto3"


  do ilon=1,lonmax-1  !longitude
     if( (lon(i).ge.lonmod(ilon)) .and. (lon(i).lt.lonmod(ilon+1)) )then
         
    do ilat=1,latmax-1  !latitude
        if ( (lat(i).ge.latmod(ilat)) .and. (lat(i).lt.latmod(ilat+1)) )then

     
             indtot(ilon,ilat)=indtot(ilon,ilat)+1
             
             nanprof=0
             do ialt=1,altmax
                if((nanfraction(ialt,i).eq.1.).or.(sefraction(ialt,i).eq.1.) &
                   .or.(rejfraction(ialt,i).eq.1.))then
                   nanprof=nanprof+1
                endif
             enddo
             
             nansat=0  
             nanlow=0
             do ialt=1,toplowl-1
                if(satfraction(ialt,i).eq.1.)then
                   nansat=nansat+1       ! fully att layer
                endif
                if(sefraction(ialt,i).eq.1.)then
                   nanlow=nanlow+1       ! under SE layer
                endif
             enddo
         
             nanmid=0
             do ialt=toplowl,topmidl-1
                if(satfraction(ialt,i).eq.1.)then
                   nanmid=nanmid+1
                endif
             enddo

          ! Cloudtop Height average (GEWEX option)
          heightday(ilat,ilon,jour)=heightday(ilat,ilon,jour)+height(i)
          heightday2(ilat,ilon,jour)=heightday2(ilat,ilon,jour)+height2(i)
          if( height(i).gt.0)then
             indheight(ilat,ilon,jour)=indheight(ilat,ilon,jour)+1
          endif

do nphase = 1,4
isccpliqday(ilat,ilon,jour,nphase)=isccpliqday(ilat,ilon,jour,nphase)+isccpliq(nphase,i)
isccpiceday(ilat,ilon,jour,nphase)=isccpiceday(ilat,ilon,jour,nphase)+isccpice(nphase,i)
  do icat=1,catmax
isccpunday(ilat,ilon,jour,nphase,icat)=isccpunday(ilat,ilon,jour,nphase,icat)+isccpun(nphase,i,icat)
  enddo
enddo

!!!! A RAJOUTER LES ISCCPblablaDAY
!!!! ENREGISTRER LES 2 TYPE DE ISCCP en occurrence ?

         ! colclearday(ilat,ilon,jour)=colclearday(ilat,ilon,jour)+colclear(i)
          
          ! select only the profil .ne. NaN & Low layer not fully attenuated
          if( (nanprof.le.(altmax-4)).and.(nansat.lt.(toplowl-1)).and. &
             (nanlow.lt.(toplowl-1)) )then
          isccplowday(ilat,ilon,jour)=isccplowday(ilat,ilon,jour)+isccplow(i)
          isccpinddaylow(ilat,ilon,jour)=isccpinddaylow(ilat,ilon,jour)+1
          endif

          ! select only the profil .ne. NaN & Mid layer not fully attenuated
          if((nanprof.le.(altmax-4)).and.(nanmid.lt.(topmidl-toplowl)))then
          isccpinddaymid(ilat,ilon,jour)=isccpinddaymid(ilat,ilon,jour)+1
          isccpmidday(ilat,ilon,jour)=isccpmidday(ilat,ilon,jour)+isccpmid(i)
          endif

          ! select only the profil .ne. NaN
          if(nanprof.le.(altmax-4))then
          isccpindday(ilat,ilon,jour)=isccpindday(ilat,ilon,jour)+1  
          isccphighday(ilat,ilon,jour)=isccphighday(ilat,ilon,jour)+         &
                                       isccphigh(i)
          colcloudday(ilat,ilon,jour)=colcloudday(ilat,ilon,jour)+colcloud(i)    
          endif

          if( (isccplowday(ilat,ilon,jour)/isccpinddaylow(ilat,ilon,jour)) .gt. 1.00001)then
             print *, file6
             print *, i,lon(i),lat(i),nansat
             print *, isccplowday(ilat,ilon,jour),isccpinddaylow(ilat,ilon,jour)
          endif

 
 do  iz=1,altmax  
     if(cloudfraction(iz,i).gt.0.)then

if(indicetemp(iz,i).gt.0)then
       coltemp(ilat,ilon,jour)=coltemp(ilat,ilon,jour)+(tempmoy(iz,i)/indicetemp(iz,i))
       indcoltemp(ilat,ilon,jour)=indcoltemp(ilat,ilon,jour)+1

! Search the isscp low cloud fraction
        if (iz.lt.toplowl) then
            lowtemp(ilat,ilon,jour)=lowtemp(ilat,ilon,jour)+(tempmoy(iz,i)/indicetemp(iz,i))
            indlowtemp(ilat,ilon,jour)=indlowtemp(ilat,ilon,jour)+1
       endif

! Search the isscp mid cloud fraction
         if ((iz.ge.toplowl).and.(iz.lt.topmidl)) then
            midtemp(ilat,ilon,jour)=midtemp(ilat,ilon,jour)+(tempmoy(iz,i)/indicetemp(iz,i))
            indmidtemp(ilat,ilon,jour)=indmidtemp(ilat,ilon,jour)+1
         endif

! Search the isscp high cloud fraction
         if (iz.ge.topmidl) then
             hightemp(ilat,ilon,jour)=hightemp(ilat,ilon,jour)+(tempmoy(iz,i)/indicetemp(iz,i))
            indhightemp(ilat,ilon,jour)=indhightemp(ilat,ilon,jour)+1
         endif
endif

      endif
 enddo


! cloudfraction day as a function of temperature
  do itemp=1,tempmax-1
     cftempday(ilat,ilon,itemp,jour)=cftempday(ilat,ilon,itemp,jour)+cftemp(itemp,i)
     cftempliqday(ilat,ilon,itemp,jour)=cftempliqday(ilat,ilon,itemp,jour)+   &
                                   cftempliq(itemp,i)
     cftempiceday(ilat,ilon,itemp,jour)=cftempiceday(ilat,ilon,itemp,jour)+   &
                                   cftempice(itemp,i)
  enddo


            do ialt=altmax,1,-1  
            indnan(ilat,ilon,ialt)=indnan(ilat,ilon,ialt)+1
            cloudfractday(ilat,ilon,ialt,jour)=                              &
            cloudfractday(ilat,ilon,ialt,jour)+cloudfraction(ialt,i)
            clearfractday(ilat,ilon,ialt,jour)=                              &
            clearfractday(ilat,ilon,ialt,jour)+clearfraction(ialt,i)
            uncertfractday(ilat,ilon,ialt,jour)=                             &
            uncertfractday(ilat,ilon,ialt,jour)+uncertfraction(ialt,i)
            icecloudfractday(ilat,ilon,ialt,jour)=                            &
            icecloudfractday(ilat,ilon,ialt,jour)+icecloud(ialt,i)
            watercloudfractday(ilat,ilon,ialt,jour)=                          &
            watercloudfractday(ilat,ilon,ialt,jour)+watercloud(ialt,i)

            do icat=1,catmax
            uncloudfractday(ilat,ilon,ialt,jour,icat)=                       &
            uncloudfractday(ilat,ilon,ialt,jour,icat)+uncloud(ialt,i,icat)
            enddo


           
            if((nanfraction(ialt,i).ne.1).and.(sefraction(ialt,i).ne.1).and.(satfraction(ialt,i).ne.1).and.(rejfraction(ialt,i).ne.1))then
            indday(ilat,ilon,ialt,jour)=indday(ilat,ilon,ialt,jour)+1

  do itemp=1,tempmax-1
               if ( (tempmoy(ialt,i).gt.tempmod(itemp)).and.     &
                   (tempmoy(ialt,i) .le.tempmod(itemp+1)) )then
                   indcftemp(ilat,ilon,itemp,jour)=indcftemp(ilat,ilon,itemp,jour)+1
                endif
  enddo
            ! indice of cloudfraction day as a function of temperature

            if(indday(ilat,ilon,ialt,jour).lt.sum(uncloudfractday(ilat,ilon,ialt,jour,:)))then
               print *, 'error indice < un_cloud',i,ialt
               stop
            endif
        


            endif

            if(srmoy(ialt,i).eq.srmod(1))then
                diagSR(ilon,ilat,ialt,1)=diagSR(ilon,ilat,ialt,1)+1                
            endif
            if(srmoy(ialt,i).eq.srmod(2))then
                diagSR(ilon,ilat,ialt,2)=diagSR(ilon,ilat,ialt,2)+1
            endif

  !       if( (srmoy(ialt,i).ge.srmod(1)) .and. (srmoy(ialt,i).lt.srmod(diagmax)) )then
                  do idiag=3,diagmax-1
                     if ( (srmoy(ialt,i).ge.srmod(idiag)).and.             &
                        (srmoy(ialt,i).lt.srmod(idiag+1)) )then
                        diagSR(ilon,ilat,ialt,idiag)=                   &
                        diagSR(ilon,ilat,ialt,idiag)+1
                     endif
                  enddo

                  do idiag=8,diagmax-1
                     if ( (srmoy(ialt,i).ge.srmod(idiag)).and.             &
                        (srmoy(ialt,i).lt.srmod(idiag+1)) )then
                       if(icecloud(ialt,i).ne.0)then
                        diagSRpha(ilon,ilat,ialt,idiag-7,2)=                   &
                        diagSRpha(ilon,ilat,ialt,idiag-7,2)+1

                        elseif(watercloud(ialt,i).ne.0)then
                        diagSRpha(ilon,ilat,ialt,idiag-7,1)=                   &
                        diagSRpha(ilon,ilat,ialt,idiag-7,1)+1

                        elseif(uncloud(ialt,i,1).ne.0)then
                        diagSRpha(ilon,ilat,ialt,idiag-7,3)=                   &
                        diagSRpha(ilon,ilat,ialt,idiag-7,3)+1
   
                        endif
                     endif
                  enddo




              enddo
            endif
         enddo
      endif
   enddo



enddo !!!!!!!!!! END IT LOOP

!print *, "toto4"

!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                   RECORD INSTANTANEOUS PHASE FILE                        *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!
if('instant_switch'.eq.'on')then
call instant_phase(file4,altmax,it,phasecloud)
endif


do ilon=1,lonmax-1  !latitude
        do ilat=1,latmax-1  !longitude
   !         print *, indtotmean(ilon,ilat)
          if(indtot(ilon,ilat).gt.0)then
           indtotmean(ilon,ilat)=indtotmean(ilon,ilat)+1
          endif
        enddo
enddo

indtot(:,:)=0;

!CASE DEFAULT
!print *, "The model you entered doesn't match, try another" 
!
!ENDSELECT


!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                           WRITE OUTPUT FILES                             *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

!----------------------------- Save output file -----------------------------!
!                                                                            !
! Sauve = "wrf" when model WRF is used, then the data are saved in ASCII       !
! format, and for every hdf calipso file. Only the matching WRF grid values  !
! are saved.                                                                 !
!     - Lat : 199values between 15° to 55°                                   !
!     - Lon : 179values between -18° to 36°                                  !
!     - alt : 48levels from 0 to 35.5km                                      !
!                                                                            !
! Sauve = "chim" when model CHIMERE is used, then the data are saved in ASCII  !
! format, and for every hdf calipso file. Only the matching CHIMERE grid     !
! values are saved.                                                          !
!     - Lat : 71values between -5° to 65° (every 1°)                         !
!     - Lon : 181values between -100° to 80° (every 1°)                      !
!     - alt : 281levels from 0 to 14km (every 50m)                           !
!                                                                            !          
! Sauve = "lmdz" when LMDZ is used, the data are saved in a single netcdf file !
! when the last hdf calipso file of the period have been read                !        
!                                                                            !
!----------------------------------------------------------------------------!
622 continue

sauve=model

SELECT CASE (sauve)


!**************************** WRF OUTPUT FORMAT *****************************! 
CASE ("wrf")

comptpf=0
do i=1,it
     if ((lat(i).gt.15.0).and.(lat(i).lt.55.0))then
           if((lon(i).gt.-18.0).and.(lon(i).lt.36.0))then
              comptpf=comptpf+1
           endif
      endif
enddo

allocate(latwrf(comptpf),lonwrf(comptpf),SRwrf(altmax,comptpf),DEPOLwrf(altmax,comptpf),CRwrf(altmax,comptpf),SEwrf(comptpf),temps2wrf(comptpf))

if(comptpf.eq.0)then
   goto 624  !track out from the domain
else

comptpf=0
do i=1,it
     if ((lat(i).gt.15.0).and.(lat(i).lt.55.0))then
           if((lon(i).gt.-18.0).and.(lon(i).lt.36.0))then
              comptpf=comptpf+1
              SEwrf(comptpf)=SE(i)
              temps2wrf(comptpf)=temps2(i)
              latwrf(comptpf)=lat(i)
              lonwrf(comptpf)=lon(i)
              do j=1,altmax
                 SRwrf(j,comptpf)=srmoy(j,i)
                 CRwrf(j,comptpf)=crmoy(j,i)
                 DEPOLwrf(j,comptpf)=depolmoy(j,i)
              enddo
           endif
      endif
enddo


  write(numfichc,'(i4)')numfich
  write(datec,'(I5)')date

! Name of instantaneous file, one by hdf calipso file  
  file4='SR_CR_DEPOL_200'//datec//'_'//trim(switch)//'_'//trim(gcm)//'_'//    &
  trim(ADJUSTL(numfichc))//'.nc'


  call SR_CR_DR_2nc(file4,altmid,altmod,resd,altmax,switch,gcm,comptpf,latwrf,lonwrf,SEwrf,temps2wrf,&
                  SRwrf,CRwrf,DEPOLwrf)

endif

624 continue



! Deallocate SDS/META variables
deallocate(atb,stat = OK_buffer)
deallocate(atb2,stat = OK_buffer)
deallocate(perp,stat = OK_buffer)
deallocate(mol,stat = OK_buffer)
deallocate(pres,stat = OK_buffer) 
deallocate(SE,SEwrf,stat = OK_buffer)
deallocate(mol2,stat = OK_buffer)
deallocate(mol3,stat = OK_buffer)
deallocate(mol4)
deallocate(temp,stat = OK_buffer)
deallocate(temp2,stat = OK_buffer)

if(alt_pres=='pressure')then
   deallocate(pres2,stat = OK_buffer)
endif
deallocate(altl,stat = OK_buffer)
deallocate(altm,stat = OK_buffer)
deallocate(temps2,temps2wrf,stat = OK_buffer)
deallocate(lonwrf,latwrf,SRwrf,CRwrf,DEPOLwrf)


!continue 


!************************** CHIMERE OUTPUT FORMAT ***************************!
CASE ("chim")

print *, 'Regrid of data Done'
print *, ''
print *, 'Save the data'

do k=1,lonmax
   do j=1,latmax 
      if(mheure(j,k).ne.0)then
      box=box+1
      if(lonmod(k)==12)then
       !  print *, k
endif

      write(11,'(4(2x,I6),2(2x,F10.2),(2x,E13.6))'),numfich,box,month,day,   &
      lonmod(k),latmod(j),mheure(j,k)/indiceh(j,k)
         do iz=1,altmax
      !       write(11,103),altmod(iz),srmoy(j,k,iz),crmoy(j,k,iz),depolmoy(j,k,iz),tempmoy(j,k,iz) &
       !                   /indicet(j,k,iz),molmoy(j,k,iz)/indicem(j,k,iz)
            enddo
      endif
      
   enddo
enddo


!**************************** LMDZ OUTPUT FORMAT ****************************!
CASE ("lmdz")
goto 666

CASE DEFAULT
print *, "error" 

ENDSELECT

666 continue
647 continue


!***** DEALLOCATE SDS OUTPUT AND INSTANT LMDZ OUTPUT VAR AND CLOSE FILES ****!

deallocate(pr2moy,molmoy, stat = OK_buffer)!
deallocate(srmoy)
deallocate(lat,stat = OK_buffer)
deallocate(lon,stat = OK_buffer)
!deallocate(indtot,indret)
!deallocate(indretlow,indretmid,indrethigh)
deallocate(indice,indicem,indice2)
deallocate(indicep,indicep2,depolmoy,parmoy,perpmoy,pr2moy2,crmoy)
deallocate(tempmoy,indicetemp)

if(model=='wrf')then
 deallocate(indice,indicem,indice2)
deallocate(indicep,indicep2,crmoy,depolmoy,parmoy,perpmoy,pr2moy2)
! deallocate(altmod,stat=OK_buffer)
endif


if(model=='chimere')then
   deallocate(latmod,lonmod,altmod,stat = OK_buffer)
   deallocate(mheure,stat = OK_buffer)
   deallocate(indiceh,stat = OK_buffer)
   deallocate(indice,indicem,indice2,stat = OK_buffer)
endif

if(model=='lmdz')then
! deallocate(indice,indicem)!,indice2,indicep,indicep2
if(alt_pres=='pressure')then
 deallocate(SEp, stat = OK_buffer)
endif
 deallocate(isccplow,isccpmid,isccphigh,colcloud)!,colclear, stat = OK_buffer)
 deallocate(cloudfraction,clearfraction,satfraction,uncertfraction,          &
            nanfraction,sefraction,rejfraction)!,fractot)
 deallocate(icecloud,watercloud,uncloud,phasecloud)
 deallocate(isccpliq,isccpice,isccpun)
 deallocate(cftemp,cftempliq,cftempice)
 deallocate(height,height2)
endif

print *, 'Deallocate buffers done'


887 continue
!print *, 'rm the file'
!command2='rm -f /tmp/'//trim(filetmp)
!call system(command2)
print *, ''

print *, 'go to the next file'


  ! Temps elapsed final
  call system_clock(count=t2, count_rate=ir)
  tempstot=real(t2 - t1,kind=4)/real(ir,kind=4)

  ! Temps CPU de calcul final
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

         
 print '(//,3X,"Temps elapsed       : ",1PE10.3," sec.",/, &
           & 3X,"Temps CPU           : ",1PE10.3," sec.",/,//)', &
           tempstot,t_cpu




goto 888          ! next Calipso file


999 continue 
 
! empty file list checking
if(file2(1:1).gt.' ')then
continue
else
print *, 'hdf file list empty'
print *, file2(1:1)
stop
endif



!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                  MONTHLY AVERAGE CLOUD DIAGNOSTICS                       *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

!----------------------------- monthly average ------------------------------!
!                                                                            !
!  If there is no more hdf CALIPSO file to read, then the monthly average    !
!  can begin. This module is operating only in the lmdz case, that is to say !
!  when we have to save the monthly diagnostic.                              !
!  In the other case (chimere & wrf), the program ends here.                 !
!                                                                            !
!----------------------------------------------------------------------------!

SELECT CASE (model)
CASE ("lmdz")


!****************************************************************************!
!*!!!!!!!!!!!!!! PART I : CLOUDY LOW MID HIGH MAP FILES !!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


! Allocation / initialization of isscp monthly variables
print *, 'allocation / initialization of isccp monthly variables'

   allocate(monthheight(lonmax-1,latmax-1),indmonthheight(lonmax-1,latmax-1),      &
            monthheight2(lonmax-1,latmax-1))
   allocate(monthisccplow(lonmax-1,latmax-1),monthisccpmid(lonmax-1,latmax-1),       &
            monthisccphigh(lonmax-1,latmax-1))
   allocate(monthcolcloud(lonmax-1,latmax-1), monthcolclear(lonmax-1,latmax-1),       &
            isccpdaypermonthlow(lonmax-1,latmax-1),                                  &
            isccpdaypermonthmid(lonmax-1,latmax-1),                                  &           
            isccpdaypermonth(lonmax-1,latmax-1))
   allocate(hlow(lonmax-1,latmax-1,histmax-1),hmid(lonmax-1,latmax-1,histmax-1), &
            hhigh(lonmax-1,latmax-1,histmax-1),hcol(lonmax-1,latmax-1,histmax-1),&
            hheight(lonmax-1,latmax-1,histmax2-1))

   allocate(hlowtemp(lonmax-1,latmax-1,histtempmax2),hmidtemp(lonmax-1,latmax-1,histtempmax2), &
            hhightemp(lonmax-1,latmax-1,histtempmax2),hcoltemp(lonmax-1,latmax-1,histtempmax2))
   allocate(monthlowtemp(lonmax-1,latmax-1),indmonthlowtemp(lonmax-1,latmax-1), &
            monthmidtemp(lonmax-1,latmax-1),indmonthmidtemp(lonmax-1,latmax-1), &
            monthhightemp(lonmax-1,latmax-1),indmonthhightemp(lonmax-1,latmax-1), &
            monthcoltemp(lonmax-1,latmax-1),indmonthcoltemp(lonmax-1,latmax-1))


allocate(monthisccpliq(lonmax-1,latmax-1,4),monthisccpice(lonmax-1,latmax-1,4))
allocate(monthisccpun(lonmax-1,latmax-1,4,catmax))
allocate(monthisccpphase(lonmax-1,latmax-1,4))
allocate(indmonthphase(lonmax-1,latmax-1,4),indmonthphase2(lonmax-1,latmax-1,4))
allocate(inddayphase(latmax-1,lonmax-1,daymax,4))

inddayphase(:,:,:,:)=0;
monthisccpliq(:,:,:)=0; monthisccpice(:,:,:)=0; 
monthisccpun(:,:,:,:)=0; 
monthisccpphase(:,:,:)=0;indmonthphase(:,:,:)=0;indmonthphase2(:,:,:)=0

histmod(:)=0; histmod2(:)=0;
hlow(:,:,:)=0;hmid(:,:,:)=0;hhigh(:,:,:)=0;hcol(:,:,:)=0;hheight(:,:,:)=0;
monthheight(:,:)=0; monthheight2(:,:)=0; indmonthheight(:,:)=0;
monthisccplow(:,:)=0;monthisccpmid(:,:)=0;monthisccphigh(:,:)=0;
monthcolcloud(:,:)=0;isccpdaypermonth(:,:)=0;monthcolclear(:,:)=0
isccpdaypermonthlow(:,:)=0;isccpdaypermonthmid(:,:)=0;
histtempmod(:)=0; histtempmod2(:)=0;

hlowtemp(:,:,:)=0;monthlowtemp(:,:)=0;indmonthlowtemp(:,:)=0;
hmidtemp(:,:,:)=0;monthmidtemp(:,:)=0;indmonthmidtemp(:,:)=0;
hhightemp(:,:,:)=0;monthhightemp(:,:)=0;indmonthhightemp(:,:)=0;
hcoltemp(:,:,:)=0;monthcoltemp(:,:)=0;indmonthcoltemp(:,:)=0;

histtempmod(1)=150
histtempmod(2)=180
histtempmod(histtempmax)=320

do ihisttemp=2,histtempmax-2
   histtempmod(ihisttemp+1)=histtempmod(ihisttemp)+5
enddo


do ihisttemp=1,histtempmax2
   histtempmod2(ihisttemp)=(histtempmod(ihisttemp)+histtempmod(ihisttemp+1))/2
enddo

!print *, 'histtempmod',histtempmod
!print *, 'histtempmod2',histtempmod2


do ihist=1,histmax-2
   histmod(ihist+1)=histmod(ihist)+0.1
enddo
 
histmod(histmax)=1.01

do ihist=1,histmax2-1
   histmod2(ihist+1)=histmod2(ihist)+1
enddo




!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!
! changement à rajouter !
! normaliser le calcul de la phase par jour


!********* CALCULATION OF DAILY DIAGNOSTIC WITH MATCHING INDEXES ************!
! moyenne journalière effectuée avant la moyenne mensuelle pour donner le même
! poid à tous les jours malgré la quantité variable de donnée d'un jour à 
! l'autre.

      do jour=1,31 
         do ilon=1,lonmax-1
           do ilat=1,latmax-1  

             if ( indlowtemp(ilat,ilon,jour).gt.0 ) then
           lowtemp(ilat,ilon,jour)=                                         &
           (lowtemp(ilat,ilon,jour)/indlowtemp(ilat,ilon,jour))+273.15
             endif
            if ( indmidtemp(ilat,ilon,jour).gt.0 ) then
           midtemp(ilat,ilon,jour)=                                         &
           (midtemp(ilat,ilon,jour)/indmidtemp(ilat,ilon,jour))+273.15
             endif
            if ( indhightemp(ilat,ilon,jour).gt.0 ) then
           hightemp(ilat,ilon,jour)=                                         &
           (hightemp(ilat,ilon,jour)/indhightemp(ilat,ilon,jour))+273.15
             endif
            if ( indcoltemp(ilat,ilon,jour).gt.0 ) then
           coltemp(ilat,ilon,jour)=                                         &
           (coltemp(ilat,ilon,jour)/indcoltemp(ilat,ilon,jour))+273.15
             endif


             if ( isccpindday(ilat,ilon,jour).gt.0 ) then
           isccphighday(ilat,ilon,jour)=                                     &
           isccphighday(ilat,ilon,jour)/isccpindday(ilat,ilon,jour)
           colcloudday(ilat,ilon,jour)=                                      &
           colcloudday(ilat,ilon,jour)/isccpindday(ilat,ilon,jour)

if (isccpiceday(ilat,ilon,jour,3)+isccpliqday(ilat,ilon,jour,3)+         &
    +sum(isccpunday(ilat,ilon,jour,3,4:5))+isccpunday(ilat,ilon,jour,3,1).gt.0)then
isccptemp=isccpiceday(ilat,ilon,jour,3)+isccpliqday(ilat,ilon,jour,3)+          &
          isccpunday(ilat,ilon,jour,3,1)+ sum(isccpunday(ilat,ilon,jour,3,4:5))

          isccpliqday(ilat,ilon,jour,3)= isccpliqday(ilat,ilon,jour,3)/            &
          isccptemp * isccphighday(ilat,ilon,jour)
          isccpiceday(ilat,ilon,jour,3)= isccpiceday(ilat,ilon,jour,3)/            &
          isccptemp * isccphighday(ilat,ilon,jour)
do icat=1,catmax
          isccpunday(ilat,ilon,jour,3,icat)= isccpunday(ilat,ilon,jour,3,icat)/            &
          isccptemp * isccphighday(ilat,ilon,jour)
enddo

endif 

if (isccpiceday(ilat,ilon,jour,4)+isccpliqday(ilat,ilon,jour,4)+         &
    +sum(isccpunday(ilat,ilon,jour,4,4:5))+isccpunday(ilat,ilon,jour,4,1).gt.0)then
isccptemp=isccpiceday(ilat,ilon,jour,4)+isccpliqday(ilat,ilon,jour,4)+          &
          isccpunday(ilat,ilon,jour,4,1)+ sum(isccpunday(ilat,ilon,jour,4,4:5))

          isccpliqday(ilat,ilon,jour,4)= isccpliqday(ilat,ilon,jour,4)/            &
          isccptemp * colcloudday(ilat,ilon,jour)
          isccpiceday(ilat,ilon,jour,4)= isccpiceday(ilat,ilon,jour,4)/            &
          isccptemp * colcloudday(ilat,ilon,jour)

do icat=1,catmax
          isccpunday(ilat,ilon,jour,4,icat)= isccpunday(ilat,ilon,jour,4,icat)/      &
          isccptemp * isccphighday(ilat,ilon,jour)
enddo

endif

             endif

             if ( indheight(ilat,ilon,jour).gt.0 ) then
                heightday(ilat,ilon,jour)=                                   &
                heightday(ilat,ilon,jour)/indheight(ilat,ilon,jour)                      
                heightday2(ilat,ilon,jour)=                                   &
                heightday2(ilat,ilon,jour)/indheight(ilat,ilon,jour)                  
             endif

             if ( isccpinddaylow(ilat,ilon,jour).gt.0 ) then
           isccplowday(ilat,ilon,jour)=                                      &
           isccplowday(ilat,ilon,jour)/isccpinddaylow(ilat,ilon,jour)

if (isccpiceday(ilat,ilon,jour,1)+isccpliqday(ilat,ilon,jour,1)+         &
    +sum(isccpunday(ilat,ilon,jour,1,4:5))+isccpunday(ilat,ilon,jour,1,1).gt.0)then
isccptemp=isccpiceday(ilat,ilon,jour,1)+isccpliqday(ilat,ilon,jour,1)+          &
          isccpunday(ilat,ilon,jour,1,1)+ sum(isccpunday(ilat,ilon,jour,1,4:5))

          isccpliqday(ilat,ilon,jour,1)= isccpliqday(ilat,ilon,jour,1)/            &
          isccptemp * isccplowday(ilat,ilon,jour)
          isccpiceday(ilat,ilon,jour,1)= isccpiceday(ilat,ilon,jour,1)/            &
          isccptemp * isccplowday(ilat,ilon,jour)

do icat=1,catmax
          isccpunday(ilat,ilon,jour,1,icat)= isccpunday(ilat,ilon,jour,1,icat)/      &
          isccptemp * isccplowday(ilat,ilon,jour)
enddo

endif


             endif
             if ( isccpinddaymid(ilat,ilon,jour).gt.0 ) then
           isccpmidday(ilat,ilon,jour)=                                      &
           isccpmidday(ilat,ilon,jour)/isccpinddaymid(ilat,ilon,jour)

if (isccpiceday(ilat,ilon,jour,2)+isccpliqday(ilat,ilon,jour,2)+         &
    +sum(isccpunday(ilat,ilon,jour,2,4:5))+isccpunday(ilat,ilon,jour,2,1).gt.0)then
isccptemp=isccpiceday(ilat,ilon,jour,2)+isccpliqday(ilat,ilon,jour,2)+          &
          isccpunday(ilat,ilon,jour,2,1)+ sum(isccpunday(ilat,ilon,jour,2,4:5))

          isccpliqday(ilat,ilon,jour,2)= isccpliqday(ilat,ilon,jour,2)/            &
          isccptemp * isccpmidday(ilat,ilon,jour)
          isccpiceday(ilat,ilon,jour,2)= isccpiceday(ilat,ilon,jour,2)/            &
          isccptemp * isccpmidday(ilat,ilon,jour)

do icat=1,catmax
          isccpunday(ilat,ilon,jour,2,icat)= isccpunday(ilat,ilon,jour,2,icat)/     &
          isccptemp * isccpmidday(ilat,ilon,jour)
enddo 

endif



if ((isccpiceday(ilat,ilon,jour,2)+isccpliqday(ilat,ilon,jour,2)).gt.isccpmidday(ilat,ilon,jour))then
print *, isccpmidday(ilat,ilon,jour),isccpiceday(ilat,ilon,jour,2)+isccpliqday(ilat,ilon,jour,2)
endif

             endif
             if ( isccplowday(ilat,ilon,jour) .gt. colcloudday(ilat,ilon,jour) )then
                print *, ilon,ilat,jour
                print *, isccplowday(ilat,ilon,jour),colcloudday(ilat,ilon,jour)
             endif


do nphase=1,4
   isccptemp=isccpiceday(ilat,ilon,jour,nphase)+isccpliqday(ilat,ilon,jour,nphase)
   if(isccptemp.gt.0.)then
      inddayphase(ilat,ilon,jour,nphase)=inddayphase(ilat,ilon,jour,nphase)+1.
      isccpphaseday(ilat,ilon,jour,nphase)=isccpiceday(ilat,ilon,jour,nphase)/isccptemp
   endif
enddo

           enddo
         enddo
      enddo


!************* INCREMENTING THE MONTHLY DIAGNOSTIC AND INDEXES **************!

      do jour=1,31 
         do ilon=1,lonmax-1
           do ilat=1,latmax-1  

            if ( indlowtemp(ilat,ilon,jour).gt.0 ) then
            do ihisttemp=1,histtempmax-1
              if( (lowtemp(ilat,ilon,jour).ge.histtempmod(ihisttemp)) .and.  &
                (lowtemp(ilat,ilon,jour).lt.histtempmod(ihisttemp+1)) ) then
              hlowtemp(ilon,ilat,ihisttemp)=                                         &
              hlowtemp(ilon,ilat,ihisttemp)+1
              endif
            enddo
            monthlowtemp(ilon,ilat)=monthlowtemp(ilon,ilat)+                  &
                                     lowtemp(ilat,ilon,jour)
            indmonthlowtemp(ilon,ilat)=indmonthlowtemp(ilon,ilat)+1            
            endif

            if ( indmidtemp(ilat,ilon,jour).gt.0 ) then
            do ihisttemp=1,histtempmax-1
              if( (midtemp(ilat,ilon,jour).ge.histtempmod(ihisttemp)) .and.  &
                (midtemp(ilat,ilon,jour).lt.histtempmod(ihisttemp+1)) ) then
              hmidtemp(ilon,ilat,ihisttemp)=                                         &
              hmidtemp(ilon,ilat,ihisttemp)+1
              endif
            enddo
            monthmidtemp(ilon,ilat)=monthmidtemp(ilon,ilat)+                  &
                                     midtemp(ilat,ilon,jour)
            indmonthmidtemp(ilon,ilat)=indmonthmidtemp(ilon,ilat)+1            
            endif

            if ( indhightemp(ilat,ilon,jour).gt.0 ) then
            do ihisttemp=1,histtempmax-1
              if( (hightemp(ilat,ilon,jour).ge.histtempmod(ihisttemp)) .and.  &
                (hightemp(ilat,ilon,jour).lt.histtempmod(ihisttemp+1)) ) then
              hhightemp(ilon,ilat,ihisttemp)=                                         &
              hhightemp(ilon,ilat,ihisttemp)+1
              endif
            enddo
            monthhightemp(ilon,ilat)=monthhightemp(ilon,ilat)+                  &
                                     hightemp(ilat,ilon,jour)
            indmonthhightemp(ilon,ilat)=indmonthhightemp(ilon,ilat)+1            
            endif

            if ( indcoltemp(ilat,ilon,jour).gt.0 ) then
            do ihisttemp=1,histtempmax-1
              if( (coltemp(ilat,ilon,jour).ge.histtempmod(ihisttemp)) .and.  &
                (coltemp(ilat,ilon,jour).lt.histtempmod(ihisttemp+1)) ) then
              hcoltemp(ilon,ilat,ihisttemp)=                                         &
              hcoltemp(ilon,ilat,ihisttemp)+1
              endif
            enddo
            monthcoltemp(ilon,ilat)=monthcoltemp(ilon,ilat)+                  &
                                     coltemp(ilat,ilon,jour)
            indmonthcoltemp(ilon,ilat)=indmonthcoltemp(ilon,ilat)+1            
            endif


          if ( indheight(ilat,ilon,jour).gt.0 ) then
             do ihist=1,histmax2-1
                if( (heightday(ilat,ilon,jour).ge.histmod2(ihist)) .and.&
                    (heightday(ilat,ilon,jour).lt.histmod2(ihist+1)) ) then
                    hheight(ilon,ilat,ihist)=hheight(ilon,ilat,ihist)+1
                 endif
              enddo
             monthheight(ilon,ilat)=monthheight(ilon,ilat)+                  &
                                    heightday(ilat,ilon,jour)
             monthheight2(ilon,ilat)=monthheight2(ilon,ilat)+                  &
                                    heightday2(ilat,ilon,jour)
             indmonthheight(ilon,ilat)=indmonthheight(ilon,ilat)+1
          endif

          if ( isccpindday(ilat,ilon,jour).gt.0 ) then

              do ihist=1,histmax-1
                 if( (isccphighday(ilat,ilon,jour).ge.histmod(ihist)) .and.  &
                     (isccphighday(ilat,ilon,jour).lt.histmod(ihist+1)) ) then
                     hhigh(ilon,ilat,ihist)=hhigh(ilon,ilat,ihist)+1
                 endif

                 if( (colcloudday(ilat,ilon,jour).ge.histmod(ihist)) .and.  &
                     (colcloudday(ilat,ilon,jour).lt.histmod(ihist+1)) )then
                     hcol(ilon,ilat,ihist)=hcol(ilon,ilat,ihist)+1
                  endif
               enddo
 
          monthisccphigh(ilon,ilat)=monthisccphigh(ilon,ilat)+               &
                                   isccphighday(ilat,ilon,jour)


          monthisccpliq(ilon,ilat,3)=monthisccpliq(ilon,ilat,3)+             &
                                     isccpliqday(ilat,ilon,jour,3)      
          monthisccpliq(ilon,ilat,4)=monthisccpliq(ilon,ilat,4)+             &
                                     isccpliqday(ilat,ilon,jour,4)

          monthisccpice(ilon,ilat,3)=monthisccpice(ilon,ilat,3)+             &
                                     isccpiceday(ilat,ilon,jour,3)      
          monthisccpice(ilon,ilat,4)=monthisccpice(ilon,ilat,4)+             &
                                     isccpiceday(ilat,ilon,jour,4)
      
do icat=1,catmax
          monthisccpun(ilon,ilat,3,icat)=monthisccpun(ilon,ilat,3,icat)+             &
                                     isccpunday(ilat,ilon,jour,3,icat)      
          monthisccpun(ilon,ilat,4,icat)=monthisccpun(ilon,ilat,4,icat)+             &
                                     isccpunday(ilat,ilon,jour,4,icat)
enddo

          monthcolcloud(ilon,ilat)=monthcolcloud(ilon,ilat)+                 &
                                   colcloudday(ilat,ilon,jour)
       !   monthcolclear(ilon,ilat)=monthcolclear(ilon,ilat)+                 &
       !                            colclearday(ilat,ilon,jour)

!isccpdaypermonth(ilon,ilat)=isccpdaypermonth(ilon,ilat)+isccpindday(ilat,ilon,jour)
          isccpdaypermonth(ilon,ilat)=isccpdaypermonth(ilon,ilat)+1
          endif

          if ( isccpinddaylow(ilat,ilon,jour).gt.0 ) then
             do ihist=1,histmax-1
                if( (isccplowday(ilat,ilon,jour).ge.histmod(ihist)) .and.   &
                    (isccplowday(ilat,ilon,jour).lt.histmod(ihist+1)) )then          
                    hlow(ilon,ilat,ihist)=hlow(ilon,ilat,ihist)+1
                endif
             enddo

          monthisccpliq(ilon,ilat,1)=monthisccpliq(ilon,ilat,1)+             &
                                     isccpliqday(ilat,ilon,jour,1)

          monthisccpice(ilon,ilat,1)=monthisccpice(ilon,ilat,1)+             &
                                     isccpiceday(ilat,ilon,jour,1)      
do icat=1,catmax
          monthisccpun(ilon,ilat,1,icat)=monthisccpun(ilon,ilat,1,icat)+             &
                                     isccpunday(ilat,ilon,jour,1,icat)
enddo

          monthisccplow(ilon,ilat)=monthisccplow(ilon,ilat)+                 &
                                   isccplowday(ilat,ilon,jour)
          isccpdaypermonthlow(ilon,ilat)=isccpdaypermonthlow(ilon,ilat)+ 1 !isccpinddaylow(ilat,ilon,jour)
          endif

          if ( isccpinddaymid(ilat,ilon,jour).gt.0 ) then
             do ihist=1,histmax-1
                if( (isccpmidday(ilat,ilon,jour).ge.histmod(ihist)) .and.    &
                    (isccpmidday(ilat,ilon,jour).lt.histmod(ihist+1)) )then           
                    hmid(ilon,ilat,ihist)=hmid(ilon,ilat,ihist)+1
                endif
             enddo

          monthisccpliq(ilon,ilat,2)=monthisccpliq(ilon,ilat,2)+             &
                                     isccpliqday(ilat,ilon,jour,2)

          monthisccpice(ilon,ilat,2)=monthisccpice(ilon,ilat,2)+             &
                                     isccpiceday(ilat,ilon,jour,2)   

do icat=1,catmax
          monthisccpun(ilon,ilat,2,icat)=monthisccpun(ilon,ilat,2,icat)+             &
                                     isccpunday(ilat,ilon,jour,2,icat)
enddo
   
          monthisccpmid(ilon,ilat)=monthisccpmid(ilon,ilat)+                 &
                                   isccpmidday(ilat,ilon,jour)
          isccpdaypermonthmid(ilon,ilat)=isccpdaypermonthmid(ilon,ilat)+ 1 !isccpinddaymid(ilat,ilon,jour)
          endif

 
do nphase=1,4
   isccptemp=isccpiceday(ilat,ilon,jour,nphase)+isccpliqday(ilat,ilon,jour,nphase)
   if(inddayphase(ilat,ilon,jour,nphase).gt.0.)then
      monthisccpphase(ilon,ilat,nphase)=monthisccpphase(ilon,ilat,nphase)+   &
                                        isccpphaseday(ilat,ilon,jour,nphase)
      indmonthphase(ilon,ilat,nphase)=indmonthphase(ilon,ilat,nphase)+1
   else
   isccptemp=sum(isccpunday(ilat,ilon,jour,nphase,:))
      if(isccptemp.gt.0.)then
         indmonthphase2(ilon,ilat,nphase)=indmonthphase2(ilon,ilat,nphase)+1   
      endif
   endif
enddo


        enddo
     enddo
 enddo


!******** CALCULATION OF MONTHLY DIAGNOSTIC WITH MATCHING INDEXES ***********!

 do ilat=1,latmax-1
    do ilon=1,lonmax-1


        if ( indmonthcoltemp(ilon,ilat).ne.0 ) then
           monthcoltemp(ilon,ilat)=                                        &
           monthcoltemp(ilon,ilat)/indmonthcoltemp(ilon,ilat)
        else
           monthcoltemp(ilon,ilat)=-999
           do ihisttemp=1,histmax-1
              hcoltemp(ilon,ilat,ihisttemp)=-999
           enddo
        endif

        if ( indmonthhightemp(ilon,ilat).ne.0 ) then
           monthhightemp(ilon,ilat)=                                        &
           monthhightemp(ilon,ilat)/indmonthhightemp(ilon,ilat)
        else
           monthhightemp(ilon,ilat)=-999
           do ihisttemp=1,histmax-1
              hhightemp(ilon,ilat,ihisttemp)=-999
           enddo
        endif

        if ( indmonthmidtemp(ilon,ilat).ne.0 ) then
           monthmidtemp(ilon,ilat)=                                        &
           monthmidtemp(ilon,ilat)/indmonthmidtemp(ilon,ilat)
        else
           monthmidtemp(ilon,ilat)=-999
           do ihisttemp=1,histmax-1
              hmidtemp(ilon,ilat,ihisttemp)=-999
           enddo
        endif

        if ( indmonthlowtemp(ilon,ilat).ne.0 ) then
           monthlowtemp(ilon,ilat)=                                        &
           monthlowtemp(ilon,ilat)/indmonthlowtemp(ilon,ilat)
        else
           monthlowtemp(ilon,ilat)=-999
           do ihisttemp=1,histmax-1
              hlowtemp(ilon,ilat,ihisttemp)=-999
           enddo
        endif


        if ( isccpdaypermonth(ilon,ilat).ne.0 ) then
           monthisccphigh(ilon,ilat)=                                        &
           monthisccphigh(ilon,ilat)/isccpdaypermonth(ilon,ilat)

          monthisccpliq(ilon,ilat,3)=                                        &
          monthisccpliq(ilon,ilat,3)/isccpdaypermonth(ilon,ilat)
          monthisccpice(ilon,ilat,3)=                                        &
          monthisccpice(ilon,ilat,3)/isccpdaypermonth(ilon,ilat)

          monthisccpliq(ilon,ilat,4)=                                        &
          monthisccpliq(ilon,ilat,4)/isccpdaypermonth(ilon,ilat)
          monthisccpice(ilon,ilat,4)=                                        &
          monthisccpice(ilon,ilat,4)/isccpdaypermonth(ilon,ilat)


           monthcolcloud(ilon,ilat)=                                         &
           monthcolcloud(ilon,ilat)/isccpdaypermonth(ilon,ilat)
           monthcolclear(ilon,ilat)=                                         &
           1-monthcolcloud(ilon,ilat)

do icat=1,catmax
           monthisccpun(ilon,ilat,3,icat)=                                        &
           monthisccpun(ilon,ilat,3,icat)/isccpdaypermonth(ilon,ilat)

           monthisccpun(ilon,ilat,4,icat)=                                        &
           monthisccpun(ilon,ilat,4,icat)/isccpdaypermonth(ilon,ilat)
enddo

           if (indmonthheight(ilon,ilat).ne.0) then
              monthheight(ilon,ilat)=                                        &
              monthheight(ilon,ilat)/indmonthheight(ilon,ilat)
              monthheight2(ilon,ilat)=                                        &
              monthheight2(ilon,ilat)/indmonthheight(ilon,ilat)
           else
              monthheight(ilon,ilat)=-9999.
              monthheight2(ilon,ilat)=-9999.          
           endif
        else
        !   monthisccplow(ilon,ilat)=-9999
        !   monthisccpmid(ilon,ilat)=-9999

           monthisccpliq(ilon,ilat,4)=-9999.
           monthisccpice(ilon,ilat,4)=-9999.
           monthisccpun(ilon,ilat,4,:)=-9999.


           monthisccpliq(ilon,ilat,3)=-9999.
           monthisccpice(ilon,ilat,3)=-9999.
           monthisccpun(ilon,ilat,3,:)=-9999.

           monthisccphigh(ilon,ilat)=-9999.
           monthcolcloud(ilon,ilat)=-9999.
           monthcolclear(ilon,ilat)=-9999.
           monthheight(ilon,ilat)=-9999.
           monthheight2(ilon,ilat)=-9999.          
 
           indtotmean(ilon,ilat)=-999. 
          do ihist=1,histmax2-1
              hheight(ilon,ilat,ihist)=-999.
           enddo
           do ihist=1,histmax-1
              hhigh(ilon,ilat,ihist)=-999.
           enddo
        endif


!        if ( isccpdaypermonthlow(ilon,ilat).eq.0 ) then
        if ( isccpdaypermonthlow(ilon,ilat).gt.0 ) then
           monthisccplow(ilon,ilat)=                                         &
           monthisccplow(ilon,ilat)/isccpdaypermonthlow(ilon,ilat)

          monthisccpliq(ilon,ilat,1)=                                        &
          monthisccpliq(ilon,ilat,1)/isccpdaypermonthlow(ilon,ilat)
          monthisccpice(ilon,ilat,1)=                                        &
          monthisccpice(ilon,ilat,1)/isccpdaypermonthlow(ilon,ilat)

do icat=1,catmax
           monthisccpun(ilon,ilat,1,icat)=                                        &
           monthisccpun(ilon,ilat,1,icat)/isccpdaypermonthlow(ilon,ilat)
enddo

        else
           monthisccpice(ilon,ilat,1)=-9999.
           monthisccpliq(ilon,ilat,1)=-9999.
           monthisccpun(ilon,ilat,1,:)=-9999.
           monthisccplow(ilon,ilat)=-9999.

           do ihist=1,histmax-1
              hlow(ilon,ilat,ihist)=-999.
           enddo
        endif

!        if ( isccpdaypermonthmid(ilon,ilat).eq.0 ) then
        if ( isccpdaypermonthmid(ilon,ilat).gt.0 ) then
           monthisccpmid(ilon,ilat)=                                         &
           monthisccpmid(ilon,ilat)/isccpdaypermonthmid(ilon,ilat)

          monthisccpliq(ilon,ilat,2)=                                        &
          monthisccpliq(ilon,ilat,2)/isccpdaypermonthmid(ilon,ilat)
          monthisccpice(ilon,ilat,2)=                                        &
          monthisccpice(ilon,ilat,2)/isccpdaypermonthmid(ilon,ilat)

do icat=1,catmax
           monthisccpun(ilon,ilat,2,icat)=                                        &
           monthisccpun(ilon,ilat,2,icat)/isccpdaypermonthmid(ilon,ilat)
enddo  

        else
           monthisccpice(ilon,ilat,2)=-9999.
           monthisccpliq(ilon,ilat,2)=-9999.
           monthisccpun(ilon,ilat,2,:)=-9999.
           monthisccpmid(ilon,ilat)=-9999.

           do ihist=1,histmax-1
              hmid(ilon,ilat,ihist)=-999.
           enddo

        endif

        do nphase=1,4
           if ( indmonthphase(ilon,ilat,nphase).gt.0)then
              monthisccpphase(ilon,ilat,nphase)=                             &
              monthisccpphase(ilon,ilat,nphase)/indmonthphase(ilon,ilat,nphase)
           else
              if (indmonthphase2(ilon,ilat,nphase).gt.0)then
                 monthisccpphase(ilon,ilat,nphase)=-777.  
              else
                 monthisccpphase(ilon,ilat,nphase)=-9999.
              endif
           endif
        enddo


                                    
    enddo
enddo





!***************************** SAVE THE MAP FILES ***************************!

file8=trim(file6)//'.nc'    ! name of output ncdf map file
file9=trim(file3(25:55))    ! period of map file (description of ncdf file)


!print *, 'titi'
 call create_mapnc(file8,file9,lonmid,latmid,resd,dimidsm,gcm,lonmax-1,latmax-1)
 call map_recvar2nc2(monthisccplow,monthisccpmid,monthisccphigh,monthcolcloud,&
                     monthcolclear,dimidsm,file8,lonmax-1,latmax-1)


file8=trim(file66)//'.nc'    ! name of output ncdf map file
file9=trim(file3(25:55))    ! period of map file (description of ncdf file)

 call create_maphighnc(file8,file9,lonmid,latmid,resd,dimidsm,gcm,lonmax-1,latmax-1)
 call maphigh(monthisccphigh,monthheight,monthheight2,dimidsm,file8,lonmax-1,latmax-1)
!print *, 'titi2'

! Change NaN value from -9999 to -999 to fit with the GEWEX standard
forall(ilon=1:lonmax-1, ilat=1:latmax-1, monthisccplow(ilon,ilat)==-9999.)
monthisccplow(ilon,ilat)=-999
endforall
forall(ilon=1:lonmax-1, ilat=1:latmax-1, monthisccpmid(ilon,ilat)==-9999.)
monthisccpmid(ilon,ilat)=-999
endforall
forall(ilon=1:lonmax-1, ilat=1:latmax-1, monthisccphigh(ilon,ilat)==-9999.)
monthisccphigh(ilon,ilat)=-999
endforall
forall(ilon=1:lonmax-1, ilat=1:latmax-1, monthcolcloud(ilon,ilat)==-9999.)
monthcolcloud(ilon,ilat)=-999
endforall





file8=trim(file11)//'.nc'    ! name of output ncdf map file

!print *, 'titi3'

 call create_mapnc_phase(file8,file9,lonmid,latmid,resd,dimidsm,dimidsm2,gcm,lonmax-1,latmax-1)
! print *, 'titi4'

call map_recvar2nc2phaseocc2(monthisccpliq,monthisccpice,monthisccpun,        &
                              monthisccpphase,dimidsm,dimidsm2,file8,          &
                              lonmax-1,latmax-1)

!print *, 'titi5'

!subroutine map_recvar2nc2phaseocc2(liq,ice,ho,un,dust,dim,fname,nlon,nlat)

                 !   monthlowtemp,monthmidtemp,monthhightemp,    &
                 !   monthcoltemp,hlowtemp,hmidtemp,hhightemp,hcoltemp,dimidhist3)



! Deallocate daily & monthly map variables
print *, 'deallocate daily & monthly map variables'

if(model=='lmdz')then
  deallocate(hlow,hmid,hhigh,hcol,hheight)
  deallocate(monthisccplow,monthisccpmid,monthisccphigh)
  deallocate(monthcolcloud,isccpdaypermonth)!,monthcolclear
  deallocate(monthisccpliq,monthisccpice,monthisccpun)
  deallocate(indmonthphase,monthisccpphase)
  deallocate(indmonthphase2)
  deallocate(inddayphase)
  deallocate(isccplowday,isccpmidday,isccphighday)
  deallocate(isccpliqday,isccpiceday,isccpunday)
  deallocate(colcloudday,isccpindday)!,colclearday
  deallocate(isccpinddaylow,isccpinddaymid)
  deallocate(isccpdaypermonthlow,isccpdaypermonthmid)
  deallocate(indtotmean,indtot)
  deallocate(indmonthheight,monthheight,indheight,heightday)
  deallocate(monthheight2,heightday2)
  deallocate(lowtemp,midtemp,hightemp,coltemp)
  deallocate(indlowtemp,indmidtemp,indhightemp,indcoltemp)
  deallocate(hlowtemp,hmidtemp,hhightemp,hcoltemp)

!  deallocate(indtotmean,indretmean)
!  deallocate(indretlowmean,indretmidmean,indrethighmean)
endif

print *, 'map file recorded'



!goto 621
!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!! PART III : CLOUDY MAP3D FILES !!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


! Allocation / initialization of MAP3D monthly variables
print *, 'allocation / initialization of MAP3D monthly variables'
   allocate(indphasepermonth(lonmax-1,latmax-1,altmax))!,                               &

   allocate(indpermonth(lonmax-1,latmax-1,altmax))!,                               &
         !   indpermonthtot(lonmax-1,latmax-1,altmax))
   allocate(monthcloudfract(lonmax-1,latmax-1,altmax),                           &
            monthclearfract(lonmax-1,latmax-1,altmax)) 
  ! allocate(monthsatfract(lonmax-1,latmax-1,altmax))
   allocate(monthuncertfract(lonmax-1,latmax-1,altmax))!,                          &
         !   monthnanfract(lonmax-1,latmax-1,altmax),                             &
        !    monthsefract(lonmax-1,latmax-1,altmax))
   allocate(monthicecloud(lonmax-1,latmax-1,altmax),                       &
            monthwatercloud(lonmax-1,latmax-1,altmax),                       &
            indphasemonth(lonmax-1,latmax-1,altmax))
   allocate(monthuncloud(lonmax-1,latmax-1,altmax,catmax))
   allocate(monthphasecloud(lonmax-1,latmax-1,altmax))
   allocate(indmonthphase3D(lonmax-1,latmax-1,altmax))
   allocate(indphasefractday(latmax-1,lonmax-1,altmax,daymax))
   allocate(monthcftemp(lonmax-1,latmax-1,tempmax-1))
   allocate(monthcftempice(lonmax-1,latmax-1,tempmax-1))
   allocate(monthcftempliq(lonmax-1,latmax-1,tempmax-1))
   allocate(monthcftempphase(lonmax-1,latmax-1,tempmax-1))
   allocate(indmonthphasetemp(lonmax-1,latmax-1,tempmax-1))
   allocate(indcftemppermonth(lonmax-1,latmax-1,tempmax-1))
allocate(cftempphaseday(latmax-1,lonmax-1,tempmax-1,daymax))

!allocate(indtest(lonmax-1,latmax-1,tempmax-1))
!indtest(:,:,:)=0;

cftempphaseday(:,:,:,:)=0;
indpermonth(:,:,:)=0;indphasepermonth(:,:,:)=0;
indmonthphase3D(:,:,:)=0; monthphasecloud(:,:,:)=0;
monthcloudfract(:,:,:)=0;monthclearfract(:,:,:)=0;!monthsatfract(:,:,:)=0;
monthuncertfract(:,:,:)=0;!monthnanfract(:,:,:)=0;monthsefract(:,:,:)=0;
indphasemonth(:,:,:)=0;monthicecloud(:,:,:)=0;monthwatercloud(:,:,:)=0;
monthuncloud(:,:,:,:)=0;
indphasefractday(:,:,:,:)=0;
monthcftemp(:,:,:)=0
monthcftempice(:,:,:)=0
monthcftempliq(:,:,:)=0
monthcftempphase(:,:,:)=0
indmonthphasetemp(:,:,:)=0
indcftemppermonth(:,:,:)=0



!! OCCURRENCES FILE IN 3D WITH TEMPERATURE
  do jour=1,31
    do itemp=1,tempmax-1
       do ilon=1,lonmax-1
         do ilat=1,latmax-1
           if (indcftemp(ilat,ilon,itemp,jour).gt.0) then
     indcftemppermonth(ilon,ilat,itemp) = indcftemppermonth(ilon,ilat,itemp)+indcftemp(ilat,ilon,itemp,jour)
     monthcftemp(ilon,ilat,itemp)=monthcftemp(ilon,ilat,itemp)   &
                                +cftempday(ilat,ilon,itemp,jour)
     monthcftempice(ilon,ilat,itemp)=monthcftempice(ilon,ilat,itemp)   &
                                +cftempiceday(ilat,ilon,itemp,jour)
     monthcftempliq(ilon,ilat,itemp)=monthcftempliq(ilon,ilat,itemp)   &
                                +cftempliqday(ilat,ilon,itemp,jour)
          endif
         enddo
       enddo
    enddo
  enddo


file8=trim(file13)//'_occ.nc'   ! name of output netcdf MAP3D file
file9=trim(file3(25:55))   ! period of MAP3D file (description of ncdf file) 

print *, 'MAP3D PHASE files recorded'


 call create_temp3d(file8,file9,lonmid,latmid,tempmid,tempmod_bound,resd,       &
                    dimidsp,tempmax-1,lonmax-1,latmax-1)


print *, 'MAP3D create temp'

 call temp_recvar2nc(monthcftemp,monthcftempliq,monthcftempice,        &
                     indcftemppermonth,dimidsp,file8,tempmax-1,lonmax-1,latmax-1)

monthcftemp(:,:,:)=0
monthcftempice(:,:,:)=0
monthcftempliq(:,:,:)=0
indcftemppermonth(:,:,:)=0


!print *, 'tata1'

!********* CALCULATION OF DAILY DIAGNOSTIC WITH MATCHING INDEXES ************!
!  do ilat=1,latmax-1
!!    do ilon=1,lonmax-1
!       do jour=1,31
!         do ialt=1,altmax

cfsumtemp=0.

  do jour=1,31
    do itemp=1,tempmax-1
       do ilon=1,lonmax-1
         do ilat=1,latmax-1
           if (indcftemp(ilat,ilon,itemp,jour).gt.0) then
!indtest(ilon,ilat,itemp)=indtest(ilon,ilat,itemp)+indcftemp(ilat,ilon,itemp,jour)

cfsumtemp=cftempliqday(ilat,ilon,itemp,jour)+  &
          cftempiceday(ilat,ilon,itemp,jour)

if(cfsumtemp.gt.0)then
   if(cftempiceday(ilat,ilon,itemp,jour).eq.0.)then
      cftempphaseday(ilat,ilon,itemp,jour)=0.
      indcftempphase(ilat,ilon,itemp,jour)=  &
      indcftempphase(ilat,ilon,itemp,jour)+1
   else
      cftempphaseday(ilat,ilon,itemp,jour)=    &
      cftempiceday(ilat,ilon,itemp,jour)/cfsumtemp
      indcftempphase(ilat,ilon,itemp,jour)=  &
      indcftempphase(ilat,ilon,itemp,jour)+1
   endif
endif

cftempday(ilat,ilon,itemp,jour)=   &
cftempday(ilat,ilon,itemp,jour)/indcftemp(ilat,ilon,itemp,jour)
cftempliqday(ilat,ilon,itemp,jour)=   &
cftempliqday(ilat,ilon,itemp,jour)/indcftemp(ilat,ilon,itemp,jour)
cftempiceday(ilat,ilon,itemp,jour)=   &
cftempiceday(ilat,ilon,itemp,jour)/indcftemp(ilat,ilon,itemp,jour)

if(cftempday(ilat,ilon,itemp,jour).gt.1.)then
print *, ilat,ilon,itemp,jour
print *, cftempday(ilat,ilon,itemp,jour)
print *, cftempliqday(ilat,ilon,itemp,jour)
print *, cftempiceday(ilat,ilon,itemp,jour)

endif

           endif
         enddo
       enddo
    enddo
  enddo


allocate(tot_ind(lonmax-1,latmax-1,altmax), &
         cloud_ind(lonmax-1,latmax-1,altmax), &
         ice_ind(lonmax-1,latmax-1,altmax), &
         water_ind(lonmax-1,latmax-1,altmax), &
         un_ind(lonmax-1,latmax-1,altmax,catmax))


print *, 'test'

    do ialt=1,altmax
       do ilon=1,lonmax-1
         do ilat=1,latmax-1

    cloud_ind(ilon,ilat,ialt)=sum(cloudfractday(ilat,ilon,ialt,:))
    tot_ind(ilon,ilat,ialt)=sum(indday(ilat,ilon,ialt,:))
    ice_ind(ilon,ilat,ialt)=sum(icecloudfractday(ilat,ilon,ialt,:))
    water_ind(ilon,ilat,ialt)=sum(watercloudfractday(ilat,ilon,ialt,:))

         do icat=1,catmax
    un_ind(ilon,ilat,ialt,icat)=sum(uncloudfractday(ilat,ilon,ialt,:,icat))
         enddo

       enddo
    enddo
  enddo



file8=trim(file10)//'_occ.nc'   ! name of output netcdf MAP3D file
file9=trim(file3(25:55))   ! period of MAP3D file (description of ncdf file) 

print *, 'titi'
 call create_ind3d(file8,file9,lonmid,latmid,altmid,altmod_bound,resd,   &
                       dimidsp,dimidsp2,altmax,lonmax-1,latmax-1)

print *, 'titi2'

 call record_ind3d(cloud_ind,tot_ind,ice_ind,water_ind,un_ind, &
                  dimidsp,dimidsp2,file8,altmax,lonmax-1,latmax-1)
print *, 'titi3'

deallocate(tot_ind,cloud_ind,ice_ind,water_ind,un_ind)


isccptemp=0.
  do jour=1,31
    do ialt=1,altmax
       do ilon=1,lonmax-1
         do ilat=1,latmax-1
           if (indday(ilat,ilon,ialt,jour).gt.0) then
! Cloud fraction / Clear fraction / Saturated fraction, daily
         cloudfractday(ilat,ilon,ialt,jour) =                                &
         cloudfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)        
         clearfractday(ilat,ilon,ialt,jour) =                                &
         clearfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)
        ! satfractday(ilat,ilon,ialt,jour) =                                  &
        ! satfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)
         uncertfractday(ilat,ilon,ialt,jour) =                               &
         uncertfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)
         !nanfractday(ilat,ilon,ialt,jour) =                                  &
        ! nanfractday(ilat,ilon,ialt,jour)/inddaytot(ilat,ilon,ialt,jour)
        ! sefractday(ilat,ilon,ialt,jour) =                                   &
        ! sefractday(ilat,ilon,ialt,jour)/inddaytot(ilat,ilon,ialt,jour)
         !  endif

        !   if (indphaseday(ilat,ilon,ialt,jour).gt.0) then
         isccptemp=icecloudfractday(ilat,ilon,ialt,jour)+                     &
                   watercloudfractday(ilat,ilon,ialt,jour)
if(isccptemp.gt.0)then
   if(icecloudfractday(ilat,ilon,ialt,jour).eq.0)then
   phasefractday(ilat,ilon,ialt,jour) = 0.
   indphasefractday(ilat,ilon,ialt,jour) = indphasefractday(ilat,ilon,ialt,jour)+1
   else
   indphasefractday(ilat,ilon,ialt,jour) = indphasefractday(ilat,ilon,ialt,jour)+1
   phasefractday(ilat,ilon,ialt,jour) =                                        &
   icecloudfractday(ilat,ilon,ialt,jour)/isccptemp
   endif
else
         isccptemp=sum(uncloudfractday(ilat,ilon,ialt,jour,:))
         if(isccptemp.gt.0)then
           phasefractday(ilat,ilon,ialt,jour)=-777. 
         endif
endif

         icecloudfractday(ilat,ilon,ialt,jour) =                                &
         icecloudfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)
         watercloudfractday(ilat,ilon,ialt,jour) =                                &
         watercloudfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)

do icat=1,catmax
         uncloudfractday(ilat,ilon,ialt,jour,icat) =                                &
         uncloudfractday(ilat,ilon,ialt,jour,icat)/indday(ilat,ilon,ialt,jour)
enddo


          endif

         enddo
       enddo
    enddo
  enddo



!************ INCREMENTING THE MONTHLY DIAGNOSTIC AND INDEXES **************!


  do jour=1,31
    do itemp=1,tempmax-1
       do ilon=1,lonmax-1
         do ilat=1,latmax-1
           if (indcftemp(ilat,ilon,itemp,jour).gt.0) then
     indcftemppermonth(ilon,ilat,itemp) = indcftemppermonth(ilon,ilat,itemp)+1
     monthcftemp(ilon,ilat,itemp)=monthcftemp(ilon,ilat,itemp)   &
                                +cftempday(ilat,ilon,itemp,jour)
     monthcftempice(ilon,ilat,itemp)=monthcftempice(ilon,ilat,itemp)   &
                                +cftempiceday(ilat,ilon,itemp,jour)
     monthcftempliq(ilon,ilat,itemp)=monthcftempliq(ilon,ilat,itemp)   &
                                +cftempliqday(ilat,ilon,itemp,jour)

!indcftempphase
if(indcftempphase(ilat,ilon,itemp,jour).gt.0.)then
    monthcftempphase(ilon,ilat,itemp) = monthcftempphase(ilon,ilat,itemp)+  &
                                      cftempphaseday(ilat,ilon,itemp,jour)
    indmonthphasetemp(ilon,ilat,itemp) = indmonthphasetemp(ilon,ilat,itemp)+1 
endif

           endif
         enddo
       enddo
    enddo
  enddo

  do jour=1,31
    do ialt=1,altmax
      do ilat=1,latmax-1
       do ilon=1,lonmax-1




  if (indday(ilat,ilon,ialt,jour).gt.0) then
     indpermonth(ilon,ilat,ialt) = indpermonth(ilon,ilat,ialt)+1
    indphasepermonth(ilon,ilat,ialt) = indphasepermonth(ilon,ilat,ialt)+ &
                                      indday(ilat,ilon,ialt,jour)

     ! Monthly Fraction 
     monthcloudfract(ilon,ilat,ialt) = monthcloudfract(ilon,ilat,ialt) +     &
                                       cloudfractday(ilat,ilon,ialt,jour)
     monthclearfract(ilon,ilat,ialt) = monthclearfract(ilon,ilat,ialt) +     &
                                       clearfractday(ilat,ilon,ialt,jour)
  !   monthsatfract(ilon,ilat,ialt) = monthsatfract(ilon,ilat,ialt) +         &
  !                                   satfractday(ilat,ilon,ialt,jour)
     monthuncertfract(ilon,ilat,ialt) = monthuncertfract(ilon,ilat,ialt) +   &
                                        uncertfractday(ilat,ilon,ialt,jour)

    monthicecloud(ilon,ilat,ialt) = monthicecloud(ilon,ilat,ialt) +     &
                                    icecloudfractday(ilat,ilon,ialt,jour)
    monthwatercloud(ilon,ilat,ialt) = monthwatercloud(ilon,ilat,ialt) +     &
                                    watercloudfractday(ilat,ilon,ialt,jour)

do icat=1,catmax
    monthuncloud(ilon,ilat,ialt,icat) = monthuncloud(ilon,ilat,ialt,icat) +     &
                                    uncloudfractday(ilat,ilon,ialt,jour,icat)
enddo

!!! changement ici 
!! a terminer la phase 3D
if(indphasefractday(ilat,ilon,ialt,jour).gt.0.)then
    monthphasecloud(ilon,ilat,ialt) = monthphasecloud(ilon,ilat,ialt) +    &
                                      phasefractday(ilat,ilon,ialt,jour)
    indmonthphase3D(ilon,ilat,ialt) = indmonthphase3D(ilon,ilat,ialt)+1 
endif

  endif



       enddo
      enddo
    enddo
  enddo

!print *, 'tata3'

    do ialt=1,altmax
      do ilat=1,latmax-1
       do ilon=1,lonmax-1

  if (indphasepermonth(ilon,ilat,ialt).lt.sum(monthuncloud(ilon,ilat,ialt,:))) then

print *, ilon,ilat,ialt,indphasepermonth(ilon,ilat,ialt),monthuncloud(ilon,ilat,ialt,:)
  endif

enddo
enddo
enddo


! Deallocate daily MAP3D variables
print *, 'deallocate daily MAP3D variables'

if(model=='lmdz')then
 deallocate(indday)!,inddaytot)
 deallocate(cloudfractday, clearfractday,uncertfractday)
 deallocate(indphaseday,icecloudfractday,watercloudfractday)
 deallocate(uncloudfractday)
 deallocate(cftempday,cftempliqday,cftempiceday,indcftemp,indcftempphase)
! deallocate(satfractday, nanfractday,sefractday)
endif


!******** CALCULATION OF MONTHLY DIAGNOSTIC WITH MATCHING INDEXES ***********!

do ilat=1,latmax-1
    do ilon=1,lonmax-1

       do itemp=1,tempmax-1
          if(indmonthphasetemp(ilon,ilat,itemp).gt.0)then
          monthcftempphase(ilon,ilat,itemp)=                                 &
          monthcftempphase(ilon,ilat,itemp)/indmonthphasetemp(ilon,ilat,itemp)
          else
          monthcftempphase(ilon,ilat,itemp)=-9999.   
          endif

          if(indcftemppermonth(ilon,ilat,itemp).gt.0 )then
          monthcftemp(ilon,ilat,itemp)=    &
          monthcftemp(ilon,ilat,itemp)/indcftemppermonth(ilon,ilat,itemp)
          monthcftempice(ilon,ilat,itemp)=    &
          monthcftempice(ilon,ilat,itemp)/indcftemppermonth(ilon,ilat,itemp)
          monthcftempliq(ilon,ilat,itemp)=    &
          monthcftempliq(ilon,ilat,itemp)/indcftemppermonth(ilon,ilat,itemp)

          else
          monthcftemp(ilon,ilat,itemp)=-9999.   
          monthcftempliq(ilon,ilat,itemp)=-9999.   
          monthcftempice(ilon,ilat,itemp)=-9999.   

          endif
       enddo


       do ialt=1,altmax    

  

        if ( indmonthphase3D(ilon,ilat,ialt).gt.0 ) then
        monthphasecloud(ilon,ilat,ialt)=                                    &
        monthphasecloud(ilon,ilat,ialt)/indmonthphase3D(ilon,ilat,ialt)  
        else
        monthphasecloud(ilon,ilat,ialt)=-9999.     
           do jour = 1,31
              if(phasefractday(ilat,ilon,ialt,jour).eq.-777.)then
              monthphasecloud(ilon,ilat,ialt)=-777.
              endif
           enddo
           
        endif

        if ( indpermonth(ilon,ilat,ialt).gt.0 ) then
        monthcloudfract(ilon,ilat,ialt)=                                     &
        monthcloudfract(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)
        monthclearfract(ilon,ilat,ialt)=                                     &
        monthclearfract(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)
     !   monthsatfract(ilon,ilat,ialt)=                                       &
     !   monthsatfract(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)
        monthuncertfract(ilon,ilat,ialt)=                                    &
        monthuncertfract(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)

        monthicecloud(ilon,ilat,ialt)=                                    &
        monthicecloud(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)
        monthwatercloud(ilon,ilat,ialt)=                                    &
        monthwatercloud(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)
do icat=1,catmax
        monthuncloud(ilon,ilat,ialt,icat)=                                    &
        monthuncloud(ilon,ilat,ialt,icat)/indpermonth(ilon,ilat,ialt)
enddo

        else
        monthcloudfract(ilon,ilat,ialt)=-9999.
        monthclearfract(ilon,ilat,ialt)=-9999.
      !  monthsatfract(ilon,ilat,ialt)=-9999
        monthuncertfract(ilon,ilat,ialt)=-9999.
        monthuncloud(ilon,ilat,ialt,:)=-9999.
        monthwatercloud(ilon,ilat,ialt)=-9999.
        monthicecloud(ilon,ilat,ialt)=-9999. 
        indphasepermonth(ilon,ilat,ialt)=-9999. 
       endif

      enddo
    enddo
  enddo


deallocate(phasefractday,indphasefractday)
! Deallocate monthly indexes

!**************************** SAVE THE MAP3D FILES **************************!

file8=trim(file5)//trim(version)//'.nc'   ! name of output netcdf MAP3D file
file9=trim(file3(25:55))   ! period of MAP3D file (description of ncdf file) 



 call create_profnc(file8,file9,lonmid,latmid,altmid,altmod_bound,resd,       &
                    dimidsp,altmax,lonmax-1,latmax-1)


 call prof_recvar2nc(monthcloudfract,monthclearfract,monthuncertfract,        &
                     dimidsp,file8,altmax,lonmax-1,latmax-1)

print *, 'MAP3D files recorded'

file8=trim(file10)//'.nc'   ! name of output netcdf MAP3D file
file9=trim(file3(25:55))   ! period of MAP3D file (description of ncdf file) 

 call create_depolnc3d(file8,file9,lonmid,latmid,altmid,altmod_bound,resd,   &
                       dimidsp,dimidsp2,altmax,lonmax-1,latmax-1)

!!$ call depol_recvar2nc(monthicecloud,monthwatercloud, & !indphasemonth, &!,monthsatfract,          &,monthnanfract,monthsefract
!!$                     dimidsp,    &
!!$                     file8,altmax,lonmax-1,latmax-1)
!subroutine depol_recvar2ncocc(ice,water,un,ho,dust,ind,dim,fname,alt,nlon,nlat)!nan,se,sat

 call depol_recvar2ncocc(monthicecloud,monthwatercloud,monthuncloud,          &
                         monthphasecloud,indphasepermonth,dimidsp,dimidsp2,   &
                         file8,altmax,lonmax-1,latmax-1)


! check the allocation 
 if (OK_buffer/=0) print *,'--- buffer allocation error '   
!print *, 'tutu'
  deallocate(monthphasecloud,stat = OK_buffer)
if (OK_buffer/=0) print *,'--- buffer allocation error ' 
!print *, 'tutu2'
!  deallocate(indmonthphase3D,stat = OK_buffer)
!if (OK_buffer/=0) print *,'--- buffer allocation error ' 
!print *, 'tutu3'
  deallocate(indphasepermonth,stat = OK_buffer)
if (OK_buffer/=0) print *,'--- buffer allocation error ' 


file8=trim(file13)//'.nc'   ! name of output netcdf MAP3D file
file9=trim(file3(25:55))   ! period of MAP3D file (description of ncdf file) 

print *, 'MAP3D PHASE files recorded'


 call create_temp3d(file8,file9,lonmid,latmid,tempmid,tempmod_bound,resd,       &
                    dimidsp,tempmax-1,lonmax-1,latmax-1)


print *, 'MAP3D create temp'

 call temp_recvar2nc(monthcftemp,monthcftempliq,monthcftempice,        &
                     monthcftempphase,dimidsp,file8,tempmax-1,lonmax-1,latmax-1)

 
print *, 'MAP3D record temp'


! Deallocate daily & monthly diagSR variables
print *, 'deallocate daily & monthly MAP3D variables'

if(model=='lmdz')then
  deallocate(monthcloudfract,monthclearfract,monthuncertfract)
  deallocate(monthicecloud,monthwatercloud,indphasemonth,indpermonth)
  deallocate(monthuncloud)
!print *, 'titi'
 ! deallocate(monthsatfract,monthnanfract,monthsefract)
!print *, 'titi'
!print *, 'titi'
  deallocate(indcftemppermonth)
deallocate(indmonthphasetemp)
deallocate(monthcftempphase)
!print *, 'titi'
  deallocate(monthcftemp)
deallocate(monthcftempliq)
deallocate(monthcftempice)
!print *, 'titi'
endif






!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!! PART II : DIAGSR FILES !!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


! Allocation / initialization of diagSR monthly variables
print *, 'allocation / initialization of diagSR monthly variables'




file8=trim(file7)//'.nc'   ! name of output netcdf diagSR file
file9=trim(file3(25:55))   ! period of diagSR file (description of ncdf file)

!!$ call create_diagnc2(trim(file7)//'.tmp',file9,lonmid,latmid,altmid,srmod,resd,dimidsd2, &!dimidsd2,     &
!!$                    altmax,lonmax-1,latmax-1)
!!$
!!$
!!$
!!$!************ CALCULATION OF MONTHLY DIAGNOSTIC WITH DAILY VAR **************!


print *, 'allocation monthdiagSR terminé'


print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,1),4),3),2),1)
print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,2),4),3),2),1)
print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,3),4),3),2),1)

forall(ilon=1:lonmax-1, ilat=1:latmax-1,  ialt=1:altmax, indnan(ilat,ilon,ialt) == 0)
 diagSR(ilon,ilat,ialt,1:diagmax-1)=-9999.
 diagSRpha(ilon,ilat,ialt,1:diagmax-8,:)=-9999.

endforall



 call create_diagnc(file8,file9,lonmid,latmid,altmid,altmod_bound,srmod,resd,dimidsd,dimidsdb, &!dimidsd2,     &
                    altmax,lonmax-1,latmax-1)

print *, 'creation fichier diag final'

call diag_recvar2nc3(diagSR,dimidsd,dimidsdb,file8,altmax,lonmax-1,latmax-1)
! call diag_recvar2nc(monthdiagSR15,monthdiagSR1,dimidsd,dimidsd2,file8,altmax,lonmax-1,latmax-1)

file8=trim(file12)//'.nc'   ! name of output netcdf diagSR file
file9=trim(file3(25:55))   ! period of diagSR file (description of ncdf file)

print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,1),4),3),2),1)
print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,2),4),3),2),1)
print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,3),4),3),2),1)

 call create_diagncpha(file8,file9,lonmid,latmid,altmid,altmod_bound,srmod,resd,dimidsd,dimidsdb, &!dimidsd2,     &
                    altmax,lonmax-1,latmax-1)

print *, 'creation fichier diag final'

call diag_recvar2nc3pha(diagSRpha,dimidsd,dimidsdb,file8,altmax,lonmax-1,latmax-1)

print *, 'deallocate daily & monthly diagSR variables'

if(model=='lmdz')then
   deallocate(diagSR,diagSRpha)!,monthdiagSR15,monthdiagSR1)
endif
! Deallocate daily & monthly diagSR variables
print *, 'deallocate daily & monthly diagSR variables'


621 continue



if(model=='lmdz')then
   deallocate(latmod,lonmod,prestop,altmod,srmod,pr2mod,atbrmod,srdepmod,depolmod,lonmid,latmid,altmid,tempmod,stat = OK_buffer) !crmod
endif

close(1)


CASE ("chimere")
continue
 
CASE ("wrf")
   deallocate(latmod,lonmod,prestop,altmod,srmod,lonmid,latmid,altmid,stat = OK_buffer) !crmod
   close(1)
 !  close(35)

CASE DEFAULT
print *, "error" 

ENDSELECT

 102  format(f10.2,3e13.5)
 103  format((2x,f10.2),5(2x,e13.5))
 104 format(2(2x,I5),(2x,e13.5),4(2x,f8.2))
 106  format(3(f8.3,1x),2(e12.5,1x),7(f8.3,1x))
 107  format(9(f8.2,1x))
 108 format(4(f8.3,1x),(1x,e13.6))

print *, 'END OF PROGRAM'
print *, '-----------------------------------------------------------------------'
print *, '-----------------------------------------------------------------------'






contains


!****************************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************!

!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                            END OF PROGRAM                                *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

!****************************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************!





!****************************************************************************!
!****************************************************************************!
!**************************    SUBROUTINE LIST    ***************************!
!****************************************************************************!
!****************************************************************************!



!----------------------------------------------------------------------------!
! **** SDSREAD8 ****  This routine allows to read 8-bit SDS-variable on hdf  !
!                     files                                                  !
!----------------------------------------------------------------------------!
! filename  : calipso file name
! varname   : SDS 8-bit name variable to extract from calipso file
!----------------------------------------------------------------------------!
! var       : SDS 8-bit variable to extract from calipso file
!----------------------------------------------------------------------------!
! i         : loop on the indice of calipso variables
! name      : SDS 8-bit name variable to extract from calipso file
! ret       : 
! start     :
! stride    :
! edges     :
! sd_id     : some id
! sds_id    : some id
! OK_buffer :
! npts      : value of one dimension of the read var
! nprofs    : value of one dimension of the read var
! nb        : some loop index
! index     :
! rank      :
! istat     :
! attributes:
! num_type  :
!                                                                            !
!----------------------------------------------------------------------------!
subroutine sdsread8(var,filename,varname,ret)

	implicit none
	include "hdf.f90"
	include "dffunc.f90"
 
	integer 	:: i, ret, start(2) = 0, stride(2) = 1 ,edges(2)
	integer 	:: sd_id, sds_id, OK_buffer
        integer         :: npts, nprofs, nb
        real*8,dimension(:,:),allocatable :: var
	integer*4	:: index, rank, istat, attributes, num_type
	integer*4 	:: dim_sizes (32)
	character(len=232) :: name, filename
        character       :: varname*100
    
	sd_id = sfstart (filename,	DFACC_READ)

! loop on the variables number (<100)
do i=0,100

   nb=i
	! Select and read the var
	sds_id = sfselect (sd_id, nb)
        istat = sfginfo (sds_id, name, rank, dim_sizes, num_type, attributes)
                       
        if(name==trim(varname))exit
          
 enddo

       ! Retrieving the dimesions & allocating memorie
        npts = dim_sizes(1)
        nprofs = dim_sizes(2)

      allocate(var(npts, nprofs),stat = OK_buffer)
      if (OK_buffer/=0) print *,'--- buffer allocation of ',trim(name),'error' 

        edges = [npts, nprofs]
        ret = sfrdata (sds_id, start, stride, edges, var)

	! Do something with the data
      if (ret.eq.-1) then
		print *,'ERROR'
goto 555
      end if

	! Close the crap
	ret = sfendacc(sds_id)
	ret = sfend(sd_id)
555 continue


end subroutine sdsread8
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! **** SDSREAD ****  This routine allows to read a 2dim SDS-variable on hdf  !
!                    files                                                   !
!----------------------------------------------------------------------------!
subroutine sdsread(var,filename,varname)

	implicit none
	include "hdf.f90"
	include "dffunc.f90"
 
	integer 	:: i, ret, start(2) = 0, stride(2) = 1 ,edges(2)
	integer 	:: sd_id, sds_id
        integer         :: npts, nprofs, nb, OK_buffer
        real,dimension(:,:),allocatable :: var
	integer*4	:: index, rank, istat, attributes, num_type
	integer*4 	:: dim_sizes (32)
	character(len=232) :: name, filename
         character       :: varname*100
 
!	print *,'Reading HDF file...'
     
	sd_id = sfstart (filename,	DFACC_READ)
  

do i=2,100

   nb=i
	! Select and read the var
	sds_id = sfselect (sd_id, nb)
        istat = sfginfo (sds_id, name, rank, dim_sizes, num_type, attributes)

        if(name==trim(varname))exit

enddo

      !  print *, name,'dimensions : ', dim_sizes(1:2)
        npts = dim_sizes(1)
        nprofs = dim_sizes(2)
      allocate(var(npts, nprofs),stat = OK_buffer)
      if (OK_buffer/=0) print *,'--- buffer allocation of ',trim(name),'error' 
        edges = [npts, nprofs]
        ret = sfrdata (sds_id, start, stride, edges, var)
!	print *,'Reading : ', name

	! Do something with the data
      if (ret.eq.-1) then
		print *,'ERROR'
      end if

	! Close the crap
	ret = sfendacc(sds_id)
	ret = sfend(sd_id)

end subroutine sdsread
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
!*****  METAREAD ***** This routine allows to read a Meta-variable on hdf    !
!                      files                                                 !
!----------------------------------------------------------------------------!
subroutine metaread(var,varname,filename)
 
	implicit none
	include "hdf.f90"
	include "dffunc.f90"
 
	integer 	:: ret, istat, OK_buffer
	integer		:: file_id, vdata_ref, vdata_id
!	integer		:: interlace, vdata_size ,n_rec
        character       :: varname*30
	character	:: filename*1024
!	character  	:: fields*550, vdata_name*30
	
 	real*4,dimension(:),allocatable  ::  var	
  
! Allocate var size
   if(varname=='Lidar_Data_Altitudes')then
      allocate(var(583), stat=OK_buffer)
   if (OK_buffer/=0) print *,'--- buffer allocation of ',trim(varname),'error' 
        else
      allocate(var(33), stat=OK_buffer)
   if (OK_buffer/=0) print *,'--- buffer allocation of ',trim(varname),'error' 
        endif

	!print *,'Reading HDF file...'

	file_id = hopen (filename,	DFACC_READ, 0)
 
	! initialize vdata interface,find vdata called 'metadata '
	! attach to this vdata
	istat = vfstart (file_id)
	vdata_ref = vsffnd (file_id, 'metadata')
	vdata_id = vsfatch (file_id, vdata_ref, 'r')

	! reads Varname
	!print *, varname,' exists : ', vsfex (vdata_id, varname)
	istat = vsfsfld (vdata_id, varname)
	ret = vsfread (vdata_id, var, 1, FULL_INTERLACE)
	if (ret.ne.-1) then
           continue
        else
      !     print *, 'ERROR READING METAVAR'
	end if
 
 	! TRES IMPORTANT : il faut appeler vsfseek entre chaque vsfsfld
	! sinon la lecture ne fonctionne pas !
 
	ret = vsfseek (vdata_id, 0)
 
	! Close all the crap
	istat = vsfdtch(vdata_id)
	istat = vfend (file_id)
	istat = hclose (file_id)
 
	if (istat.eq.0) then
		print *,'Reading HDF File OK.'
	else
		print *,'ERROR READING HDF'
        endif

end subroutine metaread
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! ***** NPROF ***** This routine find the hdf file number of profiles        !
!----------------------------------------------------------------------------!
subroutine nprof(filename,nb,nprofs)

	implicit none
	include "hdf.f90"
	include "dffunc.f90"
 
	integer 	:: ret
	integer 	:: sd_id, sds_id
        integer         :: npts, nprofs, nb
	integer*4	:: rank, istat, attributes, num_type !index
	integer*4 	:: dim_sizes (32)
	character(len=232) :: name, filename
 
     
	sd_id = sfstart (filename,	DFACC_READ)
  
	! Select and read the var
	sds_id = sfselect (sd_id, nb)
        istat = sfginfo (sds_id, name, rank, dim_sizes, num_type, attributes)
        
        ! Rtrieving the dimesions
        npts = dim_sizes(1)
        nprofs = dim_sizes(2)
print *, 'le nombre de profil du fichier est :', nprofs

	! Close the crap
	ret = sfendacc(sds_id)
	ret = sfend(sd_id)

end subroutine nprof
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! ***** INTERP ***** This routine proceed a linear interpolation of variables!
!                    entered                                                 !
!----------------------------------------------------------------------------!
! i        : loop index on the profil                                        !
! var      : non interpolated variable dim=33,nprofs                         !
! fvar     : altitude of meteo lvl in kilometer, = altm                      !
! fvar2    : altitude of lidar lvl in kilometer = altl                       !
! nprofs   : number of profil                                                !
!----------------------------------------------------------------------------!
! var2     : interpolated variable dim=583,nprofs                            !
!----------------------------------------------------------------------------!
! a,b      : coefficient for the interpolation equation                      !
! imol     : loop index of the meteo altitude (33)                           !
! ilid     : loop index of the lidar altitude (583)                          !
! altitude : number of level for the lidar variables                         !
! altitude2: number of level for the meteo variables                         !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : call interp(pres,pres2,altm,altl,i,it)                                !
!----------------------------------------------------------------------------!
subroutine interp(var,var2,fvar,fvar2,i,nprofs)
  implicit none

!     Indexes & parameters 
      real  ::  a,b                   
      integer  ::  i
      integer  ::  ilid                       
      integer  ::  imol                      
      integer(kind=2),parameter :: altitude2 = 33 
      integer(kind=2),parameter :: altitude = 583  
      integer  ::nprofs 
    
!     META variables
      real*4,dimension(583)  ::  fvar2  
      real*4,dimension(33)  ::   fvar 
      real,dimension(33,nprofs)  ::  var 
      real,dimension(583,nprofs)  ::  var2


    do imol=2,altitude2 
      ! exclude the nan value
      if( ( (var(imol,i).lt.1E+16).and.(var(imol,i).gt.1E+13).or.            &
         (var(imol,i).eq.-9999.) ).or.( (var(imol-1,i).lt.1E+16).and.         &
         (var(imol-1,i).gt.1E+13).or.(var(imol-1,i).eq.-9999.) ) )then

        do ilid=1,altitude
        if ((fvar2(ilid).ge.fvar(imol)).and.(fvar2(ilid).lt.fvar(imol-1)))then
            var2(ilid,i)=-9999.
        endif
        enddo

      else
      ! calculation of coefficient
         a=(var(imol,i)-var(imol-1,i))/(fvar(imol)-fvar(imol-1))
         b=var(imol,i)-a*fvar(imol)
        
        ! calculation of new interpolated variable
        do ilid=1,altitude
        if ((fvar2(ilid).ge.fvar(imol)).and.(fvar2(ilid).lt.fvar(imol-1)))then
           var2(ilid,i)=a*fvar2(ilid)+b
        endif
        enddo

      endif
    enddo

end subroutine interp
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! ***** ATB_MOL_INTERP ***** This subroutine extrapolate the atbmol from the !
!                            last atbmol known value to the Surface Elevation!
!                            in the limit of the 5 first kilometers          !
!----------------------------------------------------------------------------!
! i        : loop index on the profil                                        !
! nprofs   : number of profil                                                !
! alt      : altitude of lidar lvl in kilometer = altl                       ! 
!----------------------------------------------------------------------------!
! var3     : molecular interpolated in km-1 sr-1                             !
!----------------------------------------------------------------------------!
! l        : loop index on alt3 altitude levels                              !
! alt3     : selection of value of alt the nearest to 0 1 2 3 4 & 5km        !
! a,b      : coefficient for the interpolation equation                      !
! ilid     : loop index of the lidar altitude (583)                          !
! SeuilMol1km : treshold to add to the molecular to get the value of         !
!               molecular 1kilometer lower (this treshold is available under !
!               5km altitude)                                                !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call atb_mol_interp(mol3,altl,i,it)                                  !
!                                                                            !
!----------------------------------------------------------------------------!
!  Purpose : When the molecular value is missing close to the ground (e.g.   !
!            under 5km), this routine extrapolate this value with the tresh- !
!            old "SeuilMol1km" applied to the last known value of molecular  !
!            at 5 4 3 2 or 1 km. Then we have 2 values of molecular distant  !
!            from 1km. At this moment the routine perform the interpolation  !
!            between this 2 values at the other altitude. This process       !
!            don't stop until reaching the surface elevation.                !
!----------------------------------------------------------------------------!
subroutine atb_mol_interp(var3,alt,i,nprofs,seuil,SE)
  
  implicit none
  integer  ::  ilid, l, i, n
!  real*4,parameter  ::  SeuilMol1km = 0.00015 , SeuilTemp1km = -6.5
  real*4  :: seuil 
  real*8  ::  a,b
  integer  ::  nprofs  
  integer,dimension(10)  ::  alt3
  real,dimension(583,nprofs)  ::  var3    
  real,dimension(583) :: alt,SE



! altitude of 5 level from 0 to 5km each km (made from altl)
alt3(1)=275; alt3(2)=296; alt3(3)=329; alt3(4)=362; alt3(5)=395;  
alt3(6)=429; alt3(7)=462; alt3(8)=495; alt3(9)=529; alt3(10)=562;



if( sum(var3(1:275,i)).lt.0. )then
! Exclude entiere NaN profile
var3(:,i)=-777.

else

do ilid=275,563 

   ! Exclude the nan value
   if( (var3(ilid,i).eq.(-9999.)) .or. (var3(ilid,i).eq.(-777.)))then
     do l=1,9

      if( (ilid.gt.alt3(l)).and.(ilid.le.alt3(l+1)) ) then
           ! calculation of extrapolated value 
           var3(alt3(l+1),i)=var3(alt3(l)-1,i)+seuil
      ! Calculation of the parameter to interpolate the value between the
      ! extrapolated value and the last known value        
      a=(alt(alt3(l+1))-alt(alt3(l)-1))/(var3(alt3(l+1),i)-var3(alt3(l)-1,i))
      b=alt(alt3(l+1))-a*var3(alt3(l+1),i)
         
         ! Calculation of the new interpolated variable  
         do k=alt3(l),alt3(l+1)
            var3(k,i)=(alt(k)-b)/a
         enddo
      endif

     enddo
   endif

enddo

endif
 
endsubroutine atb_mol_interp
!----------------------------------------------------------------------------!

subroutine atb_temp_interp(var3,alt,i,nprofs,seuil,SE)
  
  implicit none
  integer  ::  ilid, l, i, n
!  real*4,parameter  ::  SeuilMol1km = 0.00015 , SeuilTemp1km = -6.5
  real*4  :: seuil 
  real*8  ::  a,b
  integer  ::  nprofs  
  integer,dimension(10)  ::  alt3
  real,dimension(583,nprofs)  ::  var3    
  real,dimension(583) :: alt,SE



! altitude of 5 level from 0 to 5km each km (made from altl)
alt3(1)=275; alt3(2)=296; alt3(3)=329; alt3(4)=362; alt3(5)=395;  
alt3(6)=429; alt3(7)=462; alt3(8)=495; alt3(9)=529; alt3(10)=562;


if ( ((var3(150,i).eq.(-9999.)).and.(var3(200,i).eq.(-9999.)).and.       &
     (var3(250,i).eq.(-9999.)).and.(var3(329,i).eq.(-9999.)).and.        &
     (var3(429,i).eq.(-9999.)).and.(var3(529,i).eq.(-9999.))) .or.       &
     ((var3(150,i).eq.(-777.)).and.(var3(200,i).eq.(-777.)).and.         &
     (var3(250,i).eq.(-777.)).and.(var3(329,i).eq.(-777.)).and.          &
     (var3(429,i).eq.(-777.)).and.(var3(529,i).eq.(-777.))) )then

!print *, i,var3(250,i),var3(495,i)
continue
else

do ilid=275,563 

   ! Exclude the nan value
   if( (var3(ilid,i).eq.(-9999.)) .or. (var3(ilid,i).eq.(-777.)))then
     do l=1,9

      if( (ilid.gt.alt3(l)).and.(ilid.le.alt3(l+1)) ) then
           ! calculation of extrapolated value 
           var3(alt3(l+1),i)=var3(alt3(l)-1,i)+seuil
      ! Calculation of the parameter to interpolate the value between the
      ! extrapolated value and the last known value        
      a=(alt(alt3(l+1))-alt(alt3(l)-1))/(var3(alt3(l+1),i)-var3(alt3(l)-1,i))
      b=alt(alt3(l+1))-a*var3(alt3(l+1),i)
         
         ! Calculation of the new interpolated variable  
         do k=alt3(l),alt3(l+1)
            var3(k,i)=(alt(k)-b)/a
         enddo
      endif

     enddo
   endif

enddo

endif
 
!if((i==49261).and.(seuil==0.00015))then
!do  ilid=1,583 
!print *, 'atb_mol_interp',var3(ilid,i)
!enddo
!endif
!print *, var3(495,i), "exit"
endsubroutine atb_temp_interp
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
! ***** ATB_MOL ***** This routine calculate the atbmol from 33lvl to 583lvl !
!                     including the calculation of the rapport               !
!----------------------------------------------------------------------------!
! alt1     : top calculation altitude level of the normalized ratio          !
! alt2     : toplow calculation altitude level of the normalized ratio       !
! nprofs   : number of profil                                                !
! i        : loop index on the profil                                        !
! var      : attenuated backscatter (=atb(583,nprof))                        !
! var2     : molecular from 33lvl to 583lvl in CN (=mol2(583,nprof))         !
!----------------------------------------------------------------------------!
! var3     : molecular in km-1 sr-1 from mol2 (=mol3(583,nprof))             !
!----------------------------------------------------------------------------!
! n        : loop index of the sliding average                               !
! ilid     : loop index of the lidar altitude (583)                          !
! matb     : sum of the atb                                                  !
! mmol     : sum of the molecular                                            !
! rapport3 : normalized ratio for each profil dim=(nprof)                    !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call atb_mol(atb,mol2,mol3,i,it,62,92)                               !
!                                                                            !
!----------------------------------------------------------------------------!
!  Purpose : In order to caculate the scattering ratio, we need to estimate  !
!            the molecular in sr-1 km-1 because :                            ! 
!            SR = atb(km-1 sr-1) / mol(km-1 sr-1)                            !
!            and we already know the atb value.                              !
!            The calculation is a sliding average of 66 profils on 20km      !
!            (1profil every 330m, 66*330m = 20km). This operation allows to  !
!            refine the molecular.                                           !
!----------------------------------------------------------------------------!
subroutine atb_mol(var,var2,var3,i,nprofs,alt1,alt2)

  implicit none
  integer  ::  ilid, i , n ,alt1, alt2       
  real  ::  matb, mmol !,rcompt  
  integer  ::  nprofs                    
  real,dimension(583,nprofs)  ::  var, var2, var3    
  real,dimension(nprofs)  ::  rapport3
      matb=0
      mmol=0
      rapport3(i)=0
    !  rcompt=0

     ! Average for the first profils
     if(i.lt.34)then  !!!!! loop on profil
        matb=0
        mmol=0 

        do n=0,65
           ! average between 22-25km or 20-25km
           do ilid=alt1,alt2   
               ! exclude the nan values
               if ((var(ilid,i).ne.(-9999)).and.(var(ilid,i).lt.1).and.      &
                   (var2(ilid,i).ne.-9999))then
                matb=matb+var(ilid,1+n)
                mmol=mmol+var2(ilid,1+n) 
              ! else
               ! rcompt=rcompt+1    !! check this process & the value rcompt
               endif
            enddo
         enddo
         ! calculation of the ratio
         rapport3(i)=matb/mmol

    ! Average for other profils
    else 
       if( (i.gt.32).and.(i.lt.(it-32)) )then
           matb=0
           mmol=0 

          do n=0,65
             ! average between 22-25km or 20-25km
             do ilid=alt1,alt2  
                ! exclude the nan values
                if ((var(ilid,i).ne.(-9999)).and.(var(ilid,i).lt.1).and.     &
                    (var2(ilid,i).ne.-9999))then
                   matb=matb+var(ilid,i-32+n)
                   mmol=mmol+var2(ilid,i-32+n) 
            !    else
              !     rcompt=rcompt+1   !! check this process & the value rcompt
                endif
             enddo
          enddo
          rapport3(i)=matb/mmol
          
    ! Average for the last profils      
          else
             matb=0
             mmol=0 
          
           do n=0,65
              ! average between 22-25km or 20-25km
              do ilid=alt1,alt2 
                 ! exclude the nan values
                 if ((var(ilid,i).ne.(-9999)).and.(var(ilid,i).lt.1).and.    &
                    (var2(ilid,i).ne.-9999))then
                   matb=matb+var(ilid,it-65+n);
                   mmol=mmol+var2(ilid,it-65+n) ;
                 endif
              enddo
           enddo
          rapport3(i)=matb/mmol

       endif

    endif !!!!! END loop on profil

! Calculation of the new molecular in km-1 sr-1   
do ilid=1,583
   if((var2(ilid,i).ne.(-9999)).and.(rapport3(i).ge.0.25e-28).and.(rapport3(i).le.0.97e-28))then
      var3(ilid,i)=rapport3(i)*var2(ilid,i)
   elseif(var2(ilid,i).eq.(-9999.))then
      var3(ilid,i)=-9999. ! allocation of nan values 
   else
    
     var3(ilid,i)=-777.
   endif
!if(i==49261)then
!print *, 'atb_mol',var3(ilid,i)
!endif
enddo
 
end subroutine atb_mol
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! **** ZERO_DETECT **** This routine detect the 0 values and substitute them !
!                        by -777                                            !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! iz     : loop index on the altitude grid                                   !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
!----------------------------------------------------------------------------!
! var    : verticaly averaged observed variables                             !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call zero_detect(pr2moy,i,iz,it,altmax)                              !
!----------------------------------------------------------------------------!
subroutine zero_detect(var,i,iz,nprofs,alt)

  implicit none
  integer ::  i,iz,nprofs,alt
  real*4,dimension(alt,nprofs)  ::  var
  
 if(var(iz,i).eq.0)then       ! if no value is detected then allocate -9999
            var(iz,i)=-9999.
 endif

end subroutine zero_detect
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! **** ZERO_DETECT_CHIM **** This routine detect the 0 values and substitute !
!                              them by -9999                                 !
!----------------------------------------------------------------------------!
! j      : loop index on latitude grid                                       !
! k      : loop index on longitude grid                                      ! 
! lat    : number of the latitude boxes                                      !
! lon    : number of the longitude boxes                                     !
! iz     : loop index on the altitude grid                                   !
! alt    : number of the altitude boxes                                      !
!----------------------------------------------------------------------------!
! var    : verticaly averaged observed variables                             !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call zero_detect_chim(pr2moy,iz,altmax,j,k,latmax,lonmax)            !
!----------------------------------------------------------------------------!
subroutine zero_detect_chim(var,iz,alt,j,k,lat,lon)

  implicit none
  integer ::  iz,alt,j,k,lat,lon
  real*4,dimension(lat,lon,alt)  ::  var
  
 if(var(j,k,iz).eq.0)then       ! if no value is detected then allocate -9999
            var(j,k,iz)=-9999
 endif

end subroutine zero_detect_chim
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! *** SR_CR_DEPOL_chim *** This routine do the SR CR & DEPOL average except  !
!                          for the -9999 values                              !
!----------------------------------------------------------------------------!
! j      : loop index on latitude grid                                       !
! k      : loop index on longitude grid                                      ! 
! lat    : number of the latitude boxes                                      !
! lon    : number of the longitude boxes                                     !
! iz     : loop index on the altitude grid                                   !
! alt    : number of the altitude boxes                                      !
! num    : index to select the method of calculation                         !
! var1, var2 : verticaly averaged observed variables                         !
! var4   : scattering ratio vertically averaged                              ! 
! ind1   : number of occurence of obs variables contained in var1            !
! ind2   : number of occurence of obs variables contained in var2            !
!----------------------------------------------------------------------------!
! var3   : result of var1/var2                                               !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call SR_CR_DEPOL_chim(pr2moy,molmoy,srmoy,indice,indicem,iz,altmax,  !
!                             j,k,latmax,lonmax,1,srmoy)                     !
!----------------------------------------------------------------------------!
subroutine SR_CR_DEPOL_chim(var1,var2,var3,ind1,ind2,iz,alt,j,k,lat,lon,num, &
                            var4)

  implicit none
  integer  ::  iz,alt,j,k,lat,lon,num
  real*4,dimension(lat,lon,alt)  ::  var1,var2,var3,var4,ind1,ind2

if(num.eq.1)then

   ! calculation if var3 = srmoy
   if ((var1(j,k,iz).eq.(-9999)).or.(var2(j,k,iz).EQ.(-9999))) then
      var3(j,k,iz) = -9999
   else
      var3(j,k,iz) = (var1(j,k,iz)*ind2(j,k,iz))/(var2(j,k,iz)*ind1(j,k,iz))
   endif

else 

    ! calculation if var3 is different than srmoy  
   if(var4(j,k,iz).eq.(-9999))then
      var3(j,k,iz) = -9999
   else
      if ((var1(j,k,iz).eq.(-9999)).or.(var2(j,k,iz).EQ.(-9999))) then    
       var3(j,k,iz) = -9999
      else
       var3(j,k,iz) = (var1(j,k,iz)*ind2(j,k,iz))/(var2(j,k,iz)*ind1(j,k,iz))
      endif
   endif

endif

end subroutine SR_CR_DEPOL_chim
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! *** SR_CR_DEPOL_mean *** This routine do the SR CR & DEPOL average except  !
!                          for the -9999 values                              !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! iz     : loop index on the altitude grid                                   !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1, var2 : verticaly averaged observed variables                         !
! ind1   : number of occurence of obs variables contained in var1            !
! ind2   : number of occurence of obs variables contained in var2            !
!----------------------------------------------------------------------------!
! var3   : result of var1/var2                                               ! 
!----------------------------------------------------------------------------!
!                                                                            !
! ex :   call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   !
!                              altmax)                                       !
!----------------------------------------------------------------------------!
subroutine SR_CR_DEPOL_mean(var1,var2,var3,ind1,ind2,i,iz,nprofs,alt)
 !  call SR_CR_DEPOL_mean(parmoy,perpmoy,depolmoy,indicep2,indicep,i,iz,it,   &
 !                              altmax)
! call SR_CR_DEPOL_mean(parmoy,perpmoy,depolmoy,indicep2,indicep,i,iz,it,   &
 !                              altmax)

  implicit none
  integer  ::  i,iz,nprofs,alt
  real,dimension(alt,nprofs)  ::  var1,var2,var3
  real,dimension(alt,nprofs)  ::  ind1,ind2

    if ((var1(iz,i).eq.(-9999.)).or.(var2(iz,i).EQ.(-9999.))) then    
     var3(iz,i) = -9999.
    elseif ((var1(iz,i).eq.(-888.)).or.(var2(iz,i).EQ.(-888.))) then 
     var3(iz,i) = -888.
    elseif ((var1(iz,i).eq.(-777.)).or.(var2(iz,i).EQ.(-777.))) then 
     var3(iz,i) = -777.  
    else
     var3(iz,i) = (var1(iz,i)*ind2(iz,i)) / (var2(iz,i)*ind1(iz,i))
    endif

end subroutine SR_CR_DEPOL_mean
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! *** FILTRE_2LVL *** This routine delete the cloud over 21km and do the     !
!                     the SR average execpt for the nan values.              !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! iz     : loop index on the altitude grid                                   !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1, var2 : verticaly averaged observed variables                         !
! ind1   : number of occurence of obs variables contained in var1            !
! ind2   : number of occurence of obs variables contained in var2            !
!----------------------------------------------------------------------------!
! var3   : result of var1/var2                                               ! 
!----------------------------------------------------------------------------!
!                                                                            !
! ex :   call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   !
!                              altmax)                                       !
!----------------------------------------------------------------------------!
!  Purpose : noise could be interprated as cloudy point in the program so,   !
!            we have built this routine to filter this phenomena over 21km   !
!            because cloud couldn't appear in this area execpt during the    !
!            polar winter.                                                   !
!----------------------------------------------------------------------------!
subroutine filtre_2lvl(var1,var2,var3,ind1,ind2,i,iz,nprofs,alt,grid)

 implicit none
 integer  ::  i,iz,nprofs,alt,lvl
 real*4,dimension(alt,nprofs)  ::  var1,var2,var3
 real*4,dimension(alt,nprofs)  ::  ind1,ind2
 character  ::  grid*8

if(gcm.eq.'LMDZ')then
   lvl=15   ! level to go past in order to apply the filter = 20.8km
elseif(gcm.eq.'LMDZ40')then
   lvl=30   !  level to go past in order to apply the filter = 21.1km
elseif(gcm.eq.'WRF')then
   lvl=46   !  level to go past in order to apply the filter = 21.7km
endif



  if (iz.ge.lvl) then

  ! filter the cloud over 21km + nan values
    if((((var1(iz,i)*ind2(iz,i))/(var2(iz,i)*ind1(iz,i))).ge.5.) .or.&
           (var1(iz,i).eq.(-9999)).or.(var2(iz,i).EQ.(-9999)).or.&
           ((var1(iz,i).eq.(-777)).or.(var2(iz,i).EQ.(-777))))then
       var3(iz,i) = -777       
    else
       ! calculation of the average
       var3(iz,i) = (var1(iz,i)*ind2(iz,i))/(var2(iz,i)*ind1(iz,i))
    endif

  elseif((var1(iz,i).eq.(-9999)).or.(var2(iz,i).EQ.(-9999))) then 
! filter the nan values       
     var3(iz,i) = -9999
  elseif ((var1(iz,i).eq.(-888)).or.(var2(iz,i).EQ.(-888))) then 
     var3(iz,i) = -888
  elseif ((var1(iz,i).eq.(-777)).or.(var2(iz,i).EQ.(-777))) then 
     var3(iz,i) = -777       
  else
        ! calculation of the average
     var3(iz,i) = (var1(iz,i)*ind2(iz,i))/(var2(iz,i)*ind1(iz,i))
  endif

   
end subroutine filtre_2lvl
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** FRACTION_SUBGRID2 *** This routine calculate the cloudy,clear,         !
!                           saturated... fraction flag(0/1)                  !
!----------------------------------------------------------------------------!
! nprofs   : number of profil                                                !
! alt      : number of the altitude boxes                                    !
! i        : loop index on the profil                                        !
! switch       : switch to chose cloudy or sat mode                          !
! var          : vertically averaged SR                                      !
! var1, var2   : verticaly averaged observed variables (atb & mol moy)       !
! ind1   : number of occurence of obs variables contained in var1            !
! ind2   : number of occurence of obs variables contained in var2            ! 
!----------------------------------------------------------------------------!
! frac1  : fully attenuated fraction                                         !
! frac2  : cloudy fraction                                                   !
! frac3  : clear fraction                                                    !
! frac4  : uncertain fraction                                                !
! frac5  : NaN fraction                                                      !
! frac6  : Surface_elevation fraction                                        !
!----------------------------------------------------------------------------!
! iz       : loop index on the altitude grid                                 !
! fracttot : sum of all the fraction (=1)                                    !
! toplvlcloud  : first cloudy point encountered from the top to the ground   !
! toplvlsat    : first fully attenuated point encountered from the top to    !
!                the ground                                                  ! 
! SeuilClearSr : clear treshold detection = 1.2                              !
! SeuilSatSr   : fully attenuated treshold detection = 0.01                  !
! SeuilSrCloud : cloudy treshold detection = 5                               !
! SeuilDeltAtb : delta atb treshold = 1.4e-03 km-1 sr-1                      !
! delta        : atb - mol                                                   !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : call fraction_subgrid2(srmoy,pr2moy,indice,molmoy,indicem,satfraction,! 
!                             cloudfraction,clearfraction,uncertfraction,    ! 
!                             nanfraction,sefraction,i,altmax,it,switch2)    !
!----------------------------------------------------------------------------!
 
subroutine fraction_subgrid2(var,var1,ind1,var2,ind2,frac1,frac2,frac3,frac4,&
                             frac5,frac6,frac7,i,alt,nprofs,switch)!frac7,&
                          !   frac8,i,alt,nprofs)!frac1bis,frac2bis

  implicit none
  integer  ::  i,iz,nprofs,alt
  real*4  ::  fracttot
  real  ::  toplvlcloud,toplvlsat
  real, parameter  ::  SeuilClearSr = 1.2  
  real, parameter  ::   SeuilSatSr = 0.01  
  real, parameter  ::    SeuilSrCloud = 5.    
  real,parameter  ::    SeuilDeltAtb = 2.5e-03   
  real*4  ::  delta    ! delta atb = atb-atbmol
  character  ::  switch*6
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac1,frac2,frac3,frac4,frac5,frac6,frac7!,frac8!,frac2bis,frac1bis


delta = 0

toplvlcloud=0; toplvlsat=0

B1 : do iz=alt,1,-1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2



  do iz=1,alt
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! pts nuageux
               endif
elseif(trim(switch).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).gt.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud).and.(var(iz,i).gt.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud).and.   &
            (delta.lt.SeuilDeltAtb))) then
                  frac3(iz,i)=frac3(iz,i)+1   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud).and.(var(iz,i).gt.SeuilClearSr).and.&
             (delta.gt.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif
  enddo   



 do  iz=1,alt  
 fracttot=frac1(iz,i)+frac2(iz,i)+frac3(iz,i)+frac4(iz,i)+frac5(iz,i)+       &
          frac6(iz,i)+frac7(iz,i) !check = 1
  if (fracttot.ne.1) then
     print *, "fraction error"
     print *, 'instant',fracttot,frac1(iz,i),frac2(iz,i),frac3(iz,i),        &
                        frac4(iz,i),frac5(iz,i),frac6(iz,i),frac7(iz,i)
  endif
enddo
end subroutine fraction_subgrid2
!----------------------------------------------------------------------------!

subroutine fraction_subgrid2_8km(seuilsnrlow,seuilsnrhigh,var,var1,ind1,var2,ind2,frac1,frac2,frac3,frac4,&
                             frac5,frac6,frac7,i,alt,nprofs,toplow,topmid,switch,switch2)!frac7,&
                          !   frac8,i,alt,nprofs)!frac1bis,frac2bis

  implicit none
  integer  ::  i,iz,nprofs,alt,toplow,topmid,seuilsnrhigh,seuilsnrlow
  real*4  ::  fracttot
  real  ::  toplvlcloud,toplvlsat
  real, parameter  ::  SeuilStrat = 30.  
  real, parameter  ::  SeuilClearSr = 1.2  
  real, parameter  ::   SeuilSatSr = 0.01  
  real, parameter  ::    SeuilSrCloud = 5.
  real, parameter  ::    SeuilSrCloud2 = 15
  real  ::  SeuilSrCloud3
  real,parameter  ::    SeuilDeltAtb = 2.5e-03   
  real*4  ::  delta    ! delta atb = atb-atbmol
  character  ::  switch*5,switch2*6
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac1,frac2,frac3,frac4,frac5,frac6,frac7!,frac8!,frac2bis,frac1bis


delta = 0

toplvlcloud=0; toplvlsat=0

B1 : do iz=alt,1,-1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2


!print *, 'switch2',trim(switch2)


if(trim(switch).eq.'day')then

SeuilSrCloud3 = SeuilSrCloud

!print *, 'test'


do iz=seuilsnrhigh,alt

     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif

          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ))then
                  frac3(iz,i)=frac3(iz,i)+1              
          endif

          if((var(iz,i).ge.SeuilSrCloud3).and.(delta.lt.SeuilDeltAtb)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -777 (surface elevation) 
          endif 

          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif          
enddo   

do iz=1,seuilsnrlow
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ))then
                  frac3(iz,i)=frac3(iz,i)+1              
          endif

          if((var(iz,i).ge.SeuilSrCloud3).and.(delta.lt.SeuilDeltAtb)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -777 (surface elevation) 
          endif

           if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif     
enddo   


B3:   do iz=1,toplow-1
      if(var(iz,i).gt.SeuilStrat)then
      SeuilSrCloud3 = SeuilSrCloud2 
      exit B3
      endif
   enddo B3

do iz=seuilsnrlow+1,seuilsnrhigh-1
 
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ))then
                  frac3(iz,i)=frac3(iz,i)+1              
          endif

          if((var(iz,i).ge.SeuilSrCloud3).and.(delta.lt.SeuilDeltAtb)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -777 (surface elevation) 
          endif

           if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif
enddo   


!
elseif(trim(switch).eq.'night')then


do iz=1,alt

SeuilSrCloud3 = SeuilSrCloud
     
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -777 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ))then
                  frac3(iz,i)=frac3(iz,i)+1              
          endif

          if((var(iz,i).ge.SeuilSrCloud3).and.(delta.lt.SeuilDeltAtb)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -777 (surface elevation) 
          endif

          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif
enddo   


endif

!print *, frac1(iz,i),frac2(iz,i),frac3(iz,i),        &
!                        frac4(iz,i),frac5(iz,i),frac6(iz,i),frac7(iz,i)

 do  iz=1,alt  
 fracttot=frac1(iz,i)+frac2(iz,i)+frac3(iz,i)+frac4(iz,i)+frac5(iz,i)+       &
          frac6(iz,i)+frac7(iz,i) !check = 1
  if (fracttot.ne.1) then
     print *, "fraction error"
     print *, 'instant',fracttot,frac1(iz,i),frac2(iz,i),frac3(iz,i),        &
                        frac4(iz,i),frac5(iz,i),frac6(iz,i),frac7(iz,i), i,iz,var(iz,i)
!     stop
  endif
enddo


end subroutine fraction_subgrid2_8km
!----------------------------------------------------------------------------!


subroutine fraction_subgrid2_8km_delta(var,var1,ind1,var2,ind2,frac1,frac2,frac3,frac4,&
                             frac5,frac6,frac7,i,alt,nprofs,toplow,topmid,switch,switch2)!frac7,&
                          !   frac8,i,alt,nprofs)!frac1bis,frac2bis

  implicit none
  integer  ::  i,iz,nprofs,alt,toplow,topmid
  real*4  ::  fracttot
  real  ::  toplvlcloud,toplvlsat
  real, parameter  ::  SeuilStrat = 30.  
  real, parameter  ::  SeuilClearSr = 1.2  
  real, parameter  ::   SeuilSatSr = 0.01  
  real, parameter  ::    SeuilSrCloud = 5.
  real, parameter  ::    SeuilSrCloud2 = 15
  real  ::  SeuilSrCloud3
  real,parameter  ::    SeuilDeltAtb = 2.5e-03   
  real*4  ::  delta    ! delta atb = atb-atbmol
  character  ::  switch*5,switch2*6
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac1,frac2,frac3,frac4,frac5,frac6,frac7!,frac8!,frac2bis,frac1bis


delta = 0

toplvlcloud=0; toplvlsat=0

B1 : do iz=alt,1,-1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2


!print *, 'switch2',trim(switch2)


if(trim(switch).eq.'day')then

SeuilSrCloud3 = SeuilSrCloud

!print *, 'test'


do iz=18,alt

     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3))then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then!.or. &
         !   ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
         !   (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
         !   (delta.lt.SeuilDeltAtb))) then
                  frac3(iz,i)=frac3(iz,i)+1   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then!.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif          
enddo   

do iz=1,5
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac3(iz,i)=frac3(iz,i)+1   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then!.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif     
enddo   


B3:   do iz=1,toplow-1
      if(var(iz,i).gt.SeuilStrat)then
      SeuilSrCloud3 = SeuilSrCloud2 
      exit B3
      endif
   enddo B3

do iz=6,17
 
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac3(iz,i)=frac3(iz,i)+1   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then !.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif
enddo   


!
elseif(trim(switch).eq.'night')then


do iz=1,alt

SeuilSrCloud3 = SeuilSrCloud
     
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac3(iz,i)=frac3(iz,i)+1   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then !.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif
enddo   


endif

!print *, frac1(iz,i),frac2(iz,i),frac3(iz,i),        &
!                        frac4(iz,i),frac5(iz,i),frac6(iz,i),frac7(iz,i)

 do  iz=1,alt  
 fracttot=frac1(iz,i)+frac2(iz,i)+frac3(iz,i)+frac4(iz,i)+frac5(iz,i)+       &
          frac6(iz,i)+frac7(iz,i) !check = 1
  if (fracttot.ne.1) then
     print *, "fraction error"
     print *, 'instant',fracttot,frac1(iz,i),frac2(iz,i),frac3(iz,i),        &
                        frac4(iz,i),frac5(iz,i),frac6(iz,i),frac7(iz,i), i,iz,var(iz,i)
!     stop
  endif
enddo


end subroutine fraction_subgrid2_8km_delta
!----------------------------------------------------------------------------!


subroutine fraction_subgrid3_8km(seuilsnrlow,seuilsnrhigh,altvar,var,var1,ind1,var2,ind2,frac,i,alt,nprofs,switch,switch2)

  implicit none
  integer  ::  i,iz,nprofs,alt,toplow,topmid,seuilsnrhigh,seuilsnrlow
  real*4  ::  fracttot
  real  ::  toplvlcloud,toplvlsat
  real, parameter  ::  SeuilStrat = 30.  
  real, parameter  ::  SeuilClearSr = 1.2  
  real, parameter  ::   SeuilSatSr = 0.01  
  real, parameter  ::    SeuilSrCloud = 5.
  real, parameter  ::    SeuilSrCloud2 = 15
  real  ::  SeuilSrCloud3
  real,parameter  ::    SeuilDeltAtb = 2.5e-03   
  real*4  ::  delta    ! delta atb = atb-atbmol
  character  ::  switch*5,switch2*6
  real*4,dimension(alt+1)  ::  altvar
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac

! uncert+4 nan+1 surf+6 rejec+7 clear+2 cld+3 sat+8


delta = 0
toplvlcloud=0; toplvlsat=0



B1 : do iz=alt,1,-1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2


!print *, 'switch2',trim(switch2)


if(trim(switch).eq.'day')then

SeuilSrCloud3 = SeuilSrCloud

!print *, 'test'

! uncert+4 nan+1 surf+6 rejec+7 clear+2 cld+3 sat+8

do iz=seuilsnrhigh,alt   !=18 for cfmip2

     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
            (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif          
enddo   


do iz=1,seuilsnrlow  ! =5 for cfmip2
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
            (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif     
enddo   


B3:   do iz=1,toplow-1
      if(var(iz,i).gt.SeuilStrat)then
      SeuilSrCloud3 = SeuilSrCloud2 
      exit B3
      endif
   enddo B3

do iz=seuilsnrlow+1,seuilsnrhigh-1
 
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
            (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif
enddo   


!
elseif(trim(switch).eq.'night')then


do iz=1,alt

SeuilSrCloud3 = SeuilSrCloud
     
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -777 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ))then
                  frac(iz,i)=frac(iz,i)+2   
          endif
          if((var(iz,i).ge.SeuilSrCloud3).and.(delta.lt.SeuilDeltAtb)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -777 (surface elevation)            
          endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif
enddo   


endif


end subroutine fraction_subgrid3_8km
!----------------------------------------------------------------------------!


subroutine fraction_subgrid3_8km_delta(var,var1,ind1,var2,ind2,frac,i,alt,nprofs,switch,switch2)

  implicit none
  integer  ::  i,iz,nprofs,alt,toplow,topmid
  real*4  ::  fracttot
  real  ::  toplvlcloud,toplvlsat
  real, parameter  ::  SeuilStrat = 30.  
  real, parameter  ::  SeuilClearSr = 1.2  
  real, parameter  ::   SeuilSatSr = 0.01  
  real, parameter  ::    SeuilSrCloud = 5.
  real, parameter  ::    SeuilSrCloud2 = 15
  real  ::  SeuilSrCloud3
  real,parameter  ::    SeuilDeltAtb = 2.5e-03   
  real*4  ::  delta    ! delta atb = atb-atbmol
  character  ::  switch*5,switch2*6
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac

! uncert+4 nan+1 surf+6 rejec+7 clear+2 cld+3 sat+8


delta = 0

toplvlcloud=0; toplvlsat=0

B1 : do iz=alt,1,-1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2


!print *, 'switch2',trim(switch2)


if(trim(switch).eq.'day')then

SeuilSrCloud3 = SeuilSrCloud

!print *, 'test'

! uncert+4 nan+1 surf+6 rejec+7 clear+2 cld+3 sat+8

do iz=18,alt

     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then !.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif          
enddo   


do iz=1,5
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then !.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif     
enddo   


B3:   do iz=1,toplow-1
      if(var(iz,i).gt.SeuilStrat)then
      SeuilSrCloud3 = SeuilSrCloud2 
      exit B3
      endif
   enddo B3

do iz=6,17
 
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then !.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif
enddo   


!
elseif(trim(switch).eq.'night')then


do iz=1,alt

SeuilSrCloud3 = SeuilSrCloud
     
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then!.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then!.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif
enddo   


endif


end subroutine fraction_subgrid3_8km_delta
!----------------------------------------------------------------------------!








!----------------------------------------------------------------------------!
! *** FRACTION_SUBGRID3 *** This routine calculate the cloudy,clear,         !
!                           saturated... fraction in order to be used by     !
!                           matlab (flag 1/2/3/4/5/6)                        !
!----------------------------------------------------------------------------!
! i        : loop index on the profil                                        !
! nprofs   : number of profil                                                !
! alt      : number of the altitude boxes                                    !
! var1, var2   : verticaly averaged observed variables (atb & mol moy)       !
! var    : vertically averaged SR                                            ! 
! ind1   : number of occurence of obs variables contained in var1            !
! ind2   : number of occurence of obs variables contained in var2            ! 
!----------------------------------------------------------------------------!
! frac   : fraction of fully, cloudy, clear, uncert, nan, se point           !
!          cloud=6, clear=2, fully_att=1, uncert=3, nan=4, SE=5              !
!----------------------------------------------------------------------------!
! iz       : loop index on the altitude grid                                 !
! toplvlcloud  : first cloudy point encountered from the top to the ground   !
! toplvlsat    : first fully attenuated point encountered from the top to    !
!                the ground                                                  ! 
! SeuilClearSr : clear treshold detection = 1.2                              !
! SeuilSatSr   : fully attenuated treshold detection = 0.01                  !
! SeuilSrCloud : cloudy treshold detection = 3                               !
! SeuilDeltAtb : delta atb treshold = 1.4e-03 km-1 sr-1                      !
! delta        : atb - mol                                                   !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : call fraction_subgrid3(srmoy,pr2moy,indice,molmoy,indicem,fractot,i,  !
!                             altmax,it)                                     !
!----------------------------------------------------------------------------!
subroutine fraction_subgrid3(var,var1,ind1,var2,ind2,frac,i,alt,nprofs,switch)

  implicit none
  integer  ::  i,iz,nprofs,alt
  character  ::  switch*6

  real, parameter  ::  SeuilClearSr = 1.2      ! Seuil de saturation en atb en attendant flag
  real, parameter  ::   SeuilSatSr = 0.01     ! seuil de saturation en ATB=valeur ATBmol à 30km
  real, parameter  ::    SeuilSrCloud = 5.    ! seuil detection nuageuse sr
  real,parameter  ::    SeuilCrCloud = 0.   ! seuil detection nuageuse  cr>0.6, ne filtre pas les gros aerosols (atbd calipso)
  real,parameter  ::    SeuilDeltAtb = 2.5e-03    ! seuil détection unclassify
  real*4  ::  delta    ! delta atb = atb-atbmol
  integer  ::  toplvlcloud,toplvlsat
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac!,frac7

!flag 1= satfraction, 2=clear, 3=uncert, 5=nan, 4=SE, 6=cloud



delta = 0

toplvlcloud=0; toplvlsat=0

B1 : do iz=alt,1,-1
   if((var1(iz,i).eq.(-9999)).or.(var(iz,i).eq.(-888)).or.(var(iz,i).eq.(-777))) exit B1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
   if((var1(iz,i).eq.(-9999)).or.(var(iz,i).eq.(-888)).or.(var(iz,i).eq.(-777))) exit B2

       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2



  do iz=1,alt
   if((var1(iz,i).eq.(-9999)).or.(var(iz,i).eq.(-888)).or.(var(iz,i).eq.(-777)))then
   delta=0
   else
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
   endif
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! pts nuageux
               endif
elseif(trim(switch).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).gt.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud).and.(var(iz,i).gt.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud).and.   &
            (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud).and.(var(iz,i).gt.SeuilClearSr).and.&
             (delta.gt.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif

!!$if((i.ge.35618).and.(i.le.35629))then
!!$print *,i
!!$if(iz.lt.4)then
!!$    print *, iz,i, delta,toplvlcloud,toplvlsat
!!$    print *, var(iz,i),frac(iz,i)
!!$endif
!!$
!!$endif


  enddo   
         

end subroutine fraction_subgrid3
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** VERTICAL_MEAN_CHIM *** This routine calculate the vertical mean for    !
!                             the CHIMERE grid                               !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! ilid   : loop index of the lidar altitude (583)                            !
! nprofs : number of profil                                                  !
! j      : loop index on latitude grid                                       !
! k      : loop index on longitude grid                                      ! 
! lat    : number of the latitude boxes                                      !
! lon    : number of the longitude boxes                                     !
! iz     : loop index on the altitude grid                                   !
! alt    : number of the altitude boxes                                      !
! nvar   : switch to chose cloudy or sat mode                                !
! var3   : attenuated backscatter                                            !
! var1   : observed variable to average                                      !
! ind    : number of occurence of obs variables contained in var1            !
!----------------------------------------------------------------------------!
! var2   : averaged variable                                                 !
!----------------------------------------------------------------------------!
!                                                                            !
! ex  call vertical_mean_chim(atb,pr2moy,atb,indice,i,iz,ilid,it,altmax,     !
!                             altitude,j,k,latmax,lonmax,1)                  !
!----------------------------------------------------------------------------!
subroutine vertical_mean_chim(var1,var2,var3,ind,i,iz,ilid,nprofs,alt,alt2,j,&
                              k,lat,lon,nvar) 

  implicit none
  integer  ::  nvar
  integer  ::  i,nprofs,j,k,iz,alt,ilid,lat,lon
  integer(kind=2)  ::  alt2
  real*4,dimension(alt2,nprofs)  ::  var1,var3
  real*4,dimension(lat,lon,alt)  ::  var2
  real*4,dimension(lat,lon,alt)  ::  ind
  

  if (nvar==1) then
     ! Vertical average for atb
     if ( (var1(ilid,i).ne.(-9999)).and.(var1(ilid,i).lt.1) ) then
        var2(j,k,iz)=var2(j,k,iz)+var1(ilid,i)
        ind(j,k,iz)=ind(j,k,iz)+1 
     endif

  else
      ! Vertical average for perp
     if(nvar==2) then
        if ( (var1(ilid,i).ne.(-9999)).and.(var3(ilid,i).ne.(-9999)) .and.   &
             (var3(ilid,i).lt.1) )then 
                var2(j,k,iz)=var2(j,k,iz)+(var3(ilid,i)-var1(ilid,i))
                ind(j,k,iz)=ind(j,k,iz)+1 
        endif
     else
      ! Other vertical average
        if ( var1(ilid,i).ne.(-9999) ) then                            
           var2(j,k,iz)=var2(j,k,iz)+var1(ilid,i)
           ind(j,k,iz)=ind(j,k,iz)+1 
        endif
     endif

  endif
                    
end subroutine vertical_mean_chim
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** VERTICAL_MEAN *** This routine calculate the vertical mean for other   !
!                       grid                                                 !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! ilid   : loop index of the lidar altitude (583)                            !
! nprofs : number of profil                                                  !
! iz     : loop index on the altitude grid                                   !
! alt    : number of the altitude boxes                                      !
! alt2   : number of level of the lidar variables                            !
! nvar   : switch to chose cloudy or sat mode                                !
! var3   : attenuated backscatter                                            !
! var1   : observed variable to average                                      !
! ind    : number of occurence of obs variables contained in var1            ! 
!----------------------------------------------------------------------------!
! var2   : averaged variable                                                 !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call vertical_mean(atb,pr2moy,atb,indice,i,iz,ilid,it,altmax,        !
!                          altitude,1)                                       !
!----------------------------------------------------------------------------!
subroutine vertical_mean(var1,var2,var3,ind,i,iz,ilid,nprofs,alt,alt2,nvar)
!             call vertical_mean(temp2,tempmoy,atb,indicetemp,i,iz,ilid,it,altmax,&
!                                 altitude,3)    
  
  implicit none
  integer  ::  nvar
  integer  ::  i,iz,nprofs,alt,ilid
  integer(kind=2)  ::  alt2
  real*4,dimension(alt2,nprofs)  ::  var1,var3
  real*4,dimension(alt,nprofs)  ::  var2
  real*4,dimension(alt,nprofs)  ::  ind
  

  if (nvar==1) then ! Vertical average for atb

     if ( (var1(ilid,i).ne.(-9999.)).and.(var1(ilid,i).lt.1).and.(var1(ilid,i).ne.(-888.)).and.(var1(ilid,i).ne.(-777.))) then         
        var2(iz,i)=var2(iz,i)+var1(ilid,i)
        ind(iz,i)=ind(iz,i)+1 
     endif

  elseif(nvar==2) then ! Vertical average for perp

     if ( (var1(ilid,i).ne.(-9999.)).and.(var3(ilid,i).ne.(-9999.)) .and.   &
          (var3(ilid,i).lt.1).and.(var1(ilid,i).ne.(-888.)) .and.  &
          (var1(ilid,i).ne.(-777.)) )then 
         var2(iz,i)=var2(iz,i)+(var3(ilid,i)-var1(ilid,i))
         ind(iz,i)=ind(iz,i)+1 
     endif

  elseif(nvar==3) then ! Other vertical average 
   
     if (( var1(ilid,i).ne.(-9999.).and.(var1(ilid,i).ne.(-888.)) ) .and.   &
        (var3(ilid,i).lt.1).and.(var3(ilid,i).ne.(-9999.)).and.  &
        (var1(ilid,i).ne.(-777.)) ) then                            
           var2(iz,i)=var2(iz,i)+var1(ilid,i)
           ind(iz,i)=ind(iz,i)+1 
     endif

  elseif(nvar==4) then ! TEMP vertical average 
  
     if (( var1(ilid,i).ne.(-9999.).and.(var1(ilid,i).ne.(-888.)) ) .and.   &
        (var1(ilid,i).ne.(-777.)) ) then                            
           var2(iz,i)=var2(iz,i)+var1(ilid,i)
           ind(iz,i)=ind(iz,i)+1 
     endif


  endif

                    
end subroutine vertical_mean
!----------------------------------------------------------------------------!


subroutine vertical_mean_hori(var1,var2,var3,ind,i,iz,ilid,nprofs,alt,alt2,nvar,profavg,profmax)
  
  implicit none
  integer  ::  nvar, n
  integer  ::  i,iz,nprofs,alt,ilid,profavg,profmax
  integer(kind=2)  ::  alt2
  real*4  ::  var2m, indm
  real*4,dimension(alt2,nprofs)  ::  var1,var3
  real*4,dimension(alt,nprofs)  ::  var2
  real*4,dimension(alt,nprofs)  ::  ind
  



  if (nvar==1) then
  
    if (i.le.profavg)then
       var2m=0
       indm=0
      do n = 0,profmax   
       ! Vertical average for atb
         if ( (var1(ilid,1+n).ne.(-9999)).and.(var1(ilid,1+n).lt.1) ) then         
            var2m=var2m+var1(ilid,1+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif


    elseif ( (i.gt.profavg).and.(i.le.(it-(profavg+1))) ) then   
       var2m=0
       indm=0
      do n = 0,profmax   
         ! Vertical average for atb
         if ( (var1(ilid,i-profavg+n).ne.(-9999)).and.(var1(ilid,i-profavg+n).lt.1) ) then         
            var2m=var2m+var1(ilid,i-profavg+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif

    else
       var2m=0
       indm=0
      do n = 0,profmax   
         ! Vertical average for atb
         if ( (var1(ilid,it-profmax+n).ne.(-9999)).and.(var1(ilid,it-profmax+n).lt.1) ) then         
            var2m=var2m+var1(ilid,it-profmax+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif

    endif
  

  else

     ! Vertical average for perp
     if(nvar==2) then
        if ( (var1(ilid,i).ne.(-9999)).and.(var3(ilid,i).ne.(-9999)) .and.   &
             (var3(ilid,i).lt.1) )then 
               var2(iz,i)=var2(iz,i)+(var3(ilid,i)-var1(ilid,i))
               ind(iz,i)=ind(iz,i)+1 
        endif
     else

    if (i.le.profavg)then
       var2m=0
       indm=0
      do n = 0,profmax   
       ! Vertical average for atb
         if ( (var1(ilid,1+n).ne.(-9999)) ) then         
            var2m=var2m+var1(ilid,1+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif


    elseif ( (i.gt.profavg).and.(i.le.(it-(profavg+1))) ) then   
       var2m=0
       indm=0
      do n = 0,profmax   
         ! Vertical average for atb
         if ( (var1(ilid,i-profavg+n).ne.(-9999)) ) then         
            var2m=var2m+var1(ilid,i-profavg+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif

    else
       var2m=0
       indm=0
      do n = 0,profmax   
         ! Vertical average for atb
         if ( (var1(ilid,it-profmax+n).ne.(-9999)) ) then         
            var2m=var2m+var1(ilid,it-profmax+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif

    endif
  



     endif

  endif
                  
 

                    
end subroutine vertical_mean_hori
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** SE_KM_2_PRES *** This subroutine convert the SE from km to hPa and     !
!                      allocate -888 to the variable below this threshold    !
!----------------------------------------------------------------------------!
! nprofs : number of profil                                                  !
! i      : loop index on the profil                                          !
! alt    : number of the altitude boxes                                      !
! var1   : Surface_elevation (nprofs)                                        !
! var3   : altitude of lidar lvl in kilometer, (583)                         !
! var4   : Pressure interpolated from 33lvl to 583lvl in hPa (583,nprofs)    !
! var5   : pressure of model boxes                                           !
!----------------------------------------------------------------------------!
! var2   : Surface_elevation in hPa  (nprofs)                                !
! var6   : add the surface_elevation to this variable                        !
!----------------------------------------------------------------------------!
! iz     : loop index on the altitude grid                                   !
! j      : loop index on Surface Elevation level                             !
! ilid   : loop index of the lidar altitude (583)                            !
! alt2   : number of level of the lidar variables                            !
! a,b    : coefficient for the interpolation equation                        !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : call SE_km_2_Pres(SE,SEp,altl,pres2,prestop,srmoy,altmax,it,i)        !
!----------------------------------------------------------------------------!
subroutine SE_km_2_Pres(var1,var2,var3,var4,var5,var6,alt,nprofs,i)

implicit none
integer  ::  ilid, iz, alt, nprofs,i,j
integer,parameter  ::  altitude=583
real*4  ::  a,b
real*4,dimension(altitude,nprofs)  ::  var4
real*4,dimension(alt)  ::  var5
real,dimension(nprofs)  ::  var1,var2
real*4,dimension(alt,nprofs)  ::  var6
real*4,dimension(altitude)  ::  var3
! var1 = SE, var2 = SEp, var3 = altl, var4 = pres2, var5 = prestop, var6 = srmoy
do ilid=1,altitude
  ! print *, var1(1,i), var3(ilid)
      if( var1(i).eq.var3(ilid) )then
         var2(i)=var4(ilid,i)
      else
        if( var1(i).eq.var3(ilid+1) )then
             var2(i)=var4(ilid+1,i)
        else
          if((var1(i).lt.var3(ilid)).and.(var1(i).gt.var3(ilid+1)))then
             if( var4(ilid,i).ne.-9999 )then
               a=(var4(ilid,i)-var4(ilid+1,i))/(var3(ilid)-var3(ilid+1))
               b=var4(ilid,i)-a*var3(ilid)
               var2(i)=a*var1(i)+b
             else
               var2(i)=-9999
             endif
          endif
        endif
      endif
   enddo

   do iz=alt,2,-1 
         if ( (var2(i).gt.var5(iz)).and.(var2(i).lt.var5(iz-1)) )then
            do j=1,iz
               var6(j,i)=-888
            enddo
         endif
   enddo

end subroutine SE_km_2_Pres
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** SE_KM_2_PRES *** This subroutine convert the SE from km to hPa and     !
!                      allocate -9999 to the variable below this threshold   !
!----------------------------------------------------------------------------!
! nprofs : number of profil                                                  !
! i      : loop index on the profil                                          !
! alt    : number of the altitude boxes                                      !
! var1   : Surface_elevation (nprofs)                                        !
! var3   : altitude of lidar lvl in kilometer, (583)                         !
! var4   : Pressure interpolated from 33lvl to 583lvl in hPa (583,nprofs)    !
! var5   : pressure of model boxes                                           !
!----------------------------------------------------------------------------!
! var2   : Surface_elevation in hPa  (nprofs)                                !
! var6   : add the surface_elevation to this variable                        !
!----------------------------------------------------------------------------!
! iz     : loop index on the altitude grid                                   !
! j      : loop index on Surface Elevation level                             !
! ilid   : loop index of the lidar altitude (583)                            !
! alt2   : number of level of the lidar variables                            !
! a,b    : coefficient for the interpolation equation                        !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : call SE_km_2_Pres(SE,SEp,altl,pres2,prestop,srmoy,altmax,it,i)        !
!----------------------------------------------------------------------------!
subroutine SE_km_2_Pres2(var1,var2,var3,var4,var5,var6,alt,nprofs,i)

implicit none
integer  ::  ilid, iz, alt, nprofs,i,j
integer,parameter  ::  altitude=583
real*4  ::  a,b
real*4,dimension(altitude,nprofs)  ::  var4
real*4,dimension(alt)  ::  var5
real,dimension(nprofs)  ::  var1,var2
real*4,dimension(alt,nprofs)  ::  var6
real*4,dimension(altitude)  ::  var3
! var1 = SE, var2 = SEp, var3 = altl, var4 = pres2, var5 = prestop, var6 = srmoy
do ilid=1,altitude
  ! print *, var1(1,i), var3(ilid)
      if( var1(i).eq.var3(ilid) )then
         var2(i)=var4(ilid,i)
      else
        if( var1(i).eq.var3(ilid+1) )then
             var2(i)=var4(ilid+1,i)
        else
          if((var1(i).lt.var3(ilid)).and.(var1(i).gt.var3(ilid+1)))then
             if( var4(ilid,i).ne.-9999 )then
               a=(var4(ilid,i)-var4(ilid+1,i))/(var3(ilid)-var3(ilid+1))
               b=var4(ilid,i)-a*var3(ilid)
               var2(i)=a*var1(i)+b
             else
               var2(i)=-9999
             endif
          endif
        endif
      endif
   enddo

   do iz=alt,2,-1 
         if ( (var2(i).gt.var5(iz)).and.(var2(i).lt.var5(iz-1)) )then
            do j=1,iz
               var6(j,i)=-9999
            enddo
         endif
   enddo

end subroutine SE_km_2_Pres2
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** SURF_DETECT *** This subroutine detect the surface elevation and       !
!                     allocate -888 to the variable below this threshold     !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1   : Surface Elevation                                                 !
! var2   : altitude of model boxes                                           !
!----------------------------------------------------------------------------!
! var3   : varically averaged variable                                       !
!----------------------------------------------------------------------------!
! j      : loop index on Surface Elevation level                             !
! iz     : loop index on the altitude grid                                   !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call Surf_detect(SE,altmod,pr2moy,altmax,it,i)                       !
!----------------------------------------------------------------------------!
subroutine Surf_detect(var1,var2,var3,alt,nprofs,i)

  implicit none

  integer  ::  iz, alt, nprofs, i,j
  real*4,dimension(alt+1)  ::  var2
  real,dimension(nprofs)  ::  var1
  real*4,dimension(alt,nprofs)  ::  var3

  B1: do iz=1,alt-1 
         if ( (var1(i).gt.var2(iz)).and.(var1(i).le.var2(iz+1)) )then
          B2:   do j=1,iz
                   var3(j,i)=-888
                enddo  B2
         exit B1            
         endif 
  
      enddo B1              
   
end subroutine Surf_detect
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** SURF_DETECT2 *** This subroutine detect the surface elevation and      !
!                      allocate -9999 to the variable below this threshold   !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1   : Surface Elevation                                                 !
! var2   : altitude of model boxes                                           !
!----------------------------------------------------------------------------!
! var3   : varically averaged variable                                       !
!----------------------------------------------------------------------------!
! j      : loop index on Surface Elevation level                             !
! iz     : loop index on the altitude grid                                   !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call Surf_detect2(SE,altmod,pr2moy,altmax,it,i)                      !
!----------------------------------------------------------------------------!
subroutine Surf_detect2(var1,var2,var3,alt,nprofs,i,mod)

  implicit none

  integer  ::  iz, alt, nprofs, i,j
  character  ::  mod*8
  real*4,dimension(alt+1)  ::  var2
  real,dimension(nprofs)  ::  var1
  real*4,dimension(alt,nprofs)  ::  var3

if(mod=='altitude')then
 B11:  do iz=1,alt-1 
         if ( (var1(i).gt.var2(iz)).and.(var1(i).le.var2(iz+1)) )then
            B22: do j=1,iz
                   var3(j,i)=-888
                enddo B22
                exit B11
         endif
      enddo B11             
elseif(mod=='pressure')then
  B12:  do iz=1,alt-1 
         if ( (var1(i).lt.var2(iz)).and.(var1(i).gt.var2(iz+1)) )then
         B21:    do j=1,iz
                   var3(j,i)=-888
                enddo B21
                exit B12
         endif  
      enddo  B12 
endif
  
end subroutine Surf_detect2
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** SE_ALT_CHIM *** This subroutine detect the surface elevation and       !
!                     allocate -9999 to the variable below this threshold for!
!                     the CHIMERE grid                                       !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1   : Surface Elevation                                                 !
! var2   : altitude of model boxes                                           !
!----------------------------------------------------------------------------!
! var3   : varically averaged variable                                       !
!----------------------------------------------------------------------------!
! l      : loop index on Surface Elevation level                             !
! ilid   : loop index on the lidar altitude (583)                            !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call SE_alt_atb(SE,altl,mol3,altitude,it,i)                         !
!----------------------------------------------------------------------------!
subroutine SE_alt_atb(var1,var2,var3,alt,nprofs,i)

  implicit none

  integer  ::  ilid, i, nprofs,ilid2,k,iSEdown,iSEup
  integer(kind=2)  ::  alt
  real  ::  echothres
  
  real*4,dimension(alt)  ::  var2
  real,dimension(nprofs)  ::  var1
  real*4,dimension(alt,nprofs)  ::  var3
              
iSEdown=0

! Control of SE value
! if SE=NaN ==> Forced SE to 0
if(var1(i).eq.-9999.)then
var1(i)=0
endif

! Put values below the Surface Elevation threshold to -9999.
! Find the lidar level called "iSEdown" matching with the Surface elevation
   b2 :do ilid=2,alt
         if ( (var1(i).ge.var2(ilid)) )then 
             iSEdown=ilid
              b1:  do ilid2=583,ilid,-1
                 var3(ilid2,i)=-9999.
                   enddo b1
                   exit b2
         endif
   enddo b2

echothres=1.

if(var1(i).gt.0.)then
  b4 :do ilid=2,iSEdown-33
         if(var3(ilid,i).gt.0.15)then
         echothres=0.4
         exit b4
         endif
      enddo b4
endif

iSEup=-9999

! Find the highest lidar ATB value of the lidar ground echo called "iSEup"
if(any(var3(iSEdown-30:iSEdown,i).ge.echothres))then
  b3: do ilid=iSEdown-30,iSEdown
     iSEup=ilid
        if(var3(ilid,i).eq.maxval(var3(iSEdown-30:iSEdown,i)))exit b3
      enddo  b3        
endif

if(iSEup.ne.-9999)then
! Put values from the ground to 3 levels above the iSEup level to -9999. 
   do ilid=iSEup-3,iSEdown
      var3(ilid,i)=-9999.
   enddo
endif



end subroutine SE_alt_atb
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
! *** SE_ALT_CHIM *** This subroutine detect the surface elevation and       !
!                     allocate -9999 to the variable below this threshold for!
!                     the CHIMERE grid                                       !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1   : Surface Elevation                                                 !
! var2   : altitude of model boxes                                           !
!----------------------------------------------------------------------------!
! var3   : varically averaged variable                                       !
!----------------------------------------------------------------------------!
! l      : loop index on Surface Elevation level                             !
! ilid   : loop index on the lidar altitude (583)                            !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call SE_alt_atb(SE,altl,mol3,altitude,it,i)                         !
!----------------------------------------------------------------------------!
subroutine SE_alt_mol(var1,var2,var3,alt,nprofs,i)

  implicit none

  integer  ::  ilid, i, nprofs,l
  integer(kind=2)  ::  alt
  real*4,dimension(alt)  ::  var2
  real,dimension(nprofs)  ::  var1
  real*4,dimension(alt,nprofs)  ::  var3
              
   b2 :do ilid=2,alt
         if ( (var1(i).gt.var2(ilid)) )then 
              b1:  do l=583,(ilid+1),-1
                 var3(l,i)=-9999
                   enddo b1
               exit b2
         endif 
   enddo b2

end subroutine SE_alt_mol
!----------------------------------------------------------------------------!


!****************************************************************************!



!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!! NETCDF RECORDING SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


!----------------------------------------------------------------------------!
! *** CHECK *** This subroutine check the status of the nf90 command         !
!----------------------------------------------------------------------------!
subroutine check(status)

    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if

end subroutine check  
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** CREATE_PROFNC *** This routine create a netcdf prof file, and its      !
!                       dimensions                                           !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf MAP3D file                              !
! dname      : period of MAP3D file (description of ncdf file)               !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! alt        : number of the altitude boxes                                  !
! vlon       : values of model longitude from 180 to -180                    !
! vlat       : values of model latitude from 90 to -90                       !
! vpres      : values of model altitude from 0 to 19.2km or 40.5km           !
! vtime      : number of days since 2000/01/01 for the trimonthly perdiod,   !
!              in day.                                                       !
! dim        : dimension id of the MaP3D variables recorded in the ncdf      !
!              files                                                         !
!----------------------------------------------------------------------------!
! date       : date from the real-system clock and has form yyyymmdd.        !
! time       : time from the real-system clock and has form hhmmss.sss.      !
! zone       : represente the difference with respect to Coordinated         !
!              Universal Time (UTC), and has form (+-)hhmm.                  !
! value      : 8 dimension value which contains the year,month,day,hour,     !
!              minute, seconds and milliseconds of the real-time.            !
! ndims      : dimension of dim (=4 : lon,lat,alt,time)                      !
! nc_id      : netcdf file id                                                !
! lon_varid  : variable id of longitude                                      !
! lat_varid  : variable id of latitude                                       !
! pres_varid : variable id of altitude                                       !
! time_varid : variable id of time                                           !
! lon_dimid  : dimension id of longitude                                     !
! lat_dimid  : dimension id of latitude                                      !
! pres_dimid : dimension id of altitude                                      !
! time_dimid : dimension id of time                                          !
!----------------------------------------------------------------------------!
!                                                                            !
! ex  create_profnc(file8,file9,lonmod,latmod,altmod,resd,dimidsp,altmax,    !
!                   lonmax,latmax)                                           !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine create_profnc(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vtime,dim,alt,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer, parameter ::  ndims = 4, nv=2
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat ,dim(ndims), alt,dim2(nv)
    integer  ::  lon_varid,lat_varid,alt_varid, time_varid, nc_id,alt_varid2
    integer  ::  lon_dimid,lat_dimid,alt_dimid,time_dimid, nv_dimid
    real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound

  
    call date_and_time(date,time,zone,value)
    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_Three-dimensionnal_Cloud_Fraction_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))

    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))


    dim = (/lon_dimid, lat_dimid, alt_dimid, time_dimid/)
    dim2 = (/alt_dimid, nv_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim2, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

 
   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  

    call check(nf90_enddef(nc_id))
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_profnc
!----------------------------------------------------------------------------!


subroutine create_temp3d(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,  &
                         vtime,dim,alt,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer, parameter ::  ndims = 4, nv=2
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat ,dim(ndims), alt,dim2(nv)
    integer  ::  lon_varid,lat_varid,alt_varid, time_varid, nc_id,alt_varid2
    integer  ::  lon_dimid,lat_dimid,alt_dimid,time_dimid, nv_dimid
    real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound

  
    call date_and_time(date,time,zone,value)
    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_Three-dimensionnal_Cloud_Fraction_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))

    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'temp', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))


    dim = (/lon_dimid, lat_dimid, alt_dimid, time_dimid/)
    dim2 = (/alt_dimid, nv_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'temp_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the temperature bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'temp_bound', NF90_FLOAT, dim2, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the temperature bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

 
   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  

    call check(nf90_enddef(nc_id))
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_temp3d
!----------------------------------------------------------------------------!


subroutine create_depolnc3d(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vtime,dim,dim2,alt,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer, parameter ::  ndims = 4, nv=2, ncat1=5,ncat2=25
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat ,dim(ndims), alt,dim3(nv),dim2(5),dim4(2)
    integer  ::  lon_varid,lat_varid,alt_varid, time_varid, nc_id,alt_varid2,cat_varid
    integer  ::  lon_dimid,lat_dimid,alt_dimid,time_dimid, nv_dimid,cat_dimid1,cat_dimid2
    real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound
    character(len=25),dimension(ncat1,ncat2)  ::  vcat
    

vcat(1,:)='UNDEFINED'
vcat(2,:)='FALSE LIQ'
vcat(3,:)='FALSE ICE'
vcat(4,:)='Horizontally Oriented'
vcat(5,:)='Unphysical value (NOISE)'

  
    call date_and_time(date,time,zone,value)
    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_Three-dimensionnal_CloudFraction_Phase_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))  
 
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'cat1', ncat1, cat_dimid1))
    call check(nf90_def_dim(nc_id, 'cat2', ncat2, cat_dimid2))

    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    dim = (/lon_dimid, lat_dimid, alt_dimid, time_dimid/)
    dim2 = (/lon_dimid, lat_dimid, alt_dimid, cat_dimid1, time_dimid/)
    dim3 = (/alt_dimid, nv_dimid/)
    dim4 = (/cat_dimid2, cat_dimid1/)
    
    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))
 
    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim3, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'category', NF90_CHAR, dim4, cat_varid))
    call check(nf90_put_att(nc_id, cat_varid, 'lon_name','Category'))
    call check(nf90_put_att(nc_id, cat_varid, 'units','Arbitrary unit'))

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  

    call check(nf90_enddef(nc_id))
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, cat_varid, vcat))
    call check(nf90_put_var(nc_id, time_varid, vtime))



    call check(nf90_close(nc_id))

end subroutine create_depolnc3d
!----------------------------------------------------------------------------!

subroutine create_ind3d(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vtime,dim,dim2,alt,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer, parameter ::  ndims = 4, nv=2, ncat1=5,ncat2=25
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat ,dim(ndims), alt,dim3(nv),dim2(5),dim4(2)
    integer  ::  lon_varid,lat_varid,alt_varid, time_varid, nc_id,alt_varid2,cat_varid
    integer  ::  lon_dimid,lat_dimid,alt_dimid,time_dimid, nv_dimid,cat_dimid1,cat_dimid2
    real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound
    character(len=25),dimension(ncat1,ncat2)  ::  vcat
    

vcat(1,:)='UNDEFINED'
vcat(2,:)='FALSE LIQ'
vcat(3,:)='FALSE ICE'
vcat(4,:)='Horizontally Oriented'
vcat(5,:)='Unphysical value (NOISE)'

  
    call date_and_time(date,time,zone,value)
    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_Three-dimensionnal_CloudFraction_Phase_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))  
 
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'cat1', ncat1, cat_dimid1))
    call check(nf90_def_dim(nc_id, 'cat2', ncat2, cat_dimid2))

    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    dim = (/lon_dimid, lat_dimid, alt_dimid, time_dimid/)
    dim2 = (/lon_dimid, lat_dimid, alt_dimid, cat_dimid1, time_dimid/)
    dim3 = (/alt_dimid, nv_dimid/)
    dim4 = (/cat_dimid2, cat_dimid1/)
    
    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))
 
    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim3, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'category', NF90_CHAR, dim4, cat_varid))
    call check(nf90_put_att(nc_id, cat_varid, 'lon_name','Category'))
    call check(nf90_put_att(nc_id, cat_varid, 'units','Arbitrary unit'))

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  

    call check(nf90_enddef(nc_id))
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, cat_varid, vcat))
    call check(nf90_put_var(nc_id, time_varid, vtime))



    call check(nf90_close(nc_id))

end subroutine create_ind3d
!----------------------------------------------------------------------------!




!----------------------------------------------------------------------------!
! *** INSTANTSR2NC *** This routine create & record a netcdf instant_SR file,!
!                      and its dimensions                                    !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf instant_SR file                         !
! daynight   : switch which allow to select day or night mode                !
! mod        : name of the grid selected                                     !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! alt        : number of the altitude boxes                                  !
! longi      : Longitude in degrees, has dimension (nprof)                   !
! lati       : Latitude in degrees, has dimension (nprof)                    !
! vprestop   : values of model altitude from 0 to 19.2km or 40.5km           !
! SEi        : Surface_Elevation in kilometers, has dimension (nprof)        !
! timei      : Time converted in UTC fractionned hour of a day               !
! SR         : scattering ratio calculated, dim=(altmax,nprof)               !
! vtime      : number of days since 2000/01/01 for the trimonthly perdiod,   !
!              in day.                                                       !
!----------------------------------------------------------------------------!
! fname2     : period of instant_SR file (description of ncdf file)          !
! date       : date from the real-system clock and has form yyyymmdd.        !
! time       : time from the real-system clock and has form hhmmss.sss.      !
! zone       : represente the difference with respect to Coordinated         !
!              Universal Time (UTC), and has form (+-)hhmm.                  !
! value      : 8 dimension value which contains the year,month,day,hour,     !
!              minute, seconds and milliseconds of the real-time.            !
! dim        : dimension id of the SR                                        !
! ndims      : dimension of dim (= 2 : alt,nprof)                            !
! nc_id      : netcdf file id                                                !
! varid3     : variable id of longitude                                      !
! varid2     : variable id of latitude                                       !
! alt_varid  : variable id of altitude                                       !
! varid5     : variable id of time                                           !
! varid6     : variable id of SE                                             !
! varid1     : variable id of instant_SR                                     !
! alt_dimid  : dimension id of altitude                                      !
! it_dimid   : dimension id of nprof                                         !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : instantSR2nc(file4,altmod,resd,altmax,switch,gcm,it,lat,lon,SE,temps2,!
!                   srmoy)                                                   !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine instantSR2nc(fname,vprestop_mid,vprestop_bound,vtime,alt,daynight,mod,nprof,lati,    &
                        longi,SEi,timei,SR)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname
    character  ::  fname2*14
    character  ::  date*8,time*10,zone*5,daynight*5,mod*8
    integer, parameter ::  ndims = 2, nv=2
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims), alt, nprof,dim2(ndims)
    integer  ::  alt_varid, varid1, varid2, varid3, varid5, varid6, nc_id,alt_varid2
    integer  ::  alt_dimid, it_dimid, nv_dimid, cat_varid
    real  ::  vtime
    real*8, dimension(nprof) :: timei  
    real, dimension(nprof) :: SEi    
    real, dimension(nprof) :: longi
    real, dimension(nprof) :: lati
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound

    real*4,dimension(alt,nprof)  ::  SR
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.

fname2=fname(12:25)

    call date_and_time(date,time,zone,value)
    call check(nf90_create('/bdd/CFMIP/GOCCP/instant_SR/temp/'//fname,     &
               NF90_CLOBBER, nc_id))

    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_instant_SR_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',fname2))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 


    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'it', nprof, it_dimid))

    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))


    dim = (/alt_dimid, it_dimid/)
    dim2 = (/alt_dimid, nv_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, it_dimid, varid3))
    call check(nf90_put_att(nc_id, varid3, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, varid3, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, varid3, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, it_dimid, varid2))
    call check(nf90_put_att(nc_id, varid2, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, varid2, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, varid2, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim2, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, it_dimid, varid5))
    call check(nf90_put_att(nc_id, varid5, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, varid5, 'units','Fractionned UTC hour'))
    call check(nf90_put_att(nc_id, varid5, 'axis','T')) 

    call check(nf90_def_var(nc_id, 'SE', NF90_FLOAT, it_dimid, varid6))
    call check(nf90_put_att(nc_id, varid6, 'lon_name','Surface_Elevation'))
    call check(nf90_put_att(nc_id, varid6, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'instant_SR', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(nc_id, varid1, 'lon_name',                       &
               'instantaneous Scattering Ratio'))
    call check(nf90_put_att(nc_id, varid1, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid1, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid1, 'rejected_value',rej))    
    call check(nf90_put_att(nc_id, varid1, 'Surface_value',se))

    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, varid2, lati))
    call check(nf90_put_var(nc_id, varid3, longi))
    call check(nf90_put_var(nc_id, varid5, timei))
    call check(nf90_put_var(nc_id, varid6, SEi))

    call check(nf90_put_var(nc_id, varid1, SR))
    call check(nf90_close(nc_id))

end subroutine instantSR2nc
!----------------------------------------------------------------------------!

subroutine SR_CR_DR_2nc(fname,vprestop_mid,vprestop_bound,vtime,alt,daynight,mod,nprof,lati,    &
                        longi,SEi,timei,SR,CR,DEPOL)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname
    character  ::  fname2*14
    character  ::  date*8,time*10,zone*5,daynight*5,mod*4
    integer, parameter ::  ndims = 2, nv = 2 
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims), alt, nprof, dim2(ndims)
    integer  ::  alt_varid, alt_varid2, varid1, varid2, varid3, varid5, varid6,varid7,varid8, nc_id
    integer  ::  alt_dimid, it_dimid, nv_dimid
    real  ::  vtime
    real*8, dimension(nprof) :: timei  
    real, dimension(nprof) :: SEi    
    real, dimension(nprof) :: longi
    real, dimension(nprof) :: lati
    real*4,dimension(alt,nprof)  ::  SR,CR,DEPOL
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound

fname2=fname(18:39)

    call date_and_time(date,time,zone,value)
    call check(nf90_create('/bdd/CFMIP/GOCCP/instant_SR_CR_DR/temp/'//fname,     &
               NF90_CLOBBER, nc_id))

    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_instant_SR_DR_CR_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',trim(fname2)))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 

    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'it', nprof, it_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    dim = (/alt_dimid, it_dimid/)
    dim2 = (/alt_dimid, nv_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, it_dimid, varid3))
    call check(nf90_put_att(nc_id, varid3, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, varid3, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, varid3, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, it_dimid, varid2))
    call check(nf90_put_att(nc_id, varid2, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, varid2, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, varid2, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim2, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer')) 

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, it_dimid, varid5))
    call check(nf90_put_att(nc_id, varid5, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, varid5, 'units','Fractionned UTC hour'))
    call check(nf90_put_att(nc_id, varid5, 'axis','T')) 

    call check(nf90_def_var(nc_id, 'SE', NF90_FLOAT, it_dimid, varid6))
    call check(nf90_put_att(nc_id, varid6, 'lon_name','Surface_Elevation'))
    call check(nf90_put_att(nc_id, varid6, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'instant_SR', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(nc_id, varid1, 'lon_name',                       &
               'instantaneous Scattering Ratio'))
    call check(nf90_put_att(nc_id, varid1, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid1, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid1, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid1, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'instant_CR', NF90_FLOAT, dim, varid7))
    call check(nf90_put_att(nc_id, varid7, 'lon_name',                       &
               'instantaneous Color Ratio'))
    call check(nf90_put_att(nc_id, varid7, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid7, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid7, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid7, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'instant_DR', NF90_FLOAT, dim, varid8))
    call check(nf90_put_att(nc_id, varid8, 'lon_name',                       &
               'instantaneous Depolarization Ratio'))
    call check(nf90_put_att(nc_id, varid8, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid8, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid8, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid8, 'Surface_value',se))


    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, varid2, lati))
    call check(nf90_put_var(nc_id, varid3, longi))
    call check(nf90_put_var(nc_id, varid5, timei))
    call check(nf90_put_var(nc_id, varid6, SEi))

    call check(nf90_put_var(nc_id, varid1, SR))
    call check(nf90_put_var(nc_id, varid7, CR))
    call check(nf90_put_var(nc_id, varid8, DEPOL))



    call check(nf90_close(nc_id))

end subroutine SR_CR_DR_2nc
!----------------------------------------------------------------------------!

! Same routine as SR_CR_DR_2nc including the record of ATB & ATBmol variables
subroutine SR_CR_DR_ATB_nc(fname,vprestop_mid,vprestop_bound,vtime,alt,      &
                           daynight,mod,nprof,lati,longi,SEi,timei,SR,CR,    &
                           DEPOL,atb,atbm,atbr,atbl,temp,cloudfrac)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname
    character  ::  fname2*21
    character  ::  date*8,time*10,zone*5,daynight*5,mod*4
    integer, parameter ::  ndims = 2, nv = 2 
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims), alt, nprof, dim2(ndims)
    integer  ::  alt_varid, alt_varid2, varid1, varid2, varid3, varid5, varid6,varid7,varid8,varid9,varid10, nc_id,varid11,varid12,varid13, varid14
    integer  ::  alt_dimid, it_dimid, nv_dimid
    real  ::  vtime
    real*8, dimension(nprof) :: timei  
    real, dimension(nprof) :: SEi    
    real, dimension(nprof) :: longi
    real, dimension(nprof) :: lati
    real,dimension(alt,nprof)  ::  SR,CR,DEPOL,atb,atbm,atbl,atbr,temp,cloudfrac
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound

fname2=fname(18:39)

    call date_and_time(date,time,zone,value)
    call check(nf90_create('/bdd/CFMIP/GOCCP/instant_SR_CR_DR/temp/'//fname,     &
               NF90_CLOBBER, nc_id))

    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_instant_SR_DR_CR_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',trim(fname2)))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 

    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'it', nprof, it_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    dim = (/alt_dimid, it_dimid/)
    dim2 = (/alt_dimid, nv_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, it_dimid, varid3))
    call check(nf90_put_att(nc_id, varid3, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, varid3, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, varid3, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, it_dimid, varid2))
    call check(nf90_put_att(nc_id, varid2, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, varid2, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, varid2, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim2, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer')) 

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, it_dimid, varid5))
    call check(nf90_put_att(nc_id, varid5, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, varid5, 'units','Fractionned UTC hour'))
    call check(nf90_put_att(nc_id, varid5, 'axis','T')) 

    call check(nf90_def_var(nc_id, 'SE', NF90_FLOAT, it_dimid, varid6))
    call check(nf90_put_att(nc_id, varid6, 'lon_name','Surface_Elevation'))
    call check(nf90_put_att(nc_id, varid6, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'instant_SR', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(nc_id, varid1, 'lon_name',                       &
               'instantaneous Scattering Ratio'))
    call check(nf90_put_att(nc_id, varid1, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid1, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid1, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid1, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'instant_Cloud', NF90_FLOAT, dim, varid14))
    call check(nf90_put_att(nc_id, varid14, 'lon_name',                       &
               'instantaneous Cloud Mask'))
    call check(nf90_put_att(nc_id, varid14, 'Legend','1=NaN, 2=Clear, 3=Cloud, 4=Uncertain, 6=Surface, 7=Rejected, 8=Fully Attenuated' ))
    call check(nf90_put_att(nc_id, varid14, 'units','Arbitrary unit'))

    call check(nf90_def_var(nc_id, 'instant_CR', NF90_FLOAT, dim, varid7))
    call check(nf90_put_att(nc_id, varid7, 'lon_name',                       &
               'instantaneous Color Ratio'))
    call check(nf90_put_att(nc_id, varid7, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid7, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid7, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid7, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'instant_DR', NF90_FLOAT, dim, varid8))
    call check(nf90_put_att(nc_id, varid8, 'lon_name',                       &
               'instantaneous Depolarization Ratio'))
    call check(nf90_put_att(nc_id, varid8, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid8, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid8, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid8, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'ATB', NF90_FLOAT, dim, varid9))
    call check(nf90_put_att(nc_id, varid9, 'lon_name',                       &
               'Attenuated Total Backscatter'))
    call check(nf90_put_att(nc_id, varid9, 'units','km-1 sr-1'))
    call check(nf90_put_att(nc_id, varid9, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid9, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid9, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'ATB_mol', NF90_FLOAT, dim, varid10))
    call check(nf90_put_att(nc_id, varid10, 'lon_name',                       &
               'ATB of Molecular'))
    call check(nf90_put_att(nc_id, varid10, 'units','km-1 sr-1'))
    call check(nf90_put_att(nc_id, varid10, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid10, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid10, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'ATB_per', NF90_FLOAT, dim, varid11))
    call check(nf90_put_att(nc_id, varid11, 'lon_name',                       &
               'Attenuated Perpendicular Backscatter'))
    call check(nf90_put_att(nc_id, varid11, 'units','km-1 sr-1'))
    call check(nf90_put_att(nc_id, varid11, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid11, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid11, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'ATB_par', NF90_FLOAT, dim, varid12))
    call check(nf90_put_att(nc_id, varid12, 'lon_name',                       &
               'Attenuated Parallel Backscatter'))
    call check(nf90_put_att(nc_id, varid12, 'units','km-1 sr-1'))
    call check(nf90_put_att(nc_id, varid12, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid12, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid12, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'TEMP', NF90_FLOAT, dim, varid13))
    call check(nf90_put_att(nc_id, varid13, 'lon_name',                       &
               'Temperature'))
    call check(nf90_put_att(nc_id, varid13, 'units','Celcius degree'))
    call check(nf90_put_att(nc_id, varid13, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid13, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid13, 'Surface_value',se))


    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, varid2, lati))
    call check(nf90_put_var(nc_id, varid3, longi))
    call check(nf90_put_var(nc_id, varid5, timei))
    call check(nf90_put_var(nc_id, varid6, SEi))

    call check(nf90_put_var(nc_id, varid1, SR))
    call check(nf90_put_var(nc_id, varid14, cloudfrac))
    call check(nf90_put_var(nc_id, varid7, CR))
    call check(nf90_put_var(nc_id, varid8, DEPOL))
    call check(nf90_put_var(nc_id, varid9, atb))
    call check(nf90_put_var(nc_id, varid10, atbm))
    call check(nf90_put_var(nc_id, varid11, atbr))
    call check(nf90_put_var(nc_id, varid12, atbl))
    call check(nf90_put_var(nc_id, varid13, temp))




    call check(nf90_close(nc_id))

end subroutine SR_CR_DR_ATB_nc
!----------------------------------------------------------------------------!


subroutine instant_phase(fname,nalt,nprof,phase)
implicit none
character(len=*)  ::  fname
integer,dimension(2)  ::  dimm
integer  ::  ncid,varid, alt_dimid,it_dimid,ndims
integer  ::  nalt,nprof
real*4  ::  phase(nalt,nprof)

!print *, 'instant phase routine'

!print *, trim(fname)
!print *, '/bdd/CFMIP/GOCCP/instant_SR_CR_DR/temp/'//trim(fname)

    call check(NF90_OPEN('/bdd/CFMIP/GOCCP/instant_SR_CR_DR/temp/'//trim(fname),NF90_WRITE,ncid))

call check(NF90_REDEF(ncid))

    call check(nf90_inq_dimid(ncid, 'altitude',  alt_dimid))
    call check(nf90_inq_dimid(ncid, 'it', it_dimid))

    dimm = (/alt_dimid, it_dimid/)

    call check(nf90_def_var(ncid, 'instant_Phase', NF90_FLOAT, dimm, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                       &
               'instantaneous Cloud Phase Mask'))
    call check(nf90_put_att(ncid, varid, 'Legend','1=LIQ, 2=ICE, 3=UNDEFINED, 4=FALSE LIQ, 5=FALSE ICE, 6=Horizontally Oriented, 7=Unphysical value (NOISE)' ))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))
    call check(nf90_put_att(ncid, varid, 'Surface_value','-888'))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid, phase))

    call check(nf90_close(ncid))

endsubroutine instant_phase


subroutine SR_DEPOL_2nc(fname,vprestop,vtime,alt,daynight,mod,nprof,lati,    &
                        longi,SEi,timei,SR,DEPOL)
!!$  call SR_DEPOL_2nc(file4,altmod,resd,altmax,switch,gcm,it,lat,lon,SE,temps2,&
!!$                  srmoy,depolmoy)
    use netcdf
    implicit none

    character(LEN=*)   ::  fname
    character  ::  fname2*13
    character  ::  date*8,time*10,zone*5,daynight*5,mod*8
    integer, parameter ::  ndims = 2
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims), alt, nprof
    integer  ::  alt_varid, varid1, varid2, varid3, varid5, varid6,varid7,varid8, nc_id
    integer  ::  alt_dimid, it_dimid
    real  ::  vtime
    real*8, dimension(nprof) :: timei  
    real, dimension(nprof) :: SEi    
    real, dimension(nprof) :: longi
    real, dimension(nprof) :: lati
    real,dimension(alt)  ::  vprestop
    real*4,dimension(alt,nprof)  ::  SR,DEPOL

fname2=fname(12:26)

    call date_and_time(date,time,zone,value)
    call check(nf90_create('/bdd/CFMIP_TEMP/DEPOL/'//fname,     &
               NF90_CLOBBER, nc_id))

    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_instant_SR_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',fname2))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))

    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'it', nprof, it_dimid))

    dim = (/alt_dimid, it_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, it_dimid, varid3))
    call check(nf90_put_att(nc_id, varid3, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, varid3, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, varid3, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, it_dimid, varid2))
    call check(nf90_put_att(nc_id, varid2, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, varid2, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, varid2, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'altitude', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Altitude'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, it_dimid, varid5))
    call check(nf90_put_att(nc_id, varid5, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, varid5, 'units','Fractionned UTC hour'))
    call check(nf90_put_att(nc_id, varid5, 'axis','T')) 

    call check(nf90_def_var(nc_id, 'SE', NF90_FLOAT, it_dimid, varid6))
    call check(nf90_put_att(nc_id, varid6, 'lon_name','Surface_Elevation'))
    call check(nf90_put_att(nc_id, varid6, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'instant_SR', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(nc_id, varid1, 'lon_name',                       &
               'instantaneous Scattering Ratio'))
    call check(nf90_put_att(nc_id, varid1, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid1, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid1, 'Surface_value','-888'))

    call check(nf90_def_var(nc_id, 'instant_DEPOL', NF90_FLOAT, dim, varid8))
    call check(nf90_put_att(nc_id, varid8, 'lon_name',                       &
               'instantaneous Depolarization'))
    call check(nf90_put_att(nc_id, varid8, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid8, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid8, 'Surface_value','-888'))

    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, alt_varid, vprestop))
    call check(nf90_put_var(nc_id, varid2, lati))
    call check(nf90_put_var(nc_id, varid3, longi))
    call check(nf90_put_var(nc_id, varid5, timei))
    call check(nf90_put_var(nc_id, varid6, SEi))

    call check(nf90_put_var(nc_id, varid1, SR))
    call check(nf90_put_var(nc_id, varid8, DEPOL))



    call check(nf90_close(nc_id))

end subroutine SR_DEPOL_2nc
!----------------------------------------------------------------------------!












!----------------------------------------------------------------------------!
! *** PROF_RECVAR2NC *** This routine record the prof variables in the netcdf!
!                        prof file                                           !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf MAP3D file                              !
! dim        : dimension id of the SR                                        !
! cloud      : fraction of cloudy point                                      !
! clear      : fraction of clear point                                       !
! sat        : fraction of fully attenuated point                            !
! uncer      : fraction of unclassify point                                  !
! nan        : fraction of nan point                                         !
! se         : fraction of se point                                          !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! alt        : number of the altitude boxes                                  !
!----------------------------------------------------------------------------!
! nanb       : Not a Number value                                            !
! ndims      : dimension of dim (= 4 : lon,lat,alt,time)                     !
! ncid       : netcdf file id                                                !
! varid3     : variable id of cloud fraction                                 !
! varid4     : variable id of clear fraction                                 !
! varid5     : variable id of fully attenuated fraction                      !
! varid6     : variable id of unclassify fraction                            !
! varid7     : variable id of nan fraction                                   !
! varid8     : variable id of se fraction                                    !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : prof_recvar2nc(monthcloudfract,monthclearfract,monthsatfract,         !
!                     monthuncertfract,monthnanfract,monthsefract,dimidsp,   ! 
!                     file8,altmax,lonmax,latmax)                            !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine prof_recvar2nc(cloud,clear,uncer,dim,fname,alt,nlon,nlat)!nan,se,sat
    use netcdf
    implicit none

    integer, parameter ::  ndims=4
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer  ::  nlon , nlat,dim(ndims),alt
    integer  ::  varid3,varid4,varid5,varid6,varid7,varid8, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt)  ::  cloud,clear,uncer!,sat,nan,se

    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))
    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'clcalipso', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO 3D Cloud fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clrcalipso', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO 3D Clear fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

   ! call check(nf90_def_var(ncid, 'monthsatfract', NF90_FLOAT, dim, varid5))
   ! call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
   !            'Full Attenuated fraction monthly mean'))
   ! call check(nf90_put_att(ncid, varid5, 'units','1 fraction'))
   ! call check(nf90_put_att(ncid, varid5, '_FillValue',nanb))

   call check(nf90_def_var(ncid, 'uncalipso', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO 3D Undefined fraction'))
    call check(nf90_put_att(ncid, varid6, 'units','1 fraction'))
   call check(nf90_put_att(ncid, varid6, '_FillValue',nan))

!    call check(nf90_def_var(ncid, 'monthnanfract', NF90_FLOAT, dim, varid7))
!    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
!               'Missing value fraction monthly mean'))
!    call check(nf90_put_att(ncid, varid7, 'units','1 fraction'))
 !   call check(nf90_put_att(ncid, varid7, '_FillValue',nanb))

 !   call check(nf90_def_var(ncid, 'monthsefract', NF90_FLOAT, dim, varid8))
 !   call check(nf90_put_att(ncid, varid8, 'lon_name',                        &
 !              'Surface Elevation fraction monthly mean'))
 !   call check(nf90_put_att(ncid, varid8, 'units','1 fraction'))
 !   call check(nf90_put_att(ncid, varid8, '_FillValue',nanb))

!!$    call check(nf90_def_var(ncid, 'monthindphase', NF90_FLOAT, dim, varid9))
!!$    call check(nf90_put_att(ncid, varid9, 'lon_name','indice of water Phase monthly mean'))
!!$    call check(nf90_put_att(ncid, varid9, 'units','1 fraction'))
!!$    call check(nf90_put_att(ncid, varid9, '_FillValue',nanb))
!!$  
!!$    call check(nf90_def_var(ncid, 'monthatbmoy', NF90_FLOAT, dim, varid1))
!!$    call check(nf90_put_att(ncid, varid1, 'lon_name','Total Attenuated Backscatter 532 monthly mean'))
!!$    call check(nf90_put_att(ncid, varid1, 'units','per kilometer per steradian'))
!!$    call check(nf90_put_att(ncid, varid1, '_FillValue',nanb))
!!$
!!$    call check(nf90_def_var(ncid, 'monthatbmolmoy', NF90_FLOAT, dim, varid2))
!!$    call check(nf90_put_att(ncid, varid2, 'lon_name','Molecular Total Attenuated Backscatter monthly mean'))
!!$    call check(nf90_put_att(ncid, varid2, 'units','per kilometer per steradian'))
!!$    call check(nf90_put_att(ncid, varid2, '_FillValue',nanb))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid3, cloud))
    call check(nf90_put_var(ncid, varid4, clear))
  !  call check(nf90_put_var(ncid, varid5, sat))
    call check(nf90_put_var(ncid, varid6, uncer))
   ! call check(nf90_put_var(ncid, varid7, nan))
   ! call check(nf90_put_var(ncid, varid8, se))
!!$    call check(nf90_put_var(ncid, varid1, pr2moy))
!!$    call check(nf90_put_var(ncid, varid2, molmoy))
!!$    call check(nf90_put_var(ncid, varid9, phase))

    call check(nf90_close(ncid))

endsubroutine prof_recvar2nc
!----------------------------------------------------------------------------!

subroutine temp_recvar2nc(cloud,liq,ice,phase,dim,fname,alt,nlon,nlat)!nan,se,sat
    use netcdf
    implicit none

    integer, parameter ::  ndims=4
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer  ::  nlon , nlat,dim(ndims),alt
    integer  ::  varid3,varid4,varid5,varid6,varid7,varid8, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt)  ::  cloud,ice,liq,phase !,sat,nan,se

    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))
    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cltemp', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO 3D Cloud fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltemp_liq', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO 3D Liquid fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltemp_ice', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO 3D Ice fraction'))
    call check(nf90_put_att(ncid, varid6, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltemp_phase', NF90_FLOAT, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'CALIPSO 3D Ratio fraction'))
    call check(nf90_put_att(ncid, varid7, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid3, cloud))
    call check(nf90_put_var(ncid, varid4, liq))
    call check(nf90_put_var(ncid, varid6, ice))
    call check(nf90_put_var(ncid, varid7, phase))

    call check(nf90_close(ncid))

endsubroutine temp_recvar2nc
!----------------------------------------------------------------------------!


subroutine depol_recvar2nc(ice,water,dim,fname,alt,nlon,nlat)!nan,se,sat
    use netcdf
    implicit none

    integer, parameter ::  ndims=4
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer  ::  nlon , nlat,dim(ndims),alt
    integer  ::  varid3,varid4,varid5,varid6,varid7,varid8, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt)  ::  ice,water!,ind !,sat,nan,se

    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))
    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'ice_cloud', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO 3D Ice Cloud Phase fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'water_cloud', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO 3D Water Cloud Phase fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

!!$    call check(nf90_def_var(ncid, 'ind_phase', NF90_FLOAT, dim, varid5))
!!$    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
!!$               'CALIPSO 3D Indice Phase fraction'))
!!$    call check(nf90_put_att(ncid, varid5, 'units','1 fraction'))
!!$    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))
 
    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid3, ice))
    call check(nf90_put_var(ncid, varid4, water))
!    call check(nf90_put_var(ncid, varid5, ind))

    call check(nf90_close(ncid))

endsubroutine depol_recvar2nc
!----------------------------------------------------------------------------!

subroutine depol_recvar2ncocc(ice,water,un,phase,ind,dim,dim2,fname,alt,nlon,nlat)!nan,se,sat
    use netcdf
    implicit none

    integer, parameter ::  ndims=4, ncat=5
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer  ::  nlon , nlat,dim(ndims),alt,dim2(ncat)
    integer  ::  varid3,varid4,varid5,varid6,varid7,varid8, varid9, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt)  ::  ice,water,phase,ind
    real*4,dimension(nlon,nlat,alt,ncat)  ::  un  

    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))
    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'clcalipso_ice', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO 3D Ice Cloud Phase fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_liq', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO 3D Liquid Cloud Phase fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_un', NF90_FLOAT, dim2, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO 3D UNCLASS Cloud Phase fraction'))
    call check(nf90_put_att(ncid, varid5, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_RPIC', NF90_FLOAT, dim, varid9))
    call check(nf90_put_att(ncid, varid9, 'lon_name',                        &
               'CALIPSO 3D Relative Percentage of Ice in Cloud'))
    call check(nf90_put_att(ncid, varid9, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid9, '_FillValue',nan))
    call check(nf90_put_att(ncid, varid9, 'Other_Phase',rej))

!!$    call check(nf90_def_var(ncid, 'ind', NF90_FLOAT, dim, varid8))
!!$    call check(nf90_put_att(ncid, varid8, 'lon_name',                        &
!!$               'CALIPSO 3D Indice Phase fraction'))
!!$    call check(nf90_put_att(ncid, varid8, 'units','1 fraction'))
!!$    call check(nf90_put_att(ncid, varid8, '_FillValue',nan))
!!$ 
    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid3, ice))
    call check(nf90_put_var(ncid, varid4, water))
    call check(nf90_put_var(ncid, varid5, un))
    call check(nf90_put_var(ncid, varid9, phase))
!!$     call check(nf90_put_var(ncid, varid8, ind))

    call check(nf90_close(ncid))

endsubroutine depol_recvar2ncocc
!----------------------------------------------------------------------------!
! call record_ind3d(sum(cloudfractday,4),sum(indday,4),sum(icecloudfractday,4), &
!                  sum(watercloudfractday,4),sum(uncloudfractday,4),  &
!                  dimidsp,dimidsp2,file8,altmax,lonmax-1,latmax-1)
subroutine record_ind3d(cloud,tot,ice,water,un,dim,dim2,fname,alt,nlon,nlat)!nan,se,sat
    use netcdf
    implicit none

    integer, parameter ::  ndims=4, ncat=5
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer  ::  nlon , nlat,dim(ndims),alt,dim2(ncat)
    integer  ::  varid2,varid1,varid3,varid4,varid5,varid6,varid7,varid8, varid9, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt)  ::  tot,cloud,ice,water,phase
    real*4,dimension(nlon,nlat,alt,ncat)  ::  un  

    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))
    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'clcalipso', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO 3D Cloud'))
    call check(nf90_put_att(ncid, varid2, 'units','1'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'sample', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO 3D Sample'))
    call check(nf90_put_att(ncid, varid1, 'units','1'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_ice', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO 3D Ice Cloud Phase occurrences'))
    call check(nf90_put_att(ncid, varid3, 'units','1'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_liq', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO 3D Liquid Cloud Phase occurrences'))
    call check(nf90_put_att(ncid, varid4, 'units','1 occurrences'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_un', NF90_FLOAT, dim2, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO 3D UNCLASS Cloud Phase occurrences'))
    call check(nf90_put_att(ncid, varid5, 'units','1'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))


!!$    call check(nf90_def_var(ncid, 'ind', NF90_FLOAT, dim, varid8))
!!$    call check(nf90_put_att(ncid, varid8, 'lon_name',                        &
!!$               'CALIPSO 3D Indice Phase occurrences'))
!!$    call check(nf90_put_att(ncid, varid8, 'units','1 fraction'))
!!$    call check(nf90_put_att(ncid, varid8, '_FillValue',nan))
!!$ 
    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid2, cloud))
    call check(nf90_put_var(ncid, varid1, tot))
    call check(nf90_put_var(ncid, varid3, ice))
    call check(nf90_put_var(ncid, varid4, water))
    call check(nf90_put_var(ncid, varid5, un))

    call check(nf90_close(ncid))

endsubroutine record_ind3d
!----------------------------------------------------------------------------!




!----------------------------------------------------------------------------!
! *** CREATE_MAPNC ***  This routine create a netcdf map file, and its       !
!                       dimensions                                           !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf MapLowMidHigh file                      !
! dname      : period of MapLowMidHigh file (description of ncdf file)       !
! grid       : name of the grid selected (cfmip/lmdz/nasa)                   !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! vlon       : values of model longitude from 180 to -180                    !
! vlat       : values of model latitude from 90 to -90                       !
! vtime      : number of days since 2000/01/01 for the trimonthly perdiod,   !
!              in day.                                                       !
! dim        : dimension id of the MaPLowMidHigh var recorded in the ncdf    !
!----------------------------------------------------------------------------!
! toplvl     : values of isccp top levels                                    !
!              files                                                         !
! ndims      : dimension of dim (=3 : lon,lat,time)                          !
! nc_id      : netcdf file id                                                !
! lon_varid  : variable id of longitude                                      !
! lat_varid  : variable id of latitude                                       !
! pres_varid : variable id of toplvl                                         !
! time_varid : variable id of time                                           !
! lon_dimid  : dimension id of longitude                                     !
! lat_dimid  : dimension id of latitude                                      !
! pres_dimid : dimension id of toplvl                                        !
! time_dimid : dimension id of time                                          !
! date       : date from the real-system clock and has form yyyymmdd.        !
! time       : time from the real-system clock and has form hhmmss.sss.      !
! zone       : represente the difference with respect to Coordinated         !
!              Universal Time (UTC), and has form (+-)hhmm.                  !
! value      : 8 dimension value which contains the year,month,day,hour,     !
!              minute, seconds and milliseconds of the real-time.            !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : create_mapnc(file8,file9,lonmod,latmod,resd,dimidsm,gcm,lonmax,latmax)!
!                                                                            !
!----------------------------------------------------------------------------!
subroutine create_mapnc(fname,dname,vlon,vlat,vtime,dim,grid,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5,grid*8
    integer, parameter ::  ndims = 3
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat,dim(ndims)
    integer  ::  lon_varid,lat_varid,pres_varid, time_varid,nc_id
    integer  ::  lon_dimid,lat_dimid,pres_dimid,time_dimid
	real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(3)  ::  toplvl

if(trim(grid).eq.'CFMIP')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2


elseif(trim(grid).eq.'NASA')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2

elseif(trim(grid).eq.'LMDZ40')then
toplowl = 14;       
topmidl = 18;      
tophighl = 40  
elseif(trim(grid).eq.'LMDZ')then
toplowl = 7;        
topmidl = 9       
tophighl = 19   !
endif
  
    call date_and_time(date,time,zone,value)

    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Map_Low_Mid_High_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'toplvl', 3, pres_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
 
     dim = (/lon_dimid, lat_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

   call check(nf90_def_var(nc_id, 'toplvl', NF90_FLOAT,pres_dimid,pres_varid))
   call check(nf90_put_att(nc_id, pres_varid,'lon_name','Top Altitude Level'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
	
    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, pres_varid, toplvl))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_mapnc
!----------------------------------------------------------------------------!


subroutine create_mapnc_phase(fname,dname,vlon,vlat,vtime,dim,dim2,grid,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5,grid*8
    integer, parameter ::  ndims = 3, ndims2 = 4, ncat1=5, ncat2=25
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat,dim(ndims), dim2(ndims2),dim3(2)
    integer  ::  lon_varid,lat_varid,pres_varid, time_varid,nc_id,cat_varid
    integer  ::  lon_dimid,lat_dimid,pres_dimid,time_dimid, cat_dimid1,cat_dimid2
	real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(3)  ::  toplvl
    character(len=25),dimension(ncat1,ncat2)  ::  vcat
!    character(len=10),dimension(ncat) :: vcat

vcat(1,:)='UNDEFINED'
vcat(2,:)='FALSE LIQ'
vcat(3,:)='FALSE ICE'
vcat(4,:)='Horizontally Oriented'
vcat(5,:)='Unphysical Value (NOISE)'

if(trim(grid).eq.'CFMIP')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2


elseif(trim(grid).eq.'NASA')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2

elseif(trim(grid).eq.'LMDZ40')then
toplowl = 14;       
topmidl = 18;      
tophighl = 40  
elseif(trim(grid).eq.'LMDZ')then
toplowl = 7;        
topmidl = 9       
tophighl = 19   !
endif
  
    call date_and_time(date,time,zone,value)

    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Map_Low_Mid_High_Phase_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'cat1', ncat1, cat_dimid1))
    call check(nf90_def_dim(nc_id, 'cat2', ncat2, cat_dimid2))

    call check(nf90_def_dim(nc_id, 'toplvl', 3, pres_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
 
     dim = (/lon_dimid, lat_dimid, time_dimid/)
     dim2 = (/lon_dimid, lat_dimid, cat_dimid1, time_dimid/)
     dim3=(/cat_dimid2,cat_dimid1 /)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'category', NF90_CHAR, dim3, cat_varid))
    call check(nf90_put_att(nc_id, cat_varid, 'lon_name','Category'))
    call check(nf90_put_att(nc_id, cat_varid, 'units','Arbitrary unit'))

    call check(nf90_def_var(nc_id, 'toplvl', NF90_FLOAT,pres_dimid,pres_varid))
    call check(nf90_put_att(nc_id, pres_varid,'lon_name','Top Altitude Level'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
	
    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, cat_varid, vcat))
    call check(nf90_put_var(nc_id, pres_varid, toplvl))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_mapnc_phase
!----------------------------------------------------------------------------!



subroutine create_mapnc2(fname,dname,vlon,vlat,vtime,dim,dim2,dim3,dim4,grid,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5,grid*5
    integer, parameter ::  ndims = 3, ndims2=4, histmax=10, histmax2=40,histmax3=28
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat,dim(ndims), dim2(ndims2),dim3(ndims2),dim4(ndims2)
    integer  ::  lon_varid,lat_varid,pres_varid, time_varid,nc_id
    integer  ::  lon_dimid,lat_dimid,pres_dimid,time_dimid,hist_dimid,hist2_dimid,hist3_dimid
	real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(3)  ::  toplvl

if(trim(grid).eq.'CFMIP')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2


elseif(trim(grid).eq.'NASA')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2

elseif(trim(grid).eq.'LMDZ40')then
toplowl = 14;       
topmidl = 18;      
tophighl = 40  
elseif(trim(grid).eq.'LMDZ')then
toplowl = 7;        
topmidl = 9       
tophighl = 19   !
endif
  
    call date_and_time(date,time,zone,value)

    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Map_Low_Mid_High_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'hist', histmax, hist_dimid))
    call check(nf90_def_dim(nc_id, 'hist2', histmax2, hist2_dimid))
    call check(nf90_def_dim(nc_id, 'hist3', histmax3, hist3_dimid))

    call check(nf90_def_dim(nc_id, 'toplvl', 3, pres_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
 
     dim = (/lon_dimid, lat_dimid, time_dimid/)
     dim2 = (/lon_dimid, lat_dimid, hist_dimid, time_dimid/)
     dim3 = (/lon_dimid, lat_dimid, hist2_dimid, time_dimid/)
     dim4 = (/lon_dimid, lat_dimid, hist3_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

   call check(nf90_def_var(nc_id, 'toplvl', NF90_FLOAT,pres_dimid,pres_varid))
   call check(nf90_put_att(nc_id, pres_varid,'lon_name','Top Altitude Level'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
	
    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, pres_varid, toplvl))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_mapnc2
!----------------------------------------------------------------------------!




!----------------------------------------------------------------------------!
! *** MAP_RECVAR2NC *** This routine record the map variables in the netcdf  !
!                       prof file                                            !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf MaPLowMidHigh file                      !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! dim        : dimension id of the isccp variables                           !
! low        : isccp low cloud flag                                          !
! mid        : isccp mid cloud flag                                          !
! high       : isccp high cloud flag                                         !
! colcloud   : isccp column cloud flag                                       !
! colclear   : isccp column clear flag                                       !
!----------------------------------------------------------------------------!
! nan        : Not a Number value                                            !
! ndims      : dimension of dim (= 3 : lon,lat,time)                         !
! ncid       : netcdf file id                                                !
! varid1     : variable id of cloud fraction                                 !
! varid2     : variable id of clear fraction                                 !
! varid3     : variable id of fully attenuated fraction                      !
! varid4     : variable id of unclassify fraction                            !
! varid5     : variable id of nan fraction                                   !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : map_recvar2nc(monthisccplow,monthisccpmid,monthisccphigh,monthcolcloud!
!                    monthcolclear,dimidsm,file8,lonmax,latmax)              !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine map_recvar2nc2(low,mid,high,colcloud,colclear, &
                         dim,fname,nlon,nlat)

    use netcdf
    implicit none

    integer, parameter ::  ndims=3
    real,parameter  ::  nan=-9999
    integer  ::  nlon , nlat ,dim(ndims)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10
    character(LEN=*)  ::  fname
    !integer*4,dimension(nlon,nlat)  ::  tot,ret,retlow,retmid,rethigh
    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, colclear


    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

  
    call check(nf90_def_var(ncid, 'cllcalipso', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clmcalipso', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clhcalipso', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltcalipso', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clccalipso', NF90_FLOAT, dim, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO Total Clear Fraction'))
    call check(nf90_put_att(ncid, varid5, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))

 
    call check(nf90_enddef(ncid))


    call check(nf90_put_var(ncid, varid1, low))
    call check(nf90_put_var(ncid, varid2, mid))
    call check(nf90_put_var(ncid, varid3, high))
    call check(nf90_put_var(ncid, varid4, colcloud))
    call check(nf90_put_var(ncid, varid5, colclear))


    call check(nf90_close(ncid))

endsubroutine map_recvar2nc2
!----------------------------------------------------------------------------!

subroutine create_maphighnc(fname,dname,vlon,vlat,vtime,dim,grid,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5,grid*8
    integer, parameter ::  ndims = 3
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat,dim(ndims)
    integer  ::  lon_varid,lat_varid,pres_varid, time_varid,nc_id
    integer  ::  lon_dimid,lat_dimid,pres_dimid,time_dimid
	real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat


  
    call date_and_time(date,time,zone,value)

    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Map_High_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 
    call check(nf90_def_dim(nc_id, 'lon', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'lat', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
 
     dim = (/lon_dimid, lat_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'lon', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'lat', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

     call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
	
    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_maphighnc
!----------------------------------------------------------------------------!


subroutine maphigh(high,top,base,dim,fname,nlon,nlat)

    use netcdf
    implicit none

    integer, parameter ::  ndims=3
    real,parameter  ::  nan=-9999
    integer  ::  nlon , nlat ,dim(ndims)
    integer  ::  varid3,varid5, varid6, ncid !
    character(LEN=*)  ::  fname
    !integer*4,dimension(nlon,nlat)  ::  tot,ret,retlow,retmid,rethigh
    real*4,dimension(nlon,nlat)  ::  high,top,base


    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'clhcalipso', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'topcalipso', NF90_FLOAT, dim, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO Top Cloud Height'))
    call check(nf90_put_att(ncid, varid5, 'units','Kilometer'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))
 
    call check(nf90_def_var(ncid, 'bascalipso', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Base Cloud Height'))
    call check(nf90_put_att(ncid, varid6, 'units','Kilometer'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))
 
    call check(nf90_enddef(ncid))


    call check(nf90_put_var(ncid, varid3, high))
    call check(nf90_put_var(ncid, varid5, top))
    call check(nf90_put_var(ncid, varid6, base))

    call check(nf90_close(ncid))

endsubroutine maphigh


! call map_recvar2nc2phase(monthisccpliq,monthisccpice,monthisccpho,monthisccpun,monthisccpdust, &
!                          isccpdaypermonthlow,isccpdaypermonthmid,isccpdaypermonth,dimidsm,     &
!                          file8,lonmax-1,latmax-1)

subroutine map_recvar2nc2phase(liq,ice,dim,fname,nlon,nlat)

    use netcdf
    implicit none

    integer, parameter ::  ndims=3
    real,parameter  ::  nan=-9999
    integer  ::  nlon , nlat ,dim(ndims)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10
    character(LEN=*)  ::  fname
    !integer*4,dimension(nlon,nlat)  ::  tot,ret,retlow,retmid,rethigh
    real*4,dimension(nlon,nlat,4)  ::  liq,ice


    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

  
    call check(nf90_def_var(ncid, 'cll_liq', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clm_liq', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clh_liq', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clt_liq', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

   call check(nf90_def_var(ncid, 'cll_ice', NF90_FLOAT, dim, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid5, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clm_ice', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid6, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clh_ice', NF90_FLOAT, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid7, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clt_ice', NF90_FLOAT, dim, varid8))
    call check(nf90_put_att(ncid, varid8, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid8, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid8, '_FillValue',nan))


    call check(nf90_enddef(ncid))


    call check(nf90_put_var(ncid, varid1, liq(:,:,1)))
    call check(nf90_put_var(ncid, varid2, liq(:,:,2)))
    call check(nf90_put_var(ncid, varid3, liq(:,:,3)))
    call check(nf90_put_var(ncid, varid4, liq(:,:,4)))

    call check(nf90_put_var(ncid, varid5, ice(:,:,1)))
    call check(nf90_put_var(ncid, varid6, ice(:,:,2)))
    call check(nf90_put_var(ncid, varid7, ice(:,:,3)))
    call check(nf90_put_var(ncid, varid8, ice(:,:,4)))



    call check(nf90_close(ncid))

endsubroutine map_recvar2nc2phase
!----------------------------------------------------------------------------!

subroutine map_recvar2nc2phaseocc(liq,ice,indlow,indmid,indtot,dim,fname,nlon,nlat)

    use netcdf
    implicit none

    integer, parameter ::  ndims=3
    real,parameter  ::  nan=-9999
    integer  ::  nlon , nlat ,dim(ndims)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11
    character(LEN=*)  ::  fname
    !integer*4,dimension(nlon,nlat)  ::  tot,ret,retlow,retmid,rethigh
    real*4,dimension(nlon,nlat,4)  ::  liq,ice
    real*4,dimension(nlon,nlat)  ::  indtot,indlow,indmid



    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

  
    call check(nf90_def_var(ncid, 'cll_liq', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clm_liq', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clh_liq', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clt_liq', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

   call check(nf90_def_var(ncid, 'cll_ice', NF90_FLOAT, dim, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid5, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clm_ice', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid6, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clh_ice', NF90_FLOAT, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid7, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clt_ice', NF90_FLOAT, dim, varid8))
    call check(nf90_put_att(ncid, varid8, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid8, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid8, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'ind_low', NF90_FLOAT, dim, varid9))
    call check(nf90_put_att(ncid, varid9, 'lon_name',                        &
               'indice low'))
    call check(nf90_put_att(ncid, varid9, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid9, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'ind_mid', NF90_FLOAT, dim, varid10))
    call check(nf90_put_att(ncid, varid10, 'lon_name',                        &
               'indice mid'))
    call check(nf90_put_att(ncid, varid10, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid10, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'ind_tot', NF90_FLOAT, dim, varid11))
    call check(nf90_put_att(ncid, varid11, 'lon_name',                        &
               'indice high/tot'))
    call check(nf90_put_att(ncid, varid11, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid11, '_FillValue',nan))

    call check(nf90_enddef(ncid))


    call check(nf90_put_var(ncid, varid1, liq(:,:,1)))
    call check(nf90_put_var(ncid, varid2, liq(:,:,2)))
    call check(nf90_put_var(ncid, varid3, liq(:,:,3)))
    call check(nf90_put_var(ncid, varid4, liq(:,:,4)))

    call check(nf90_put_var(ncid, varid5, ice(:,:,1)))
    call check(nf90_put_var(ncid, varid6, ice(:,:,2)))
    call check(nf90_put_var(ncid, varid7, ice(:,:,3)))
    call check(nf90_put_var(ncid, varid8, ice(:,:,4)))

    call check(nf90_put_var(ncid, varid9, indlow(:,:)))
    call check(nf90_put_var(ncid, varid10, indmid(:,:)))
    call check(nf90_put_var(ncid, varid11, indtot(:,:)))


    call check(nf90_close(ncid))

endsubroutine map_recvar2nc2phaseocc
!----------------------------------------------------------------------------!


subroutine map_recvar2nc2phaseocc2(liq,ice,un2,phase,dim,dim2,fname,nlon,nlat)

!call map_recvar2nc2phaseocc2(monthisccpho,monthisccpun,monthisccpdust,isccpdaypermonthlow, &
!                             isccpdaypermonthmid,isccpdaypermonth,dimidsm,     &
!                             file8,lonmax-1,latmax-1)
!
    use netcdf
    implicit none

    integer, parameter ::  ndims=3, ndims2=4
    real,parameter  ::  nan=-9999.
    integer  ::  nlon , nlat ,dim(ndims),dim2(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid21,varid22,varid23,varid24,varid25, varid26,varid27, varid28
    integer  ::  varid8,varid9,varid10,varid11, varid12,varid13,varid14,varid15,varid16
    character(LEN=*)  ::  fname
    real*4,dimension(nlon,nlat,4)  ::  ice,liq,phase
    real*4,dimension(nlon,nlat,4,catmax)  ::  un2


    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

 
    call check(nf90_def_var(ncid, 'cllcalipso_liq', NF90_FLOAT, dim, varid21))
    call check(nf90_put_att(ncid, varid21, 'lon_name',                        &
               'Low Level liquid Cloud'))
    call check(nf90_put_att(ncid, varid21, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid21, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clmcalipso_liq', NF90_FLOAT, dim, varid22))
    call check(nf90_put_att(ncid, varid22, 'lon_name',                        &
               'Middle Level liq Cloud'))
    call check(nf90_put_att(ncid, varid22, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid22, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clhcalipso_liq', NF90_FLOAT, dim, varid23))
    call check(nf90_put_att(ncid, varid23, 'lon_name',                        &
               'High Level liq Cloud'))
    call check(nf90_put_att(ncid, varid23, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid23, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltcalipso_liq', NF90_FLOAT, dim, varid24))
    call check(nf90_put_att(ncid, varid24, 'lon_name',                        &
               'Total Column Level liq Cloud'))
    call check(nf90_put_att(ncid, varid24, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid24, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cllcalipso_ice', NF90_FLOAT, dim, varid25))
    call check(nf90_put_att(ncid, varid25, 'lon_name',                        &
               'Low-level ice Cloud'))
    call check(nf90_put_att(ncid, varid25, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid25, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clmcalipso_ice', NF90_FLOAT, dim, varid26))
    call check(nf90_put_att(ncid, varid26, 'lon_name',                        &
               'CALIPSO Mid-level ice Cloud'))
    call check(nf90_put_att(ncid, varid26, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid26, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clhcalipso_ice', NF90_FLOAT, dim, varid27))
    call check(nf90_put_att(ncid, varid27, 'lon_name',                        &
               'CALIPSO High-level ice Cloud'))
    call check(nf90_put_att(ncid, varid27, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid27, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltcalipso_ice', NF90_FLOAT, dim, varid28))
    call check(nf90_put_att(ncid, varid28, 'lon_name',                        &
               'CALIPSO Total ice Cloud'))
    call check(nf90_put_att(ncid, varid28, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid28, '_FillValue',nan))



 
    call check(nf90_def_var(ncid, 'cllcalipso_un', NF90_FLOAT, dim2, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'Low Level unclassified Cloud'))
    call check(nf90_put_att(ncid, varid1, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clmcalipso_un', NF90_FLOAT, dim2, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'Middle Level unclassified Cloud'))
    call check(nf90_put_att(ncid, varid2, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clhcalipso_un', NF90_FLOAT, dim2, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'High Level unclassified Cloud'))
    call check(nf90_put_att(ncid, varid3, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltcalipso_un', NF90_FLOAT, dim2, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'Total Column Level unclassified Cloud'))
    call check(nf90_put_att(ncid, varid4, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

!!$    call check(nf90_def_var(ncid, 'cll_ho', NF90_FLOAT, dim, varid5))
!!$    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
!!$               'Low-level HO Cloud'))
!!$    call check(nf90_put_att(ncid, varid5, 'units','Fraction'))
!!$    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))
!!$
!!$    call check(nf90_def_var(ncid, 'clm_ho', NF90_FLOAT, dim, varid6))
!!$    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
!!$               'CALIPSO Mid-level HO Cloud'))
!!$    call check(nf90_put_att(ncid, varid6, 'units','Fraction'))
!!$    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))
!!$
!!$    call check(nf90_def_var(ncid, 'clh_ho', NF90_FLOAT, dim, varid7))
!!$    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
!!$               'CALIPSO High-level HO Cloud'))
!!$    call check(nf90_put_att(ncid, varid7, 'units','Fraction'))
!!$    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))
!!$
!!$    call check(nf90_def_var(ncid, 'clt_ho', NF90_FLOAT, dim, varid8))
!!$    call check(nf90_put_att(ncid, varid8, 'lon_name',                        &
!!$               'CALIPSO Total HO Cloud'))
!!$    call check(nf90_put_att(ncid, varid8, 'units','Fraction'))
!!$    call check(nf90_put_att(ncid, varid8, '_FillValue',nan))
!!$
!!$    call check(nf90_def_var(ncid, 'cll_dust', NF90_FLOAT, dim, varid12))
!!$    call check(nf90_put_att(ncid, varid12, 'lon_name',                        &
!!$               'CALIPSO Low-level dust Cloud'))
!!$    call check(nf90_put_att(ncid, varid12, 'units','Fraction'))
!!$    call check(nf90_put_att(ncid, varid12, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cllcalipso_RPIC', NF90_FLOAT, dim, varid13))
    call check(nf90_put_att(ncid, varid13, 'lon_name',                        &
               'CALIPSO Low-level Relative Percentage of Ice in Cloud'))
    call check(nf90_put_att(ncid, varid13, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid13, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clmcalipso_RPIC', NF90_FLOAT, dim, varid14))
    call check(nf90_put_att(ncid, varid14, 'lon_name',                        &
               'CALIPSO Mid-level Relative Percentage of Ice in Cloud'))
    call check(nf90_put_att(ncid, varid14, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid14, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clhcalipso_RPIC', NF90_FLOAT, dim, varid15))
    call check(nf90_put_att(ncid, varid15, 'lon_name',                        &
               'CALIPSO High-level Relative Percentage of Ice in Cloud'))
    call check(nf90_put_att(ncid, varid15, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid15, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltcalipso_RPIC', NF90_FLOAT, dim, varid16))
    call check(nf90_put_att(ncid, varid16, 'lon_name',                        &
               'CALIPSO Total Relative Percentage of Ice in Cloud'))
    call check(nf90_put_att(ncid, varid16, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid16, '_FillValue',nan))



!!$    call check(nf90_def_var(ncid, 'ind_low', NF90_FLOAT, dim, varid9))
!!$    call check(nf90_put_att(ncid, varid9, 'lon_name',                        &
!!$               'indice low'))
!!$    call check(nf90_put_att(ncid, varid9, 'units','Fraction'))
!!$    call check(nf90_put_att(ncid, varid9, '_FillValue',nan))
!!$
!!$    call check(nf90_def_var(ncid, 'ind_mid', NF90_FLOAT, dim, varid10))
!!$    call check(nf90_put_att(ncid, varid10, 'lon_name',                        &
!!$               'indice mid'))
!!$    call check(nf90_put_att(ncid, varid10, 'units','Fraction'))
!!$    call check(nf90_put_att(ncid, varid10, '_FillValue',nan))
!!$
!!$    call check(nf90_def_var(ncid, 'ind_tot', NF90_FLOAT, dim, varid11))
!!$    call check(nf90_put_att(ncid, varid11, 'lon_name',                        &
!!$               'indice high/tot'))
!!$    call check(nf90_put_att(ncid, varid11, 'units','Fraction'))
!!$    call check(nf90_put_att(ncid, varid11, '_FillValue',nan))
    call check(nf90_enddef(ncid))




    call check(nf90_put_var(ncid, varid21, liq(:,:,1)))
    call check(nf90_put_var(ncid, varid22, liq(:,:,2)))
    call check(nf90_put_var(ncid, varid23, liq(:,:,3)))
    call check(nf90_put_var(ncid, varid24, liq(:,:,4)))

    call check(nf90_put_var(ncid, varid25, ice(:,:,1)))
    call check(nf90_put_var(ncid, varid26, ice(:,:,2)))
    call check(nf90_put_var(ncid, varid27, ice(:,:,3)))
    call check(nf90_put_var(ncid, varid28, ice(:,:,4)))

    call check(nf90_put_var(ncid, varid1, un2(:,:,1,:)))
    call check(nf90_put_var(ncid, varid2, un2(:,:,2,:)))
    call check(nf90_put_var(ncid, varid3, un2(:,:,3,:)))
    call check(nf90_put_var(ncid, varid4, un2(:,:,4,:)))

    call check(nf90_put_var(ncid, varid13, phase(:,:,1)))
    call check(nf90_put_var(ncid, varid14, phase(:,:,2)))
    call check(nf90_put_var(ncid, varid15, phase(:,:,3)))
    call check(nf90_put_var(ncid, varid16, phase(:,:,4)))

 

!    call check(nf90_put_var(ncid, varid9, indlow(:,:)))
!    call check(nf90_put_var(ncid, varid10, indmid(:,:)))
!    call check(nf90_put_var(ncid, varid11, indtot(:,:)))


    call check(nf90_close(ncid))

endsubroutine map_recvar2nc2phaseocc2
!----------------------------------------------------------------------------!


subroutine map_recvar2nc3(low,mid,high,colcloud,colclear,height,indtot,&
                         hlow,hmid,hhigh,hcol,hheight,dim,dim2,dim3,fname,nlon,nlat)

! h_CA utile car 31 valeur au max...
    use netcdf
    implicit none

    integer, parameter ::  ndims=3,histmax=10,histmax2=40,ndims2=4
    real,parameter  ::  nan=-999
    integer  ::  nlon , nlat ,dim(ndims),ihist,dim2(ndims2),dim3(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11,varid12,varid13,varid14
    integer  ::  varid15,varid16,varid17
    character(LEN=*)  ::  fname
    !integer*4,dimension(nlon,nlat)  ::  tot,ret,retlow,retmid,rethigh
    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, colclear, height
    integer,dimension(nlon,nlat)  ::  indtot,f_CA,f_CAL,f_CAM,f_CAH,f_CZ
    integer,dimension(nlon,nlat,histmax)  ::  hlow,hmid,hhigh,hcol
    real*4,dimension(nlon,nlat,histmax2)  ::  hheight    
    integer,dimension(histmax2)  ::  histmod2
    real,dimension(histmax)  ::  histmod


histmod(:)=0;
do ihist=1,histmax
   histmod(ihist+1)=histmod(ihist)+0.1
enddo

histmod2(:)=0;
do ihist=1,histmax2
   histmod2(ihist+1)=histmod2(ihist)+1
enddo


!print *, tot(5,5),ret(5,5)

do ilon=1,nlon
  do ilat=1,nlat

if (colcloud(ilon,ilat).eq.0)then
f_CA(ilon,ilat)=0
elseif(colcloud(ilon,ilat).eq.-999)then
f_CA(ilon,ilat)=-999
else
f_CA(ilon,ilat)=100
endif

if (low(ilon,ilat).eq.0)then
f_CAL(ilon,ilat)=0
elseif(low(ilon,ilat).eq.-999)then
f_CAL(ilon,ilat)=-999
else
f_CAL(ilon,ilat)=100
endif

if (mid(ilon,ilat).eq.0)then
f_CAM(ilon,ilat)=0
elseif(mid(ilon,ilat).eq.-999)then
f_CAM(ilon,ilat)=-999
else
f_CAM(ilon,ilat)=100
endif

if (high(ilon,ilat).eq.0)then
f_CAH(ilon,ilat)=0
elseif(high(ilon,ilat).eq.-999)then
f_CAH(ilon,ilat)=-999
else
f_CAH(ilon,ilat)=100
endif

if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif


enddo
enddo

!!$do ilon=1,nlon
!!$  do ilat=1,nlat
!!$     do ihist=1,ihistmax-1
!!$
!!$        if(h_CA(ilon,ilat,ihist)
!!$enddo
!!$enddo
!!$enddo

    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'n_tot', NF90_INT4, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'Number of orbit passages'))
    call check(nf90_put_att(ncid, varid7, 'units','Arbitrary Unit'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))


    call check(nf90_def_var(ncid, 'a_CAL', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAM', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAH', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CA', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CZ', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Cloud Top Height level'))
    call check(nf90_put_att(ncid, varid6, 'units','1-40'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))
 
    call check(nf90_def_var(ncid, 'f_CAL', NF90_INT4, dim, varid8))
    call check(nf90_def_var(ncid, 'f_CAM', NF90_INT4, dim, varid9))
    call check(nf90_def_var(ncid, 'f_CAH', NF90_INT4, dim, varid10))
    call check(nf90_def_var(ncid, 'f_CA', NF90_INT4, dim, varid11))
    call check(nf90_def_var(ncid, 'f_CZ', NF90_INT4, dim, varid12))

    call check(nf90_def_var(ncid, 'h_CAL', NF90_INT4, dim2, varid13))
    call check(nf90_def_var(ncid, 'h_CAM', NF90_INT4, dim2, varid14))
    call check(nf90_def_var(ncid, 'h_CAH', NF90_INT4, dim2, varid15))
    call check(nf90_def_var(ncid, 'h_CA', NF90_INT4, dim2, varid16))
    call check(nf90_def_var(ncid, 'h_CZ', NF90_FLOAT, dim3, varid17))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid1, low))
    call check(nf90_put_var(ncid, varid2, mid))
    call check(nf90_put_var(ncid, varid3, high))
    call check(nf90_put_var(ncid, varid4, colcloud))
 !   call check(nf90_put_var(ncid, varid5, colclear))
    call check(nf90_put_var(ncid, varid6, height)) 
    call check(nf90_put_var(ncid, varid7, indtot)) 
    call check(nf90_put_var(ncid, varid8, f_CAL)) 
 call check(nf90_put_var(ncid, varid9, f_CAM)) 
 call check(nf90_put_var(ncid, varid10, f_CAH)) 
 call check(nf90_put_var(ncid, varid11, f_CA)) 
 call check(nf90_put_var(ncid, varid12, f_CZ))
 call check(nf90_put_var(ncid, varid13, hlow ))
 call check(nf90_put_var(ncid, varid14, hmid )) 
 call check(nf90_put_var(ncid, varid15, hhigh)) 
 call check(nf90_put_var(ncid, varid16, hcol)) 
 call check(nf90_put_var(ncid, varid17, hheight))
 

    call check(nf90_close(ncid))

endsubroutine map_recvar2nc3
!----------------------------------------------------------------------------!

subroutine map_recvar2nc7(low,mid,high,colcloud,height,indtot,&
                         hlow,hmid,hhigh,hcol,hheight,dim,dim2,dim3,dim4,fname,nlon,nlat, &
                         lowtemp,midtemp,hightemp,coltemp,hlowtemp,hmidtemp,hhightemp, &
                         hcoltemp)
! h_CA utile car 31 v!aleur au max...
    use netcdf
    implicit none

    integer, parameter ::  ndims=3,histmax=10,histmax2=40,ndims2=4, histmax3=28
    real,parameter  ::  nan=-999
    integer  ::  nlon , nlat ,dim(ndims),ihist,dim2(ndims2),dim3(ndims2),dim4(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11,varid12,varid13,varid14
    integer  ::  varid15,varid16,varid17
    character(LEN=*)  ::  fname
    !integer*4,dimension(nlon,nlat)  ::  tot,ret,retlow,retmid,rethigh
    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, height
    integer,dimension(nlon,nlat)  ::  indtot,f_CA,f_CAL,f_CAM,f_CAH,f_CZ
    integer,dimension(nlon,nlat)  ::  f_CT,f_CTL,f_CTM,f_CTH

    integer,dimension(nlon,nlat,histmax)  ::  hlow,hmid,hhigh,hcol
    real*4,dimension(nlon,nlat,histmax2)  ::  hheight    
    integer,dimension(histmax2)  ::  histmod2
    real,dimension(histmax)  ::  histmod
    integer,dimension(nlon,nlat,histmax3)  ::  hlowtemp,hmidtemp,hhightemp,hcoltemp

    real*4,dimension(nlon,nlat)  ::  lowtemp,midtemp,hightemp,coltemp
    integer  ::  varid27, varid28, varid29, varid30, varid18, varid19, varid20, varid21
    integer  ::  varid31, varid32, varid33, varid34


histmod(:)=0;
do ihist=1,histmax
   histmod(ihist+1)=histmod(ihist)+0.1
enddo

histmod2(:)=0;
do ihist=1,histmax2
   histmod2(ihist+1)=histmod2(ihist)+1
enddo


do ilon=1,nlon
  do ilat=1,nlat

if (colcloud(ilon,ilat).eq.0)then
f_CA(ilon,ilat)=0
elseif(colcloud(ilon,ilat).eq.-999)then
f_CA(ilon,ilat)=-999
else
f_CA(ilon,ilat)=100
endif

if (low(ilon,ilat).eq.0)then
f_CAL(ilon,ilat)=0
elseif(low(ilon,ilat).eq.-999)then
f_CAL(ilon,ilat)=-999
else
f_CAL(ilon,ilat)=100
endif

if (mid(ilon,ilat).eq.0)then
f_CAM(ilon,ilat)=0
elseif(mid(ilon,ilat).eq.-999)then
f_CAM(ilon,ilat)=-999
else
f_CAM(ilon,ilat)=100
endif

if (high(ilon,ilat).eq.0)then
f_CAH(ilon,ilat)=0
elseif(high(ilon,ilat).eq.-999)then
f_CAH(ilon,ilat)=-999
else
f_CAH(ilon,ilat)=100
endif

if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif

if (lowtemp(ilon,ilat).gt.0)then
f_CTL(ilon,ilat)=100
else
f_CTL(ilon,ilat)=0
endif

if (midtemp(ilon,ilat).gt.0)then
f_CTM(ilon,ilat)=100
else
f_CTM(ilon,ilat)=0
endif

if (hightemp(ilon,ilat).gt.0)then
f_CTH(ilon,ilat)=100
else
f_CTH(ilon,ilat)=0
endif


if (coltemp(ilon,ilat).gt.0)then
f_CT(ilon,ilat)=100
else
f_CT(ilon,ilat)=0
endif

enddo
enddo

!!$do ilon=1,nlon
!!$  do ilat=1,nlat
!!$     do ihist=1,ihistmax-1
!!$
!!$        if(h_CA(ilon,ilat,ihist)
!!$enddo
!!$enddo
!!$enddo

    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'n_tot', NF90_INT4, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'Number of orbit passages'))
    call check(nf90_put_att(ncid, varid7, 'units','Arbitrary Unit'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))


    call check(nf90_def_var(ncid, 'a_CAL', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAM', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAH', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CA', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CZ', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Cloud Top Height level'))
    call check(nf90_put_att(ncid, varid6, 'units','1-40'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))
 
    call check(nf90_def_var(ncid, 'f_CAL', NF90_INT4, dim, varid8))
    call check(nf90_def_var(ncid, 'f_CAM', NF90_INT4, dim, varid9))
    call check(nf90_def_var(ncid, 'f_CAH', NF90_INT4, dim, varid10))
    call check(nf90_def_var(ncid, 'f_CA', NF90_INT4, dim, varid11))
    call check(nf90_def_var(ncid, 'f_CZ', NF90_INT4, dim, varid12))

    call check(nf90_def_var(ncid, 'h_CAL', NF90_INT4, dim2, varid13))
    call check(nf90_def_var(ncid, 'h_CAM', NF90_INT4, dim2, varid14))
    call check(nf90_def_var(ncid, 'h_CAH', NF90_INT4, dim2, varid15))
    call check(nf90_def_var(ncid, 'h_CA', NF90_INT4, dim2, varid16))
    call check(nf90_def_var(ncid, 'h_CZ', NF90_FLOAT, dim3, varid17))

    call check(nf90_def_var(ncid, 'a_CTL', NF90_FLOAT, dim, varid27))
    call check(nf90_def_var(ncid, 'a_CTM', NF90_FLOAT, dim, varid28))
    call check(nf90_def_var(ncid, 'a_CTH', NF90_FLOAT, dim, varid29))
    call check(nf90_def_var(ncid, 'a_CT', NF90_FLOAT, dim, varid30))

    call check(nf90_def_var(ncid, 'f_CTL', NF90_INT4, dim, varid18))
    call check(nf90_def_var(ncid, 'f_CTM', NF90_INT4, dim, varid19))
    call check(nf90_def_var(ncid, 'f_CTH', NF90_INT4, dim, varid20))
    call check(nf90_def_var(ncid, 'f_CT', NF90_INT4, dim, varid21))

    call check(nf90_def_var(ncid, 'h_CTL', NF90_INT4, dim4, varid31))
    call check(nf90_def_var(ncid, 'h_CTM', NF90_INT4, dim4, varid32))
    call check(nf90_def_var(ncid, 'h_CTH', NF90_INT4, dim4, varid33))
    call check(nf90_def_var(ncid, 'h_CT', NF90_INT4, dim4, varid34))


    call check(nf90_enddef(ncid))

  
    call check(nf90_put_var(ncid, varid1, low))
    call check(nf90_put_var(ncid, varid2, mid))
    call check(nf90_put_var(ncid, varid3, high))
    call check(nf90_put_var(ncid, varid4, colcloud))
    call check(nf90_put_var(ncid, varid6, height)) 
    call check(nf90_put_var(ncid, varid7, indtot)) 
    call check(nf90_put_var(ncid, varid8, f_CAL)) 
    call check(nf90_put_var(ncid, varid9, f_CAM)) 
    call check(nf90_put_var(ncid, varid10, f_CAH)) 
    call check(nf90_put_var(ncid, varid11, f_CA)) 
    call check(nf90_put_var(ncid, varid12, f_CZ))
    call check(nf90_put_var(ncid, varid13, hlow ))
    call check(nf90_put_var(ncid, varid14, hmid )) 
    call check(nf90_put_var(ncid, varid15, hhigh)) 
    call check(nf90_put_var(ncid, varid16, hcol)) 
    call check(nf90_put_var(ncid, varid17, hheight))
    call check(nf90_put_var(ncid, varid18, f_CTL)) 
    call check(nf90_put_var(ncid, varid19, f_CTM)) 
    call check(nf90_put_var(ncid, varid20, f_CTH)) 
    call check(nf90_put_var(ncid, varid21, f_CT)) 
    call check(nf90_put_var(ncid, varid27, lowtemp))
    call check(nf90_put_var(ncid, varid28, midtemp))
    call check(nf90_put_var(ncid, varid29, hightemp))
    call check(nf90_put_var(ncid, varid30, coltemp))
    call check(nf90_put_var(ncid, varid31, hlowtemp))
    call check(nf90_put_var(ncid, varid32, hmidtemp))
    call check(nf90_put_var(ncid, varid33, hhightemp))
    call check(nf90_put_var(ncid, varid34, hcoltemp))

 

    call check(nf90_close(ncid))

endsubroutine map_recvar2nc7



subroutine map_recvar2nc6(low,mid,high,colcloud,height,indtot,&
                         hlow,hmid,hhigh,hcol,hheight,dim,dim2,dim3,fname,nlon,nlat)!,&
!                         lowtemp,midtemp! &
!                         hightemp,coltemp,hlowtemp,hmidtemp,hhightemp, &
!                         hcoltemp,dim4)

! h_CA utile car 31 valeur au max...
    use netcdf
    implicit none

    integer, parameter ::  ndims=3,histmax=10,histmax2=40,ndims2=4,histmax3=28
    real,parameter  ::  nan=-999
    integer  ::  nlon , nlat ,dim(ndims),ihist,dim2(ndims2),dim3(ndims2),dim4(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11,varid12,varid13,varid14
    integer  ::  varid15,varid16,varid17
    character(LEN=*)  ::  fname
    !integer*4,dimension(nlon,nlat)  ::  tot,ret,retlow,retmid,rethigh
    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, height
    integer,dimension(nlon,nlat)  ::  indtot,f_CA,f_CAL,f_CAM,f_CAH,f_CZ
    integer,dimension(nlon,nlat,histmax)  ::  hlow,hmid,hhigh,hcol
    integer,dimension(nlon,nlat,histmax2)  ::  hheight    
    integer,dimension(histmax2)  ::  histmod2
    real,dimension(histmax)  ::  histmod

   
    real*4,dimension(nlon,nlat)  ::  lowtemp,midtemp,hightemp,coltemp
    integer,dimension(nlon,nlat)  ::  f_CT,f_CTL,f_CTM,f_CTH
    integer,dimension(nlon,nlat,histmax3)  ::  hlowtemp,hmidtemp,hhightemp,hcoltemp
    integer,dimension(histmax3)  ::  histmod3





histmod(:)=0;
do ihist=1,histmax
   histmod(ihist+1)=histmod(ihist)+0.1
enddo

histmod2(:)=0;
do ihist=1,histmax2
   histmod2(ihist+1)=histmod2(ihist)+1
enddo


!print *, tot(5,5),ret(5,5)

do ilon=1,nlon
  do ilat=1,nlat

if (colcloud(ilon,ilat).eq.0)then
f_CA(ilon,ilat)=0
elseif(colcloud(ilon,ilat).eq.-999)then
f_CA(ilon,ilat)=-999
else
f_CA(ilon,ilat)=100
endif

if (low(ilon,ilat).eq.0)then
f_CAL(ilon,ilat)=0
elseif(low(ilon,ilat).eq.-999)then
f_CAL(ilon,ilat)=-999
else
f_CAL(ilon,ilat)=100
endif

if (mid(ilon,ilat).eq.0)then
f_CAM(ilon,ilat)=0
elseif(mid(ilon,ilat).eq.-999)then
f_CAM(ilon,ilat)=-999
else
f_CAM(ilon,ilat)=100
endif

if (high(ilon,ilat).eq.0)then
f_CAH(ilon,ilat)=0
elseif(high(ilon,ilat).eq.-999)then
f_CAH(ilon,ilat)=-999
else
f_CAH(ilon,ilat)=100
endif

if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif


enddo
enddo

!!$do ilon=1,nlon
!!$  do ilat=1,nlat
!!$     do ihist=1,ihistmax-1
!!$
!!$        if(h_CA(ilon,ilat,ihist)
!!$enddo
!!$enddo
!!$enddo

    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'n_tot', NF90_INT4, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'Number of orbit passages'))
    call check(nf90_put_att(ncid, varid7, 'units','Arbitrary Unit'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))


    call check(nf90_def_var(ncid, 'a_CAL', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAM', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAH', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CA', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CZ', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Cloud Top Height level'))
    call check(nf90_put_att(ncid, varid6, 'units','1-40'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))
 
call check(nf90_def_var(ncid, 'f_CAL', NF90_INT4, dim, varid8))
call check(nf90_def_var(ncid, 'f_CAM', NF90_INT4, dim, varid9))
call check(nf90_def_var(ncid, 'f_CAH', NF90_INT4, dim, varid10))
call check(nf90_def_var(ncid, 'f_CA', NF90_INT4, dim, varid11))
call check(nf90_def_var(ncid, 'f_CZ', NF90_INT4, dim, varid12))

call check(nf90_def_var(ncid, 'h_CAL', NF90_INT4, dim2, varid13))
call check(nf90_def_var(ncid, 'h_CAM', NF90_INT4, dim2, varid14))
call check(nf90_def_var(ncid, 'h_CAH', NF90_INT4, dim2, varid15))
call check(nf90_def_var(ncid, 'h_CA', NF90_INT4, dim2, varid16))
call check(nf90_def_var(ncid, 'h_CZ', NF90_INT4, dim3, varid17))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid1, low))
    call check(nf90_put_var(ncid, varid2, mid))
    call check(nf90_put_var(ncid, varid3, high))
    call check(nf90_put_var(ncid, varid4, colcloud))
    call check(nf90_put_var(ncid, varid6, height)) 
    call check(nf90_put_var(ncid, varid7, indtot)) 
    call check(nf90_put_var(ncid, varid8, f_CAL)) 
 call check(nf90_put_var(ncid, varid9, f_CAM)) 
 call check(nf90_put_var(ncid, varid10, f_CAH)) 
 call check(nf90_put_var(ncid, varid11, f_CA)) 
 call check(nf90_put_var(ncid, varid12, f_CZ))
 call check(nf90_put_var(ncid, varid13, hlow ))
 call check(nf90_put_var(ncid, varid14, hmid )) 
 call check(nf90_put_var(ncid, varid15, hhigh)) 
 call check(nf90_put_var(ncid, varid16, hcol)) 
 call check(nf90_put_var(ncid, varid17, hheight))
 

    call check(nf90_close(ncid))

endsubroutine map_recvar2nc6
!----------------------------------------------------------------------------!











subroutine map_recvar2nc4(low,mid,high,colcloud,colclear,height,indtot,&
                         hlow,hmid,hhigh,hcol,hheight,lowtemp,midtemp, &
                         hightemp,coltemp,hlowtemp,hmidtemp,hhightemp, &
                         hcoltemp,dim,dim2,dim3,dim4,fname,nlon,nlat)

! h_CA utile car 31 valeur au max...
    use netcdf
    implicit none

    integer, parameter ::  ndims=3,histmax=10,histmax2=40,ndims2=4,histmax3=28
    real,parameter  ::  nan=-999
    integer  ::  nlon , nlat ,dim(ndims),ihist,dim2(ndims2),dim3(ndims2),dim4(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11,varid12,varid13,varid14
    integer  ::  varid15,varid16,varid17,varid18,varid19,varid20,varid21,varid22,varid23,varid24,varid25,varid26,varid27,varid28,varid29,varid30
    character(LEN=*)  ::  fname
    !integer*4,dimension(nlon,nlat)  ::  tot,ret,retlow,retmid,rethigh
    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, colclear, height
    real*4,dimension(nlon,nlat)  ::  lowtemp,midtemp,hightemp,coltemp

    integer,dimension(nlon,nlat)  ::  indtot,f_CA,f_CAL,f_CAM,f_CAH,f_CZ
    integer,dimension(nlon,nlat)  ::  f_CT,f_CTL,f_CTM,f_CTH
    integer,dimension(nlon,nlat,histmax)  ::  hlow,hmid,hhigh,hcol
    integer,dimension(nlon,nlat,histmax2)  ::  hheight   
    integer,dimension(nlon,nlat,histmax3)  ::  hlowtemp,hmidtemp,hhightemp,hcoltemp

    integer,dimension(histmax2)  ::  histmod2
    real,dimension(histmax)  ::  histmod
    integer,dimension(histmax3)  ::  histmod3


histmod(:)=0;
do ihist=1,histmax
   histmod(ihist+1)=histmod(ihist)+0.1
enddo

histmod2(:)=0;
do ihist=1,histmax2
   histmod2(ihist+1)=histmod2(ihist)+1
enddo



print *, 'f_ begin'

!print *, tot(5,5),ret(5,5)

do ilon=1,nlon
  do ilat=1,nlat

if (colcloud(ilon,ilat).eq.0)then
f_CA(ilon,ilat)=0
elseif(colcloud(ilon,ilat).eq.-999)then
f_CA(ilon,ilat)=-999
else
f_CA(ilon,ilat)=100
endif

if (low(ilon,ilat).eq.0)then
f_CAL(ilon,ilat)=0
elseif(low(ilon,ilat).eq.-999)then
f_CAL(ilon,ilat)=-999
else
f_CAL(ilon,ilat)=100
endif

if (mid(ilon,ilat).eq.0)then
f_CAM(ilon,ilat)=0
elseif(mid(ilon,ilat).eq.-999)then
f_CAM(ilon,ilat)=-999
else
f_CAM(ilon,ilat)=100
endif

if (high(ilon,ilat).eq.0)then
f_CAH(ilon,ilat)=0
elseif(high(ilon,ilat).eq.-999)then
f_CAH(ilon,ilat)=-999
else
f_CAH(ilon,ilat)=100
endif

if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif

if (lowtemp(ilon,ilat).gt.0)then
f_CTL(ilon,ilat)=100
else
f_CTL(ilon,ilat)=0
endif

if (midtemp(ilon,ilat).gt.0)then
f_CTM(ilon,ilat)=100
else
f_CTM(ilon,ilat)=0
endif

if (hightemp(ilon,ilat).gt.0)then
f_CTH(ilon,ilat)=100
else
f_CTH(ilon,ilat)=0
endif


if (coltemp(ilon,ilat).gt.0)then
f_CT(ilon,ilat)=100
else
f_CT(ilon,ilat)=0
endif


enddo
enddo

print *, 'f_ ok'


    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

     call check(nf90_def_var(ncid, 'a_CTL', NF90_FLOAT, dim, varid27))
     call check(nf90_def_var(ncid, 'a_CTM', NF90_FLOAT, dim, varid28))
     call check(nf90_def_var(ncid, 'a_CTH', NF90_FLOAT, dim, varid29))
     call check(nf90_def_var(ncid, 'a_CT', NF90_FLOAT, dim, varid30))

call check(nf90_def_var(ncid, 'f_CTL', NF90_INT4, dim, varid18))
call check(nf90_def_var(ncid, 'f_CTM', NF90_INT4, dim, varid19))
call check(nf90_def_var(ncid, 'f_CTH', NF90_INT4, dim, varid20))
call check(nf90_def_var(ncid, 'f_CT', NF90_INT4, dim, varid21))

call check(nf90_def_var(ncid, 'h_CTL', NF90_INT4, dim4, varid23))
call check(nf90_def_var(ncid, 'h_CTM', NF90_INT4, dim4, varid24))
call check(nf90_def_var(ncid, 'h_CTH', NF90_INT4, dim4, varid25))
call check(nf90_def_var(ncid, 'h_CT', NF90_INT4, dim4, varid26))


print *, 'def var ok'

    call check(nf90_enddef(ncid))
print *, 'f_ ok'

 call check(nf90_put_var(ncid, varid18, f_CTL)) 
 call check(nf90_put_var(ncid, varid19, f_CTM)) 
 call check(nf90_put_var(ncid, varid20, f_CTH)) 
 call check(nf90_put_var(ncid, varid21, f_CT)) 
print *, 'f_ ok'

 call check(nf90_put_var(ncid, varid23, hlowtemp ))
 call check(nf90_put_var(ncid, varid24, hmidtemp )) 
 call check(nf90_put_var(ncid, varid25, hhightemp)) 
 call check(nf90_put_var(ncid, varid26, hcoltemp)) 
print *, 'f_ ok'

    call check(nf90_put_var(ncid, varid27, lowtemp))
    call check(nf90_put_var(ncid, varid28, midtemp))
    call check(nf90_put_var(ncid, varid29, hightemp))
    call check(nf90_put_var(ncid, varid30, coltemp))



print *, 'f_ ok'

 

    call check(nf90_close(ncid))

endsubroutine map_recvar2nc4
!----------------------------------------------------------------------------!


subroutine map_recvar2nc5(low,mid,high,colcloud,height,indtot,&
                         hlow,hmid,hhigh,hcol,hheight,lowtemp,midtemp, &
                         hightemp,coltemp,hlowtemp,hmidtemp,hhightemp, &
                         hcoltemp,dim,dim2,dim3,dim4,fname,nlon,nlat)

! h_CA utile car 31 valeur au max...
    use netcdf
    implicit none

    integer, parameter ::  ndims=3,histmax=10,histmax2=40,ndims2=4,histmax3=28
    real,parameter  ::  nan=-999
    integer  ::  nlon , nlat ,dim(ndims),ihist,dim2(ndims2),dim3(ndims2),dim4(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11,varid12,varid13,varid14
    integer  ::  varid15,varid16,varid17
    character(LEN=*)  ::  fname
    !integer*4,dimension(nlon,nlat)  ::  tot,ret,retlow,retmid,rethigh
    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, height
    integer,dimension(nlon,nlat)  ::  indtot,f_CA,f_CAL,f_CAM,f_CAH,f_CZ
    integer,dimension(nlon,nlat,histmax)  ::  hlow,hmid,hhigh,hcol
    integer,dimension(nlon,nlat,histmax2)  ::  hheight    
    integer,dimension(histmax2)  ::  histmod2
    real,dimension(histmax)  ::  histmod

    integer  ::  varid18,varid19,varid20,varid21,varid22,varid23,varid24,varid25,varid26,varid27,varid28,varid29,varid30

    real*4,dimension(nlon,nlat)  ::  lowtemp,midtemp,hightemp,coltemp

    integer,dimension(nlon,nlat)  ::  f_CT,f_CTL,f_CTM,f_CTH
    integer,dimension(nlon,nlat,histmax3)  ::  hlowtemp,hmidtemp,hhightemp,hcoltemp

histmod(:)=0;
do ihist=1,histmax
   histmod(ihist+1)=histmod(ihist)+0.1
enddo

histmod2(:)=0;
do ihist=1,histmax2
   histmod2(ihist+1)=histmod2(ihist)+1
enddo


!print *, tot(5,5),ret(5,5)

do ilon=1,nlon
  do ilat=1,nlat

if (colcloud(ilon,ilat).eq.0)then
f_CA(ilon,ilat)=0
elseif(colcloud(ilon,ilat).eq.-999)then
f_CA(ilon,ilat)=-999
else
f_CA(ilon,ilat)=100
endif

if (low(ilon,ilat).eq.0)then
f_CAL(ilon,ilat)=0
elseif(low(ilon,ilat).eq.-999)then
f_CAL(ilon,ilat)=-999
else
f_CAL(ilon,ilat)=100
endif

if (mid(ilon,ilat).eq.0)then
f_CAM(ilon,ilat)=0
elseif(mid(ilon,ilat).eq.-999)then
f_CAM(ilon,ilat)=-999
else
f_CAM(ilon,ilat)=100
endif

if (high(ilon,ilat).eq.0)then
f_CAH(ilon,ilat)=0
elseif(high(ilon,ilat).eq.-999)then
f_CAH(ilon,ilat)=-999
else
f_CAH(ilon,ilat)=100
endif

if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif


if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif

if (lowtemp(ilon,ilat).gt.0)then
f_CTL(ilon,ilat)=100
else
f_CTL(ilon,ilat)=0
endif

if (midtemp(ilon,ilat).gt.0)then
f_CTM(ilon,ilat)=100
else
f_CTM(ilon,ilat)=0
endif

if (hightemp(ilon,ilat).gt.0)then
f_CTH(ilon,ilat)=100
else
f_CTH(ilon,ilat)=0
endif


if (coltemp(ilon,ilat).gt.0)then
f_CT(ilon,ilat)=100
else
f_CT(ilon,ilat)=0
endif


enddo
enddo

print *, 'creation fichier'

    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'n_tot', NF90_INT4, dim, varid7))
    call check(nf90_def_var(ncid, 'a_CAL', NF90_FLOAT, dim, varid1))
    call check(nf90_def_var(ncid, 'a_CAM', NF90_FLOAT, dim, varid2))
    call check(nf90_def_var(ncid, 'a_CAH', NF90_FLOAT, dim, varid3))
    call check(nf90_def_var(ncid, 'a_CA', NF90_FLOAT, dim, varid4))
    call check(nf90_def_var(ncid, 'a_CZ', NF90_FLOAT, dim, varid6))
 
call check(nf90_def_var(ncid, 'f_CAL', NF90_INT4, dim, varid8))
call check(nf90_def_var(ncid, 'f_CAM', NF90_INT4, dim, varid9))
call check(nf90_def_var(ncid, 'f_CAH', NF90_INT4, dim, varid10))
call check(nf90_def_var(ncid, 'f_CA', NF90_INT4, dim, varid11))
call check(nf90_def_var(ncid, 'f_CZ', NF90_INT4, dim, varid12))

call check(nf90_def_var(ncid, 'h_CAL', NF90_INT4, dim2, varid13))
call check(nf90_def_var(ncid, 'h_CAM', NF90_INT4, dim2, varid14))
call check(nf90_def_var(ncid, 'h_CAH', NF90_INT4, dim2, varid15))
call check(nf90_def_var(ncid, 'h_CA', NF90_INT4, dim2, varid16))
call check(nf90_def_var(ncid, 'h_CZ', NF90_INT4, dim3, varid17))


!!$call check(nf90_def_var(ncid, 'f_CTL', NF90_INT4, dim, varid18))
!!$call check(nf90_def_var(ncid, 'f_CTM', NF90_INT4, dim, varid19))
!!$call check(nf90_def_var(ncid, 'f_CTH', NF90_INT4, dim, varid20))
!!$call check(nf90_def_var(ncid, 'f_CT', NF90_INT4, dim, varid21))

!!$call check(nf90_def_var(ncid, 'h_CTL', NF90_INT4, dim4, varid23))
!!$call check(nf90_def_var(ncid, 'h_CTM', NF90_INT4, dim4, varid24))
!!$call check(nf90_def_var(ncid, 'h_CTH', NF90_INT4, dim4, varid25))
!!$call check(nf90_def_var(ncid, 'h_CT', NF90_INT4, dim4, varid26))
!!$
!!$call check(nf90_def_var(ncid, 'a_CTL', NF90_FLOAT, dim, varid27))
!!$call check(nf90_def_var(ncid, 'a_CTM', NF90_FLOAT, dim, varid28))
!!$call check(nf90_def_var(ncid, 'a_CTH', NF90_FLOAT, dim, varid29))
!!$call check(nf90_def_var(ncid, 'a_CT', NF90_FLOAT, dim, varid30))


    call check(nf90_enddef(ncid))
print *, 'end def'

    call check(nf90_put_var(ncid, varid1, low))
    call check(nf90_put_var(ncid, varid2, mid))
    call check(nf90_put_var(ncid, varid3, high))
print *, 'file gewex1'

    call check(nf90_put_var(ncid, varid4, colcloud))
    call check(nf90_put_var(ncid, varid6, height)) 
    call check(nf90_put_var(ncid, varid7, indtot)) 
    call check(nf90_put_var(ncid, varid8, f_CAL)) 
print *, 'file gewex1'

 call check(nf90_put_var(ncid, varid9, f_CAM)) 
 call check(nf90_put_var(ncid, varid10, f_CAH)) 
 call check(nf90_put_var(ncid, varid11, f_CA)) 
 call check(nf90_put_var(ncid, varid12, f_CZ))
print *, 'file gewex1'

 call check(nf90_put_var(ncid, varid13, hlow ))
 call check(nf90_put_var(ncid, varid14, hmid )) 
 call check(nf90_put_var(ncid, varid15, hhigh)) 
 call check(nf90_put_var(ncid, varid16, hcol)) 
print *, 'file gewex1'

 call check(nf90_put_var(ncid, varid17, hheight))
 
print *, 'file gewex1'


!!$ call check(nf90_put_var(ncid, varid18, f_CTL)) 
!!$ call check(nf90_put_var(ncid, varid19, f_CTM)) 
!!$ call check(nf90_put_var(ncid, varid20, f_CTH)) 
!!$ call check(nf90_put_var(ncid, varid21, f_CT)) 

! call check(nf90_put_var(ncid, varid23, hlowtemp ))
! call check(nf90_put_var(ncid, varid24, hmidtemp )) 
! call check(nf90_put_var(ncid, varid25, hhightemp)) 
! call check(nf90_put_var(ncid, varid26, hcoltemp)) 

!    call check(nf90_put_var(ncid, varid27, lowtemp))
!    call check(nf90_put_var(ncid, varid28, midtemp))
 !   call check(nf90_put_var(ncid, varid29, hightemp))
 !   call check(nf90_put_var(ncid, varid30, coltemp))



    call check(nf90_close(ncid))

endsubroutine map_recvar2nc5
!----------------------------------------------------------------------------!












!----------------------------------------------------------------------------!
! *** CREATE_MAPNC ***  This routine create a netcdf diag file, and its      !
!                       dimensions                                           !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf diagSR file                             !
! dname      : period of diagSR file (description of ncdf file)              !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! vlon       : values of model longitude from 180 to -180                    !
! vlat       : values of model latitude from 90 to -90                       !
! vprestop   : values of model altitude                                      !
! vsrmod     : values of model diag boxes                                    !
! alt        : values of model altitude                                      !
! dim        : dimension id of the diagSR var recorded in the ncdf           !
!              files                                                         !
! vtime      : number of days since 2000/01/01 for the trimonthly perdiod,   !
!              in day.                                                       !
!----------------------------------------------------------------------------!
! date       : date from the real-system clock and has form yyyymmdd.        !
! time       : time from the real-system clock and has form hhmmss.sss.      !
! zone       : represente the difference with respect to Coordinated         !
!              Universal Time (UTC), and has form (+-)hhmm.                  !
! value      : 8 dimension value which contains the year,month,day,hour,     !
!              minute, seconds and milliseconds of the real-time.            !
! diagmax    : values of SR diagbox                                          !
! ndims      : dimension of dim (=5 : lon,lat,alt,diag,time)                 !
! nc_id      : netcdf file id                                                !
! lon_varid  : variable id of longitude                                      !
! lat_varid  : variable id of latitude                                       !
! pres_varid : variable id of altitude                                       !
! srmod_varid: variable id of diagsr                                         !
! time_varid : variable id of time                                           !
! lon_dimid  : dimension id of longitude                                     !
! lat_dimid  : dimension id of latitude                                      !
! pres_dimid : dimension id of altitude                                      !
! srmod_dimid: dimension id of srbox                                         !
! time_dimid : dimension id of time                                          !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : create_diagnc(file8,file9,lonmod,latmod,altmod,srmod,resd,dimidsd,    !
!                    altmax,lonmax,latmax)                                   !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine create_diagnc(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vsrmod,vtime,dim,dim2,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon , nlat 
    integer, parameter :: diagmax = 19, ndims = 5, nv=2
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims),dim2(ndims),alt,dim3(nv),dim4(nv),dim5(nv)
    integer  ::  lon_varid,lat_varid,alt_varid, srmod_varid,time_varid,srmod_varid2,srmod_varid3,srmod_varid4,nc_id,alt_varid2,srmod_varid33,srmod_varid44
    integer  ::  lon_dimid,lat_dimid,alt_dimid,srmod_dimid,time_dimid,srmod_dimid2,srmod_dimid3,srmod_dimid4, nv_dimid,sr_dimid,sr_dimid2
    real  ::  vtime  
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(diagmax)  ::  vsrmod
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound
    real, dimension(diagmax-1)  ::  vsrmod2
    real, dimension(diagmax-1,2)  ::  vsrmod_bound


vsrmod2(:)=0;
vsrmod_bound(:,:)=0;

   
vsrmod2(1)=-888
vsrmod2(2)=-777

      do iz=3,diagmax-1
            vsrmod2(iz) = (vsrmod(iz)+vsrmod(iz+1))/2
       enddo

vsrmod_bound(1,:)=-888
vsrmod_bound(2,:)=-777
vsrmod_bound(3,1)=-776
vsrmod_bound(3,2)=0

      do iz=4,diagmax-1
         vsrmod_bound(iz,1)=vsrmod(iz)
         vsrmod_bound(iz,2)=vsrmod(iz+1)
      enddo



    call date_and_time(date,time,zone,value)

    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Histogram_of_Scattering_Ratio'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    call check(nf90_def_dim(nc_id, 'box', diagmax-4, srmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'box2', 3, srmod_dimid4 ))   

    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
    dim = (/lon_dimid, lat_dimid, alt_dimid, srmod_dimid2, time_dimid/)
    dim2 = (/lon_dimid, lat_dimid, alt_dimid, srmod_dimid4, time_dimid/)
    dim3 = (/alt_dimid, nv_dimid/)
    dim4 = (/srmod_dimid2, nv_dimid/)
    dim5 = (/srmod_dimid4, nv_dimid/)



    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim3, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'srbox_mid', NF90_FLOAT, srmod_dimid2,         &
               srmod_varid2))
    call check(nf90_put_att(nc_id, srmod_varid2, 'lon_name',                  &
               'Middle Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, srmod_varid2, 'axis','I')) 
  
    call check(nf90_def_var(nc_id, 'srbox_mid2', NF90_FLOAT, srmod_dimid4,         &
               srmod_varid4))
    call check(nf90_put_att(nc_id, srmod_varid4, 'lon_name',                  &
               'Middle Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid4, 'units','Arbitrary unit'))
 
    call check(nf90_def_var(nc_id, 'srbox_bound', NF90_FLOAT, dim4,         &
               srmod_varid33))
    call check(nf90_put_att(nc_id, srmod_varid33, 'lon_name',                  &
               'Boundarie Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid33, 'units','Arbitrary unit'))

    call check(nf90_def_var(nc_id, 'srbox_bound2', NF90_FLOAT, dim5,         &
               srmod_varid44))
    call check(nf90_put_att(nc_id, srmod_varid44, 'lon_name',                  &
               'Boundarie Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid44, 'units','Arbitrary unit'))
 

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  



    call check(nf90_enddef(nc_id))
	
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))

    call check(nf90_put_var(nc_id, srmod_varid2, vsrmod2(4:18)))
    call check(nf90_put_var(nc_id, srmod_varid4, vsrmod2(1:3)))
    call check(nf90_put_var(nc_id, srmod_varid33, vsrmod_bound(4:18,:)))
    call check(nf90_put_var(nc_id, srmod_varid44, vsrmod_bound(1:3,:)))



    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_diagnc
!----------------------------------------------------------------------------!

subroutine create_diagncpha(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vsrmod,vtime,dim,dim2,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon , nlat 
    integer, parameter :: diagmax = 19, ndims = 5, nv=2
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims),dim2(ndims),alt,dim3(nv),dim4(nv),dim5(nv)
    integer  ::  lon_varid,lat_varid,alt_varid, srmod_varid,time_varid,srmod_varid2,srmod_varid3,srmod_varid4,nc_id,alt_varid2,srmod_varid33,srmod_varid44
    integer  ::  lon_dimid,lat_dimid,alt_dimid,srmod_dimid,time_dimid,srmod_dimid2,srmod_dimid3,srmod_dimid4, nv_dimid,sr_dimid,sr_dimid2
    real  ::  vtime  
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(diagmax)  ::  vsrmod
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound
    real, dimension(diagmax-8)  ::  vsrmod2
    real, dimension(diagmax-8,2)  ::  vsrmod_bound


vsrmod2(:)=0;
vsrmod_bound(:,:)=0;

      do iz=8,diagmax-1
            vsrmod2(iz-7) = (vsrmod(iz)+vsrmod(iz+1))/2
       enddo

      do iz=8,diagmax-1
         vsrmod_bound(iz-7,1)=vsrmod(iz)
         vsrmod_bound(iz-7,2)=vsrmod(iz+1)
      enddo

print *, "srmod creer"

    call date_and_time(date,time,zone,value)

    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Phase_Histogram_of_Scattering_Ratio'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    call check(nf90_def_dim(nc_id, 'box', diagmax-8, srmod_dimid2 ))   

    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
    dim = (/lon_dimid, lat_dimid, alt_dimid, srmod_dimid2, time_dimid/)
    dim3 = (/alt_dimid, nv_dimid/)
    dim4 = (/srmod_dimid2, nv_dimid/)



    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim3, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'srbox_mid', NF90_FLOAT, srmod_dimid2,         &
               srmod_varid2))
    call check(nf90_put_att(nc_id, srmod_varid2, 'lon_name',                  &
               'Middle Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, srmod_varid2, 'axis','I')) 
  
    call check(nf90_def_var(nc_id, 'srbox_bound', NF90_FLOAT, dim4,         &
               srmod_varid33))
    call check(nf90_put_att(nc_id, srmod_varid33, 'lon_name',                  &
               'Boundarie Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid33, 'units','Arbitrary unit'))

 

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  



    call check(nf90_enddef(nc_id))
	
print *, 'titi'
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))

 print *, 'titi'
   call check(nf90_put_var(nc_id, srmod_varid2, vsrmod2))
    call check(nf90_put_var(nc_id, srmod_varid33, vsrmod_bound))

print *, 'titi'

    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_diagncpha
!----------------------------------------------------------------------------!

subroutine create_diagnc2(fname,dname,vlon,vlat,vprestop,vsrmod,vtime,dim2,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon , nlat 
    integer, parameter :: diagmax = 19, ndims = 4
    integer,dimension(8)  ::  value
    integer  ::  dim2(ndims),alt
    integer  ::  lon_varid,lat_varid,pres_varid, srmod_varid,time_varid,nc_id
    integer  ::  lon_dimid,lat_dimid,pres_dimid,srmod_dimid,time_dimid,srmod_dimid2
    real  ::  vtime  
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(diagmax)  ::  vsrmod
    real,dimension(alt)  ::  vprestop
    real, dimension(diagmax-1)  ::  vsrmod2


vsrmod2=vsrmod(2:19)
    call date_and_time(date,time,zone,value)

    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Histogram_of_Scattering_Ratio'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, pres_dimid))
    call check(nf90_def_dim(nc_id, 'box2', diagmax-1, srmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'box', diagmax, srmod_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
   ! dim = (/lon_dimid, lat_dimid, pres_dimid, srmod_dimid, time_dimid/)
    dim2 = (/lon_dimid, lat_dimid, pres_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'altitude', NF90_FLOAT, pres_dimid, pres_varid))
    call check(nf90_put_att(nc_id, pres_varid, 'lon_name','Altitude'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

    call check(nf90_def_var(nc_id, 'srbox', NF90_FLOAT, srmod_dimid,         &
               srmod_varid))
    call check(nf90_put_att(nc_id, srmod_varid, 'lon_name',                  &
               'Scattering Ratio Value'))
    call check(nf90_put_att(nc_id, srmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, srmod_varid, 'axis','I'))  

   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
    call check(nf90_enddef(nc_id))
	
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, pres_varid, vprestop))
    call check(nf90_put_var(nc_id, srmod_varid, vsrmod))
    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_diagnc2
!----------------------------------------------------------------------------!


subroutine create_depolnc2(fname,dname,vlon,vlat,vprestop,vsrmod,vdepolmod,vtime,dim2,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon , nlat 
    integer, parameter :: diagmax = 17, ndims = 5, depolmax = 21
    integer,dimension(8)  ::  value
    integer  ::  dim2(ndims),alt
    integer  ::  lon_varid,lat_varid,pres_varid, srmod_varid,time_varid,nc_id, depolmod_varid
    integer  ::  lon_dimid,lat_dimid,pres_dimid,srmod_dimid,time_dimid,srmod_dimid2, depolmod_dimid, depolmod_dimid2
    real  ::  vtime  
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(diagmax)  ::  vsrmod
    real, dimension(depolmax)  ::  vdepolmod
    real,dimension(alt)  ::  vprestop
    real, dimension(diagmax-1)  ::  vsrmod2
    real, dimension(depolmax-1)  ::  vdepolmod2


vsrmod2=vsrmod(2:17)
vdepolmod2=vdepolmod(2:21)

    call date_and_time(date,time,zone,value)

    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Histogram_of_Depolarization_Ratio'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, pres_dimid))
    call check(nf90_def_dim(nc_id, 'srbox2', diagmax-1, srmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'srbox', diagmax, srmod_dimid))
    call check(nf90_def_dim(nc_id, 'depolbox2', depolmax-1, depolmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'depolbox', depolmax, depolmod_dimid))


    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
   ! dim = (/lon_dimid, lat_dimid, pres_dimid, srmod_dimid, time_dimid/)
    dim2 = (/lon_dimid, lat_dimid, pres_dimid, srmod_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'altitude', NF90_FLOAT, pres_dimid, pres_varid))
    call check(nf90_put_att(nc_id, pres_varid, 'lon_name','Altitude'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

    call check(nf90_def_var(nc_id, 'srbox', NF90_FLOAT, srmod_dimid,         &
               srmod_varid))
    call check(nf90_put_att(nc_id, srmod_varid, 'lon_name',                  &
               'Scattering Ratio Value'))
    call check(nf90_put_att(nc_id, srmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, srmod_varid, 'axis','I'))  

    call check(nf90_def_var(nc_id, 'depolbox', NF90_FLOAT, depolmod_dimid,         &
               depolmod_varid))
    call check(nf90_put_att(nc_id, depolmod_varid, 'lon_name',                  &
               'Depolarization Ratio Value'))
    call check(nf90_put_att(nc_id, depolmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, depolmod_varid, 'axis','I'))  


   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
    call check(nf90_enddef(nc_id))
	
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, pres_varid, vprestop))
    call check(nf90_put_var(nc_id, srmod_varid, vsrmod))
    call check(nf90_put_var(nc_id, depolmod_varid, vdepolmod))

    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_depolnc2
!----------------------------------------------------------------------------!


subroutine create_depolnc(fname,dname,vlon,vlat,vprestop,vsrmod,vdepolmod,vtime,dim,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon , nlat 
    integer, parameter :: diagmax = 17, ndims = 6, depolmax = 21
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims),dim2(ndims-1),alt
    integer  ::  lon_varid,lat_varid,pres_varid, srmod_varid,time_varid,nc_id, depolmod_varid
    integer  ::  lon_dimid,lat_dimid,pres_dimid,srmod_dimid,time_dimid,srmod_dimid2, depolmod_dimid2, depolmod_dimid
    real  ::  vtime  
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(diagmax)  ::  vsrmod
    real, dimension(depolmax)  ::  vdepolmod
    real,dimension(alt)  ::  vprestop
    real, dimension(diagmax-1)  ::  vsrmod2
    real, dimension(depolmax-1)  ::  vdepolmod2



vsrmod2=vsrmod(2:17)
vdepolmod2=vdepolmod(2:21)

    call date_and_time(date,time,zone,value)

    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Histogram_of_Depolarization_Ratio'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, pres_dimid))
    call check(nf90_def_dim(nc_id, 'srbox2', diagmax-1, srmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'srbox', diagmax, srmod_dimid))
    call check(nf90_def_dim(nc_id, 'depolbox2', depolmax-1, depolmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'depolbox', depolmax, depolmod_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
    dim = (/lon_dimid, lat_dimid, pres_dimid, srmod_dimid2, depolmod_dimid2, time_dimid/)
    
   ! dim2 = (/lon_dimid, lat_dimid, pres_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'altitude', NF90_FLOAT, pres_dimid, pres_varid))
    call check(nf90_put_att(nc_id, pres_varid, 'lon_name','Altitude'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

    call check(nf90_def_var(nc_id, 'srbox', NF90_FLOAT, srmod_dimid,         &
               srmod_varid))
    call check(nf90_put_att(nc_id, srmod_varid, 'lon_name',                  &
               'Scattering Ratio Value'))
    call check(nf90_put_att(nc_id, srmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, srmod_varid, 'axis','I'))  

    call check(nf90_def_var(nc_id, 'depolbox', NF90_FLOAT, depolmod_dimid,         &
               depolmod_varid))
    call check(nf90_put_att(nc_id, depolmod_varid, 'lon_name',                  &
               'Depolarization Ratio Value'))
    call check(nf90_put_att(nc_id, depolmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, depolmod_varid, 'axis','I'))  

   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
    call check(nf90_enddef(nc_id))
	
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, pres_varid, vprestop))
    call check(nf90_put_var(nc_id, srmod_varid, vsrmod))
    call check(nf90_put_var(nc_id, depolmod_varid, vdepolmod))

    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_depolnc
!----------------------------------------------------------------------------!

subroutine create_diagPHAnc(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vtempmod,vtime,dim,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon , nlat, iz 
    integer, parameter :: tempmax = 35, ndims = 6, nv=2, ncat=3
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims),alt,dim3(nv),dim4(nv),dim5(nv)
    integer  ::  lon_varid,lat_varid,alt_varid, tempmod_varid,time_varid,tempmod_varid2,srmod_varid3,srmod_varid4,nc_id,alt_varid2,srmod_varid33,srmod_varid44,cat_varid
    integer  ::  lon_dimid,lat_dimid,alt_dimid,tempmod_dimid,time_dimid,tempmod_dimid2,srmod_dimid3,srmod_dimid4, nv_dimid,sr_dimid,sr_dimid2, cat_dimid
    real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(tempmax-1)  ::  vtempmid
    real, dimension(tempmax)  ::  vtempmod
    real, dimension(tempmax-1,2)  ::  vtempmod2
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound
    character(len=3),dimension(ncat,ncat) :: vcat

vcat(1,:)='LIQ'
vcat(2,:)='ICE'
vcat(3,:)='VAP'



      do iz=1,tempmax-1
         vtempmod2(iz,1)=vtempmod(iz)
         vtempmod2(iz,2)=vtempmod(iz+1)
         vtempmid(iz)=(vtempmod(iz)+vtempmod(iz+1))/2
      enddo



    call date_and_time(date,time,zone,value)

    call check(nf90_create('/homedata/mreverdy/GOCCP/'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Histogram_of_Phase_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'tempbox', tempmax-1, tempmod_dimid )) 
    call check(nf90_def_dim(nc_id, 'cat', ncat, cat_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
    dim = (/lon_dimid, lat_dimid, alt_dimid, tempmod_dimid, cat_dimid, time_dimid/)
    dim3 = (/alt_dimid, nv_dimid/)
    dim4 = (/tempmod_dimid, nv_dimid/)
    dim5 = (/cat_dimid, cat_dimid /)


    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim3, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'temp_mid', NF90_FLOAT, tempmod_dimid,         &
               tempmod_varid))
    call check(nf90_put_att(nc_id, tempmod_varid, 'lon_name',                  &
               'Middle Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, tempmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, tempmod_varid, 'axis','I')) 

    call check(nf90_def_var(nc_id, 'temp_bound', NF90_FLOAT, dim4,         &
               tempmod_varid2))
    call check(nf90_put_att(nc_id, tempmod_varid2, 'lon_name',                  &
               'Middle Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, tempmod_varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, tempmod_varid2, 'axis','I')) 
  
    call check(nf90_def_var(nc_id, 'category', NF90_CHAR, dim5, cat_varid))
    call check(nf90_put_att(nc_id, cat_varid, 'lon_name','Category'))
    call check(nf90_put_att(nc_id, cat_varid, 'units','Arbitrary unit'))


    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  



    call check(nf90_enddef(nc_id))
	
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))

    call check(nf90_put_var(nc_id, tempmod_varid, vtempmid))
    call check(nf90_put_var(nc_id, tempmod_varid2, vtempmod2))

    call check(nf90_put_var(nc_id, cat_varid, vcat))

    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_diagPHAnc
!----------------------------------------------------------------------------!

subroutine diagPHA_recvar2nc3(diagpha,dim,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999
    integer, parameter :: tempmax = 34, ndims = 6
    integer  ::  dim(ndims),alt
    integer  ::  varid, varid2, ncid
    character(LEN=*)   ::  fname
    real,dimension(nlon,nlat,alt,tempmax,3)  ::  diagpha


    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarPHASE', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar ICE/LIQUID/VAPOR Phase occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))


    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diagpha(:,:,:,:,:)))

    call check(nf90_close(ncid))
    
end subroutine diagPHA_recvar2nc3

!----------------------------------------------------------------------------!
! *** DIAG_RECVAR2NC *** This routine record the diag variables in the netcdf!
!                        prof file                                           !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf diagSR file                             !
! diag       : diagSR variable                                               !
! dim        : dimension id of the isccp variables                           !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! alt        : number of model altitude                                      !
!----------------------------------------------------------------------------!
! nan        : Not a Number value                                            !
! diagmax    : number of SR diagbox                                          !
! ndims      : dimension of dim (= 5 : lon,lat,alt,diagbox,time)             !
! ncid       : netcdf file id                                                !
! varid      : variable id of diagSR                                         !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : diag_recvar2nc(monthdiagSR,dimidsd,file8,altmax,lonmax,latmax)        !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine diag_recvar2nc(diag,diag1,dim,dim2,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 15, ndims = 5
    integer  ::  dim(ndims), dim2(ndims-1),alt
    integer  ::  varid,varid2, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt,diagmax)  ::  diag
    real*4,dimension(nlon,nlat,alt)  ::  diag1


    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarsr532_Occ', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Number of Scattering Ratio occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'negdiagSR_Occ', NF90_FLOAT, dim2, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                         &
             'Number of Scattering Ratio occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag))
    call check(nf90_put_var(ncid, varid2, diag1))

    call check(nf90_close(ncid))
    
end subroutine diag_recvar2nc
!----------------------------------------------------------------------------!



subroutine diag_recvar2nc2(diag,dim,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 16, ndims = 5
    integer  ::  dim(ndims), alt
    integer  ::  varid, ncid
    character(LEN=*)   ::  fname
    real*8,dimension(nlon,nlat,alt,diagmax)  ::  diag
 !   real*4,dimension(nlon,nlat,alt)  ::  diag1

    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarsr532_Frac', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
           'Lidar Scattering Ratio CFAD (532nm) Fraction of occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag))

    call check(nf90_close(ncid))
    
end subroutine diag_recvar2nc2
!----------------------------------------------------------------------------!

subroutine diag_recvar2nc3(diag,dim,dim2,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 18, ndims = 5
    integer  ::  dim(ndims),dim2(ndims),alt
    integer  ::  varid, varid2, ncid
    character(LEN=*)   ::  fname
    real,dimension(nlon,nlat,alt,diagmax)  ::  diag


    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarsr532_Occ', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))

   call check(nf90_def_var(ncid, 'cfad_lidarsr532_Occ2', NF90_FLOAT, dim2, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))



    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag(:,:,:,4:18)))
    call check(nf90_put_var(ncid, varid2, diag(:,:,:,1:3)))

    call check(nf90_close(ncid))
    
end subroutine diag_recvar2nc3
!----------------------------------------------------------------------------!

subroutine diag_recvar2nc3pha(diag,dim,dim2,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999.
    integer, parameter :: diagmax = 11, ndims = 5
    integer  ::  dim(ndims),dim2(ndims),alt
    integer  ::  varid, varid2, varid3, ncid
    character(LEN=*)   ::  fname
    real,dimension(nlon,nlat,alt,diagmax,3)  ::  diag


    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarsr532_liq', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))

   call check(nf90_def_var(ncid, 'cfad_lidarsr532_ice', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

   call check(nf90_def_var(ncid, 'cfad_lidarsr532_un', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid3, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))



    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag(:,:,:,:,1)))
    call check(nf90_put_var(ncid, varid2, diag(:,:,:,:,2)))
    call check(nf90_put_var(ncid, varid3, diag(:,:,:,:,3)))

    call check(nf90_close(ncid))
    
end subroutine diag_recvar2nc3pha
!----------------------------------------------------------------------------!

subroutine depol_recvar2nc3(diag,dim,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 17, ndims = 6, depolmax = 21
    integer  ::  dim(ndims),alt
    integer  ::  varid, ncid
    character(LEN=*)   ::  fname
    real*8,dimension(nlon,nlat,alt,diagmax-1,depolmax-1)  ::  diag


    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidardepol532_Occ', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar Depolarization Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag))


    call check(nf90_close(ncid))
    
end subroutine depol_recvar2nc3
!----------------------------------------------------------------------------!


subroutine diag_recvar2nc4(diag,dim,fname,alt,nlon,nlat,ndiag)

    use netcdf
    implicit none

    integer :: nlon , nlat
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 18, ndims = 4
    integer  ::  dim(ndims),alt
    integer  ::  varid, ncid
    character(LEN=*)   ::  fname
    character  ::  ndiag*2
    integer,dimension(nlon,nlat,alt)  ::  diag

print *, 'monthdiagSR_Occ'// trim(ADJUSTL(ndiag))
    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarsr532_Occ'//trim(ADJUSTL(ndiag)), NF90_INT4, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag))


    call check(nf90_close(ncid))
    
end subroutine diag_recvar2nc4
!----------------------------------------------------------------------------!


subroutine depol_recvar2nc4(diag,dim,fname,alt,nlon,nlat,ndiag)

    use netcdf
    implicit none

    integer :: nlon , nlat
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 16, ndims = 5
    integer  ::  dim(ndims),alt
    integer  ::  varid, ncid
    character(LEN=*)   ::  fname
    character  ::  ndiag*2
    real*8,dimension(nlon,nlat,alt,diagmax)  ::  diag

print *, 'monthdepolSR_Occ'// trim(ADJUSTL(ndiag))
    call check(nf90_open('/homedata/mreverdy/GOCCP/'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidardepol532_Occ'//trim(ADJUSTL(ndiag)), NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar Depolarization Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag))


    call check(nf90_close(ncid))
    
end subroutine depol_recvar2nc4
!----------------------------------------------------------------------------!







subroutine rdnc3(fname,var,alt,nlon,nlat,nvar)
  use netcdf
  implicit none
character(len=*)  ::  nvar, fname
integer  ::  ncid,varid
integer :: nlon , nlat, alt
integer,dimension(nlon,nlat,alt) :: var

print *, nvar

call check(NF90_OPEN(fname,NF90_NOWRITE,ncid))
call check(NF90_inq_varid(ncid,nvar,varid))
call check(NF90_get_var(ncid,varid,var))
call check(NF90_CLOSE(ncid))
 !  call rdnc3('/homedata/mreverdy/GOCCP/'//trim(file7)//'.tmp',monthdiagSR1,'monthdiagSR_Occ'//trim(adjustl(idiagc)))

end subroutine rdnc3

subroutine rdnc4(fname,var,alt,nlon,nlat,nvar)
  use netcdf
  implicit none
character(len=*)  ::  nvar, fname
integer  ::  ncid,varid
integer, parameter  :: diagmax = 17
integer :: nlon , nlat, alt
real*8,dimension(nlon,nlat,alt,diagmax-1) :: var

print *, nvar

call check(NF90_OPEN(fname,NF90_NOWRITE,ncid))
call check(NF90_inq_varid(ncid,nvar,varid))
call check(NF90_get_var(ncid,varid,var))
call check(NF90_CLOSE(ncid))
 !  call rdnc3('/homedata/mreverdy/GOCCP/'//trim(file7)//'.tmp',monthdiagSR1,'monthdiagSR_Occ'//trim(adjustl(idiagc)))

end subroutine rdnc4



!****************************************************************************!


end program A-GCM

!****************************************************************************!
!                                                                            !
!                               END PROGRAM                                  !
!                                                                            !
!****************************************************************************!



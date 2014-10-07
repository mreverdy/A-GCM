!****************************************************************************!
!                                                                            !
!                            GOECP 1.00 			             !
!                            						     !
!									     !
!                    Mathieu Reverdy / LMD / ESA                             !
!                      last update 11/19/2013                                !
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

!************************* SUBROUTINE SCHEME ***************************!
!                                                                       !
! 									!
!                                                                       !
!-----------------------------------------------------------------------!

program AGCM

	use netcdf
	implicit none
  
	integer i
	character*100 paramfile, filename, rep, filename2,filename3,filename4,filename5,filename6,filename7,filename8,filename9,filename10,filename11,filename12,filename13,filename14,paramdaynight,filename16, filename17
	real, dimension(:,:), allocatable ::rayParaValues, mieParaValues, totPerpValues, ATBatlid, ATBatlid20, ATBatlid480,SRvert, SRthres, SRfiltered,ATBmolatlid20, ATBmol480, ATBvert, ATBmolvert, ATBmolatlid, SRmoy, filtreelevation,ATBatlid20elevation,ATBmolatlid20elevation,ATBatlidmoy,ATBmolatlidmoy,ATBmolatlidmoyelevation,ATBatlidmoyelevation, p_mol_355, pmol20,pmolatlidmoy
	real, dimension(:), allocatable:: altValues, altValues20, lidAltValues, clidValues,lonValues, rejectionflag, elevation,altValuesaverage

	integer numLon,numAlt,numzscene 
	real, dimension(:,:,:), allocatable :: pressureValues, temperatureValues
	real, dimension(:), allocatable :: zsceneValues,altaverageValues,z480Values
	real, dimension(:), allocatable :: pressureValues2,temperatureValues2
	real vertTemp

!Lecture des fichiers
	call paramread(paramfile,rep,filename,filename2,filename3,filename4,filename5,filename6,filename7,filename8,filename9,filename10,filename11,filename12,filename13,filename14,paramdaynight,filename16,filename17)

	call lidarread(rep,filename,rayParaValues,mieParaValues,totPerpValues,clidValues,lidAltValues,altValues,numLon,numAlt,lonValues)

	call environmentread(rep,filename2,pressureValues,temperatureValues,numzscene,zsceneValues)

!Calcul de l'ATBmol et de l'ATB. VR: 100 m méthode 1

	call ATBcalc(rayParaValues,mieParaValues,totPerpValues,clidValues,lidAltValues,altValues,ATBatlid,numLon,numAlt)
	call ATBmolcalc(rayParaValues,clidValues,lidAltValues,altValues,ATBmolatlid,numLon,numAlt)

!Calcul de l'ATBmol : 100 m méthode 2
	call ATBmolcalc2(pressureValues,temperatureValues,numzscene,zsceneValues,p_mol_355)

!Ecriture des résultats
	call lidarwrite2D(rep,filename6,ATBmolatlid,lonValues,altValues)
	call lidarwrite2D(rep,filename5,ATBatlid,lonValues,altValues)

! Choix du moyennage vertical
	call verticaldimension(altValues,altValues20,altValuesaverage,vertTemp)

!interpolation sur 20m
	call interp2D(ATBatlid,ATBatlid20,altValues,altValues20,vertTemp)
	call interp2D(ATBmolatlid,ATBmolatlid20,altValues,altValues20,vertTemp)
	call interp2D(p_mol_355,pmol20,zsceneValues,altValues20,vertTemp)	

!Moyennage sur xm (240,480,...)
	call average2D(ATBatlid20,vertTemp,altValues20,ATBatlidmoy)
	call average2D(ATBmolatlid20,vertTemp,altValues20,ATBmolatlidmoy)	
	call average2D(pmol20,vertTemp,altValues20,pmolatlidmoy)		

	call lidarwrite2D(rep,filename3,ATBatlidmoy,lonValues,altValuesaverage)
	call lidarwrite2D(rep,filename4,ATBmolatlidmoy,lonValues,altValuesaverage)
	call lidarwrite2D(rep,filename17,pmolatlidmoy,lonValues,altValuesaverage)

!Rejection des profils trop bruités
	call rejectionfiltering(ATBatlidmoy,ATBmolatlidmoy,altValuesaverage,ATBvert,ATBmolvert)
	call SRcalc(ATBmolvert,ATBvert,SRvert)

!Ecriture des résultats
        call lidarwrite2D(rep,filename10,ATBvert,lonValues,altValuesaverage)
        call lidarwrite2D(rep,filename11,ATBmolvert,lonValues,altValuesaverage)
	call lidarwrite2D(rep,filename7,SRvert,lonValues,altValuesaverage)

!Détection de la surface. VR: 100 m 
	call surfaceelevation(altValues,ATBatlid,filtreelevation,elevation)
	call interp2Delevation(filtreelevation,ATBatlid,altValues,altValues20,ATBatlid20elevation)
	call interp2Delevation(filtreelevation,ATBmolatlid,altValues,altValues20,ATBmolatlid20elevation)

	call lidarwrite2D(rep,filename12,ATBatlid20elevation,lonValues,altValues20)

	call average2Delevation(ATBatlid20elevation,vertTemp,altValues20,altValuesaverage,elevation,ATBatlidmoyelevation)
	call average2Delevation(ATBmolatlid20elevation,vertTemp,altValues20,altValuesaverage,elevation,ATBmolatlidmoyelevation)		

	call lidarwrite2D(rep,filename13,ATBatlidmoyelevation,lonValues,altValuesaverage)
	call lidarwrite2D(rep,filename14,ATBmolatlidmoyelevation,lonValues,altValuesaverage)

	call SRcalc(ATBmolatlidmoyelevation,ATBatlidmoyelevation,SRmoy)
	call SRthrescalc(paramdaynight,SRmoy,SRthres)

	call lidarwrite2D(rep,filename16,SRmoy,lonValues,altValuesaverage)
	call lidarwrite2D(rep,filename8,SRthres,lonValues,altValuesaverage)
!
!	call SNR(SRthres,SRfiltered)
!
!	call lidarwrite2D(rep,filename9,SRfiltered,lonValues,z480Values)




!	

!	call SRcalc(ATBmolvert,ATBvert,SRvert)

!	call lidarwrite1D(rep,filename4,z480Values)

	


!	




contains

!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                            END OF PROGRAM                                *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!



!****************************************************************************!
!****************************************************************************!
!**************************    SUBROUTINE LIST    ***************************!
!****************************************************************************!
!****************************************************************************!
!----------------------------------------------------------------------------!
subroutine paramread(paramfile,rep,filename,filename2,filename3,filename4,filename5,filename6,filename7,filename8,filename9,filename10,filename11,filename12,filename13,filename14,paramdaynight,filename16,filename17)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Read the parameters in the parameter file. 
!ATLID Level 1 directories, filenames. 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

	implicit none
	

	character*100 paramfile,rep,filename,filename2,filename3,filename4,filename5,filename6,filename7
	character*100 filename8,filename9,filename10,filename11,filename12,filename13,filename14,paramdaynight,filename16,filename17
	integer i,io
	real, dimension(:),allocatable :: var1

 	write(*, '(A)', advance="no") 'parameter file name? ='
	read(*,*) paramfile



open(unit=55,file=trim(paramfile)//'.txt')
print *,'Read parameter file'
   
     read(55,'(A)') rep 
     read(55,'(A)') filename
     read(55,'(A)') filename2
     read(55,'(A)') filename3
     read(55,'(A)') filename4
     read(55,'(A)') filename5
     read(55,'(A)') filename6
     read(55,'(A)') filename7
     read(55,'(A)') filename8
     read(55,'(A)') filename9
     read(55,'(A)') filename10
     read(55,'(A)') filename11
     read(55,'(A)') filename12
     read(55,'(A)') filename13
     read(55,'(A)') filename14
     read(55,'(A)') paramdaynight
     read(55,'(A)') filename16
     read(55,'(A)') filename17

end subroutine paramread
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine lidarread(rep,filename,rayParaValues,mieParaValues,totPerpValues,clidValues,lidAltValues,altValues,numLon,numAlt,lonValues)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Read the different data in the level 1 ATLID files. 
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!


	use netcdf
	implicit none
	
	integer, dimension(nf90_max_var_dims) :: dimIDs
	character*100 filename, repdata, rep
	integer lon, lat,alt, lidAlt, clid, status, ncid
	integer numLat, numLidAlt, numClid
	integer numLon, numAlt
	real, dimension(:), allocatable :: lonValues, latValues
	real, dimension(:), allocatable:: altValues, lidAltValues, clidValues    	
	integer rayPara, miePara, totPerp
	real, dimension(:,:), allocatable :: rayParaValues, mieParaValues, totPerpValues 

repdata=trim(rep)//filename

print *,'Read lidar netcdf files'
print *,trim(repdata)

status=nf90_open(repdata, nf90_nowrite, ncid)
status=nf90_inq_varid(ncid, "x_scene", lon)
status=nf90_inquire_variable(ncid, lon, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numLon)
allocate(lonValues(numLon))
status=nf90_get_var(ncid, lon, lonValues)

status=nf90_inq_varid(ncid, "y_scene", lat)
status=nf90_inquire_variable(ncid, lat, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numLat)
allocate(latValues(numLat))
status=nf90_get_var(ncid, lat, latValues)

status=nf90_inq_varid(ncid, "height", alt)
status=nf90_inquire_variable(ncid, alt, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numAlt)
allocate(altValues(numAlt))
status=nf90_get_var(ncid, alt, altValues)

status=nf90_inq_varid(ncid, "Lid_alt", lidAlt)
status=nf90_inquire_variable(ncid, lidAlt, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numLidAlt)
allocate(lidAltValues(numLidAlt))
status=nf90_get_var(ncid, lidAlt, lidAltValues)

status=nf90_inq_varid(ncid, "C_lid", clid)
status=nf90_inquire_variable(ncid, clid, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numClid)
allocate(clidValues(numClid))
status=nf90_get_var(ncid, clid, clidValues)

status=nf90_inq_varid(ncid, "Ray_para", rayPara)
status=nf90_inquire_variable(ncid, rayPara, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numLon)
status=nf90_inquire_dimension(ncid, dimIDs(2), len = numAlt)
allocate(rayParaValues(numLon,numAlt))
status=nf90_get_var(ncid, rayPara, rayParaValues)

status=nf90_inq_varid(ncid, "Mie_para", miePara)
status=nf90_inquire_variable(ncid, miePara, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numLon)
status=nf90_inquire_dimension(ncid, dimIDs(2), len = numAlt)
allocate(mieParaValues(numLon,numAlt))
status=nf90_get_var(ncid, miePara, mieParaValues)

status=nf90_inq_varid(ncid, "Tot_perp", totPerp)
status=nf90_inquire_variable(ncid, totPerp, dimids = dimIDs)
status=nf90_inquire_dimension(ncid, dimIDs(1), len = numLon)
status=nf90_inquire_dimension(ncid, dimIDs(2), len = numAlt)
allocate(totperpValues(numLon,numAlt))
status=nf90_get_var(ncid, totPerp, totPerpValues)

end subroutine lidarread
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine environmentread(rep,filename2,pressureValues,temperatureValues,numzscene,zsceneValues)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Read the environment data in the ATLID level 1 files 
!NB: may change to ECMWF files
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!


	use netcdf
	implicit none
	
 
	integer, dimension(nf90_max_var_dims) :: dimIDs
   	character*100 filename2, repdata, rep
	integer xscene, yscene,zscene, status, ncid2,pressure,temperature
	integer numyscene, numxscene, numzscene
	real, dimension(:), allocatable :: xsceneValues, ysceneValues, zsceneValues
	real, dimension(:,:,:), allocatable :: pressureValues, temperatureValues

repdata=trim(rep)//filename2

print *,'Read environment netcdf files'
print *,trim(repdata)

status=nf90_open(repdata, nf90_nowrite, ncid2)
status=nf90_inq_varid(ncid2, "x", xscene)
status=nf90_inquire_variable(ncid2, xscene, dimids = dimIDs)
status=nf90_inquire_dimension(ncid2, dimIDs(1), len = numxscene)
allocate(xsceneValues(numxscene))
status=nf90_get_var(ncid2, xscene, xsceneValues)

status=nf90_inq_varid(ncid2, "y", yscene)
status=nf90_inquire_variable(ncid2, yscene, dimids = dimIDs)
status=nf90_inquire_dimension(ncid2, dimIDs(1), len = numyscene)
allocate(ysceneValues(numyscene))
status=nf90_get_var(ncid2, yscene, ysceneValues)

status=nf90_inq_varid(ncid2, "z", zscene)
status=nf90_inquire_variable(ncid2, zscene, dimids = dimIDs)
status=nf90_inquire_dimension(ncid2, dimIDs(1), len = numzscene)
allocate(zsceneValues(numzscene))
status=nf90_get_var(ncid2, zscene, zsceneValues)

status=nf90_inq_varid(ncid2, "Pressure", pressure)
status=nf90_inquire_variable(ncid2, pressure, dimids = dimIDs)
status=nf90_inquire_dimension(ncid2, dimIDs(1), len = numxscene)
status=nf90_inquire_dimension(ncid2, dimIDs(2), len = numyscene)
status=nf90_inquire_dimension(ncid2, dimIDs(3), len = numzscene)
allocate(pressureValues(numxscene,numyscene,numzscene))
status=nf90_get_var(ncid2, pressure, pressureValues)

status=nf90_inq_varid(ncid2, "Temperature", temperature)
status=nf90_inquire_variable(ncid2, temperature, dimids = dimIDs)
status=nf90_inquire_dimension(ncid2, dimIDs(1), len = numxscene)
status=nf90_inquire_dimension(ncid2, dimIDs(2), len = numyscene)
status=nf90_inquire_dimension(ncid2, dimIDs(3), len = numzscene)
allocate(temperatureValues(numxscene,numyscene,numzscene))
status=nf90_get_var(ncid2, temperature, temperatureValues)

end subroutine environmentread
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine ATBcalc(rayParaValues,mieParaValues,totPerpValues,clidValues,lidAltValues,altValues,ATBatlid,numLon,numAlt)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Calculate the ATB at native vertical resolution: 100m 
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
	implicit none

	integer i,j,numLon,numAlt
	real, dimension(:), allocatable:: altValues, lidAltValues, clidValues,BB 
	real, dimension(:,:), allocatable :: rayParaValues, mieParaValues, totPerpValues  
	real, dimension(:,:), allocatable:: ATBatlid
	real, dimension(:,:), allocatable :: AA

print *,'Compute ATB'


allocate(AA(numLon,numAlt),BB(numAlt),ATBatlid(numLon,numAlt))

AA=(rayParaValues+mieParaValues+totPerpValues)/clidValues(1)

BB=((lidAltValues(1)-altValues)*1000)**2

do i = 1 , numLon
ATBatlid(i,:)=((AA(i,:)*BB))*1000
end do

print *,'ATB Computed'
end subroutine ATBcalc
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine ATBmolcalc(rayParaValues,clidValues,lidAltValues,altValues,ATBmolatlid,numLon,numAlt)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Calculate the ATBmol at native vertical resolution: 100m 
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
	implicit none

	integer i,j,numLon,numAlt
	real, dimension(:), allocatable:: altValues, lidAltValues, clidValues,BB 
	real, dimension(:,:), allocatable :: rayParaValues
	real, dimension(:,:), allocatable:: ATBmolatlid
	real, dimension(:,:), allocatable :: AA

print *,'Compute ATB'


allocate(AA(numLon,numAlt),BB(numAlt),ATBmolatlid(numLon,numAlt))

AA=(rayParaValues)/clidValues(1)

BB=((lidAltValues(1)-altValues)*1000)**2

do i = 1 , numLon
ATBmolatlid(i,:)=((AA(i,:)*BB))*1000
end do

print *,'ATB mol Computed'
end subroutine ATBmolcalc
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
subroutine ATBmolcalc2(var,var1,var2,var3,fvar)
!ATBmolcalc2(pressureValues,temperatureValues,numzscene,zsceneValues,p_mol_355)

	use netcdf
	implicit none

	integer i,j,var2,u
	real Cmol_355,kb
	real, dimension(:,:,:), allocatable :: var1, var
	real, dimension(:), allocatable :: beta_mol_355,alpha_mol_355,tau_mol_355,tau_mol_lay_355
	real, parameter :: pi = 3.141592653589793
	 
	real, dimension(:), allocatable :: var3
	real, dimension(:,:), allocatable :: Press,Temp,fvar

	integer x_dimid,ncid,status,varid1,varid2
	double precision A
	integer,dimension(1) :: dimids



print *,'Compute molecular ATB'


	
allocate(Press(size(var,1),var2),Temp(size(var,1),var2),beta_mol_355(var2),alpha_mol_355(var2),tau_mol_355(var2),tau_mol_lay_355(var2),fvar(size(var,1),var2))


Press=var1(:,101,:)
Temp=var(:,101,:)


	Cmol_355=3.2662e-31
	kb=1.38e-23

do u=1,size(var,1)

beta_mol_355(:)=Press(u,:)*100/kb/Temp(u,:)*Cmol_355
alpha_mol_355(:)=8*pi/3*beta_mol_355

do i=1,999

  tau_mol_355(i)=alpha_mol_355(i)*(zsceneValues(i+1)-zsceneValues(i))*1000
end do

tau_mol_355(1000)=alpha_mol_355(1000)*(99.9749-99.875)*1000


do i=999,1,-1
    tau_mol_355(i)=tau_mol_355(i)+tau_mol_355(i+1)
end do

A=-2.0*tau_mol_355(1000)
fvar(u,1000)=beta_mol_355(1000)/(2*tau_mol_355(1000))*(1.-exp(A))

do i=999,1,-1
    tau_mol_lay_355(i)=tau_mol_355(i)-tau_mol_355(i+1);

   if (tau_mol_lay_355(i) .gt. 0) then

      fvar(u,i)=beta_mol_355(i)*exp(-2*tau_mol_355(i+1))/(2*tau_mol_lay_355(i))*(1-exp(-2*tau_mol_lay_355(i)));
 else

        fvar(u,i)=beta_mol_355(i)*exp(-2*tau_mol_355(i+1));
    end if
end do
end do

print *,'molecular ATB Computed'
end subroutine ATBmolcalc2
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
subroutine SRcalc(var1,var2,fvar)
!SRcalc(ATBmolvert,ATBvert,SRvert)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Calculate the Scattering Ratio
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
	implicit none

	integer i,j
	real, dimension(:,:), allocatable:: var1,var2,fvar
	

print *,'Compute SR'

allocate(fvar(size(var2,1),size(var2,2)))

do i=1,size(var2,1)

	do j=1,size(var2,2)
	fvar(i,j)=var2(i,j)/var1(i,j)
	end do
end do

print *,'SR Computed'

end subroutine SRcalc
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine SRthrescalc(var,var2,fvar)
!SRthrescalc(paramdaynight,SR,SRthres)
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Apply a cloud detection threshold tto he Scattering Ratio
!NB: Modified SR threshold. SR>2.92 => clouds. May change
!NB: Day threshold may be modified. 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!


	implicit none

	integer i,j
	real, dimension(:,:), allocatable::var2,fvar
	character*100 var


print *,'Compute SR threshold'

allocate(fvar(size(var2,1),size(var2,2)))

if (trim(var) .eq. 'night' )then

do i=1,size(var2,1)

	do j=1,size(var2,2)
		if (var2(i,j) .lt. 2.92) then
		fvar(i,j)=-999
		else
		fvar(i,j)=var2(i,j)
		endif
	end do
end do

elseif (trim(var) .eq. 'day' )then

do i=1,size(var2,1)

	do j=1,size(var2,2)
		if (var2(i,j) .lt. 2.92) then
		fvar(i,j)=-999
		else
		fvar(i,j)=var2(i,j)
		endif
	end do
end do

endif

print *,'SR threshold Computed'

end subroutine SRthrescalc
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine SNR(SRthres,SRfiltered)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

	implicit none

	integer i,j
	real, dimension(:,:), allocatable::SRthres,SRfiltered


print *,'Compute SNR threshold'

allocate(SRfiltered(size(SRthres,1),size(SRthres,2)))

do i=1,size(SRthres,1)

	do j=1,size(SRthres,2)
		if (SRthres(i,j) .gt. 150) then
		SRfiltered(i,j)=-999
		else
		SRfiltered(i,j)=SRthres(i,j)
		endif
	end do
end do

print *,'SNR threshold Computed'

end subroutine SNR
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine verticaldimension(fvar,fvar2,fvar3,vertTemp)
!verticaldimension(altValues,altValues20,altValuesaverage,vertTemp)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Define vertical step to process with. 
!Define altitude array
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

	real vertTemp
	real, dimension(:), allocatable :: fvar,fvar2,fvar3
	integer a,b
	120 format (A,E)
 	write(*, '(A)', advance="no") 'vertical resolution (km)? ='
	read(*,*) vertTemp

	a=floor((fvar(size(fvar))+0.02)/0.02)
	b=floor((fvar(size(fvar))+0.1)/vertTemp)

allocate(fvar2(a))
allocate(fvar3(b))

fvar2(1)=0

do i=2,size(fvar2)
fvar2(i)=fvar2(i-1)+0.02
enddo

fvar3(1)=0

do i=2,size(fvar3)
fvar3(i)=fvar3(i-1)+vertTemp
enddo

end subroutine verticaldimension
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine interp2D(var,var2,fvar,fvar2,fvar3)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Vertical interpolation over 20 meters
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

  implicit none

	real, dimension(:,:), allocatable :: var, var2
	real, dimension(:), allocatable :: fvar, fvar2
	integer i,j,k
	real a,b,fvar3
	


print *,"data interpolation over 20 meters"
allocate(var2(size(var,1),size(fvar2)))


do k=1,size(var,1)
	do i=2,size(fvar)
	a=(var(k,i)-var(k,i-1))/(fvar(i)-fvar(i-1))
	b=var(k,i)-a*fvar(i)

		do j=1,size(fvar2)

			if ((fvar2(j).ge.fvar(i-1)).and.(fvar2(j).lt.fvar(i)))then
		      var2(k,j)=a*fvar2(j)+b

			end if
		end do
	end do
end do


end subroutine interp2D
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine average2D(var,var2,var3,fvar)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Vertical averaging over x meters. 
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

  implicit none

	real, dimension(:,:), allocatable :: var,fvar
	real, dimension(:), allocatable :: var3
	integer i,j,k,a,pas
	real b,var2

print *,'vertical averaging over',var2,'km'

a=floor((var3(size(var3))+var2)/var2)

allocate(fvar(size(var,1),a))

pas=(var2*1000/20)

a=1

do i=1,size(var,1)

	do j=1,size(var,2),pas
	fvar(i,a)=sum(var(i,j:j+pas:1))/(pas+1)
	a=a+1
	enddo
a=1

enddo

end subroutine average2D
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
subroutine interp2Delevation(var,var1,var2,var3,fvar)
!interp2Delevation(filtreelevation,ATBatlid,altValues,altValues20,ATBatlid20elevation)
 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Vertical interpolation over 20 meters without elevation points. 
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

 implicit none

	real, dimension(:,:), allocatable :: var, var1,fvar
	real, dimension(:), allocatable :: var2, var3
	integer i,j,k
	real a,b

print *,"data interpolation"
allocate(fvar(size(var,1),size(var3)))

var1(:,:)=var1(:,:)*var(:,:)


do k=1,size(var,1)

	do i=1,size(var,2)
	
		if (var1(k,i) .ne. 0) then
	
		j=i
		exit
		endif
	enddo

	do i=1,j
	fvar(k,i)=0
	enddo

	do i=j+1,size(var,2)
	a=(var1(k,i)-var1(k,i-1))/(var2(i)-var2(i-1))
	b=var1(k,i)-a*var2(i)

		do j=1,size(var3)

			if ((var3(j).ge.var2(i-1)).and.(var3(j).lt.var2(i)))then
			fvar(k,j)=a*var3(j)+b
			end if
		end do
	end do
end do

end subroutine interp2Delevation
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
subroutine average2Delevation(var,var2,var3,var4,var5,fvar)
!average2Delevation(ATBatlid20elevation,vertTemp,altValues20,altValuesaverage,elevation,ATBatlidmoyelevation)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Vertical averaging over x meters without elevation points. 
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

  implicit none

	real, dimension(:,:), allocatable :: var,fvar
	real, dimension(:), allocatable :: var3, var4,var5,tmpvar
	integer i,j,k,ii,jj,a,pas,pas2
	real b,var2

print *,'vertical averaging over',var2,'km'


allocate(fvar(size(var,1),size(var4)))
allocate(tmpvar(size(var3)))

pas=(var2*1000)
pas2=(var2*1000)/20
tmpvar=var3*1000

do i=1,size(var,1)
a=1	


	do k=1,size(var,2)
		if (var(i,k) .ne. 0) then
		jj=k
		exit
		endif
	enddo

	b=floor(tmpvar(jj)/pas)

	if (tmpvar(jj) .lt. pas) then
		do k=1,size(var,2)
			if (var(i,k) .ne. 0) then
			j=k
			exit
			endif
		enddo

		do ii=j,pas2,pas2
			fvar(i,a)=sum(var(i,ii:pas2:1))/((pas-tmpvar(ii))/20)
			a=a+1
			
		enddo

		do j=pas2+1,size(var,2),pas2
			fvar(i,a)=sum(var(i,j:j+pas2:1))/(pas2)
			a=a+1
			
		enddo

	endif

	if (tmpvar(jj) .gt. pas) then
		do j=1,b*pas2,pas2
			fvar(i,a)=0
			a=a+1
			
		enddo
		do j=b*pas2+1,(b+1)*pas2,pas2
			fvar(i,a)=sum(var(i,j:(b+1)*pas2:1))/(((b+1)*pas-tmpvar(j))/20)
			a=a+1
		enddo

		do j=(b+1)*pas2+1,size(var,2),pas2
			fvar(i,a)=sum(var(i,j:j+pas2:1))/(pas2)
			a=a+1
		enddo
	endif
enddo


end subroutine average2Delevation
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine surfaceelevation(var,var2,fvar,fvar2)
!surfaceelevation(altValues,ATBatlid,filtreelevation,elevation)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Reject points below the surface elevation 
!NB: may change in futur version
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

	 implicit none

	real, dimension(:), allocatable :: var, surface,fvar2
	real, dimension(:,:), allocatable :: var2,fvar
	integer i,j,k,kk
	real a,b

print *,"surface elevation detection"
allocate(fvar(size(var2,1),size(var2,2)))
allocate(surface(size(var2,1)))
allocate(fvar2(size(var2,1)))

!PARTIE A SUPPRIMER DANS LE FUTUR
!CONSTRUCTION D'UNE ELEVATION FICTIVE

surface(:)=var2(:,1)*1000.

do k=1,size(surface,1)

if (surface(k) .lt. 10) then
surface(k)=surface(k)+100
endif

do while (surface(k) .gt. 2000)
surface(k)=surface(k)/10
enddo

end do

!FIN PARTIE A SUPPRIMER DANS LE FUTUR
!CONSTRUCTION D'UNE ELEVATION FICTIVE


fvar2=surface

do k=1,size(surface,1)    

	do kk=1,size(var,1)

    
		if (surface(k).gt.var(kk)*1000) then
		fvar(k,kk)=0
		endif

		if (surface(k).lt.var(kk)*1000) then
		fvar(k,kk)=1
		endif

	end do

end do


end subroutine surfaceelevation
!----------------------------------------------------------------------------!

subroutine rejectionfiltering(var,var2,var3,fvar,fvar2)
!rejectionfiltering(ATBatlidmoy,ATBmolatlidmoy,altValuesaverage,ATBvert,ATBmolvert)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Slide window averaging over horizontal profiles 
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

implicit none

	real, dimension(:), allocatable :: var3, res, res2
	real, dimension(:,:), allocatable :: var, var2, fvar, fvar2
	integer i,j
	real n

print *,"Filtering signal in stratosphere"


 allocate(fvar(size(var,1),size(var,2)))
 allocate(fvar2(size(var,1),size(var,2)))



n=18

do j=1,size(var,2)
 do i=1,size(var,1)

fvar(i,j)=0.
fvar2(i,j)=0.
	

! A changer pour se mettre dans la strato. 
	if (var3(j) .ge. 11.9) then


			if (i .gt. n .and. i .lt. (size(var,1)-n+1)) then

			fvar(i,j)=sum(var(i-n:i+n:1,j))/(2*n+1)
			fvar2(i,j)=sum(var2(i-n:i+n:1,j))/(2*n+1)

			endif
	endif

	enddo


		do i=1,n

if (var3(j) .ge. 11.9) then
fvar(i,j)=fvar(n+1,j)
fvar2(i,j)=fvar2(n+1,j)
end if
		enddo


		do i=size(var,1)-n,size(var,1)
if (var3(j) .ge. 11.9) then
fvar(i,j)=fvar(i-1,j)
fvar2(i,j)=fvar2(i-1,j)
end if
		enddo


enddo

end subroutine rejectionfiltering
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine lidarwrite2D(rep,filenamewrite,var,xvar,yvar)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Write 2D results in files. 
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

	use netcdf
	implicit none

	real, dimension(:), allocatable :: xvar,yvar
 	real, dimension(:,:), allocatable :: var
	integer, dimension(2) :: numatbid
   	character*100 filenamewrite, repdata, rep
	integer status, ncid,numx,numy
	integer numxid,numyid,varid,varid2,varid3

repdata=trim(rep)//filenamewrite

print *,'Write lidar netcdf files'
print *,trim(repdata)

numx=size(xvar)
numy=size(yvar)

status=nf90_create(repdata, nf90_clobber, ncid)
status=nf90_def_dim(ncid,"x",numx,numxid)
status=nf90_def_dim(ncid,"y",numy,numyid)

numatbid=(/ numxid, numyid /)

status=nf90_def_var(ncid, "lon", nf90_float,numxid, varid)
status=nf90_def_var(ncid, "alt", nf90_float,numyid, varid2)
status=nf90_def_var(ncid,"ATBatlid",nf90_float,numatbid,varid3)

status=nf90_enddef(ncid)

status=nf90_put_var(ncid,varid,xvar)
status=nf90_put_var(ncid,varid2,yvar)
status=nf90_put_var(ncid,varid3,var)

status=nf90_close(ncid)

end subroutine lidarwrite2D
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
subroutine lidarwrite1D(rep,filenamewrite,var,xvar)

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!Write 1D results in files. 
!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

	use netcdf
	implicit none

	real, dimension(:),allocatable :: xvar,var
   	character*100 filenamewrite, repdata, rep
	integer status, ncid,numx
	integer numxid,varid,varid2

repdata=trim(rep)//filenamewrite	


print *,'Write lidar netcdf files'
print *,trim(repdata)
status=nf90_create(repdata, nf90_clobber, ncid)

numx=size(xvar)
status=nf90_def_dim(ncid,"x",numx,numxid)


status=nf90_def_var(ncid, "lon", NF90_FLOAT,numxid, varid)
status=nf90_def_var(ncid,"pmol",NF90_FLOAT,numxid,varid2)

status=nf90_enddef(ncid)


status=nf90_put_var(ncid,varid,xvar)
status=nf90_put_var(ncid,varid2,var)

status=nf90_close(ncid)


end subroutine lidarwrite1D
!----------------------------------------------------------------------------!

!****************************************************************************!

!----------------------------------------------------------------------------!
!subroutine interp1D(var,var2,fvar,fvar2)
!  implicit none
!
!	real, dimension(:), allocatable :: var
!	real, dimension(:), allocatable :: fvar, fvar2, var2
!	integer i,j
!	real a,b
!
!print *,"data interpolation"
!allocate(var2(size(fvar2)))
!
!do i=2,size(fvar)
!	a=(var(i)-var(i-1))/(fvar(i)-fvar(i-1))
!	b=var(i)-a*fvar(i)
!
!	do j=1,size(fvar2)
!
!		if ((fvar2(j).ge.fvar(i-1)).and.(fvar2(j).lt.fvar(i)))then
!		var2(j)=a*fvar2(j)+b
!		end if
!	end do
!end do
!
!end subroutine interp1D
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
!subroutine check(istatus)
!
!	use netcdf
!	implicit none
!	
! 
!	integer, intent (in) :: istatus
!
!if (istatus /= nf90_noerr) then
!write (*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
!end if
!
!
!end subroutine check
!----------------------------------------------------------------------------!

end program AGCM

!****************************************************************************!
!                                                                            !
!                               END PROGRAM                                  !
!                                                                            !
!****************************************************************************!



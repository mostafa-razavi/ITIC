program	main
implicit double precision(A-H,O-Z)


character :: zures_datafile*88,ITIC_conditions*88
character(len=100) :: arg1,arg2
parameter(maxICs=11)
parameter(maxICTs=5)
parameter(nPointsOnIC=3)
parameter(nMaxData=150)


double precision :: Tfile(nMaxData),rhofile(nMaxData),Zfile(nMaxData),Uresfile(nMaxData)

double precision :: Ures_IT(20)
double precision :: Ures_IT_vr(20), Ures_IT_vr2(20)
double precision :: Ures_IC(maxICs,maxICTs)	
Double Precision :: T_IC_calc(maxICs,maxICTs)


double precision :: T_IT(20),rho_IT(0:20),Z_IT(20)
double precision :: T_IT_vr(20),rho_IT_vr(0:20),Z_IT_vr(20)
double precision :: T_IT2_vr(20),rho_IT2_vr(0:20),Z_IT2_vr(20), uDep_over_rho_vr(4)


double precision :: T_IC(maxICs,maxICTs),rho_IC(maxICs),Z_IC(maxICs,maxICTs)
integer :: convergeStatus(maxICs)

logical calculateVirialSupercritical,calculateVirialSubcritical

double precision :: Zmin1OverRho_IT(0:15),aDep_IT(0:15),uDepT_IT(15),ThousandOverT_IT(15)
double precision :: Zmin1OverRho_IC(maxICs,maxICTs),aDep_IC(maxICs,maxICTs),uDepT_IC(maxICs,maxICTs), &
					ThousandOverT_IC(maxICs,maxICTs)
double precision :: multi_IT(15), multi_IC(maxICs,maxICTs)
double precision :: Tsat(maxICs),uDepTsat(maxICs),aDepSat(maxICs), TsatLnDev(MaxICs)
double precision :: rhoV(MaxICs),rhoVp(MaxICs),rhoVpp(MaxICs),zLiq(MaxICs),Psat(MaxICs),Hvap(MaxICs)
double precision :: B2sat(MaxICs),B3sat(MaxICs)
double precision :: dB2_dBeta(maxICs)
double precision :: T_IT_calc(50),rho_IT_calc(50),rho_IC_calc(50)
double precision :: rho_IT_calc_vr(15)

real,dimension (3)::bmatrix
real,dimension (3)::xmatrix
real,dimension (3,3)::amatrix

character :: dumString1*50,dumString2*50,dumString3*50,String*50




maxTolerance = 1.0d0


call getarg(1,arg1)
call getarg(2,arg2)
zures_datafile = arg1
ITIC_conditions = arg2


open(2231,file=ITIC_conditions)
nItPts = 9
nICs = 5
do loopline=1,200	
	read(2231,*,ioStat=ioErr) dumString1,dumString2
	String = trim(dumString1)
	if(String.eq."RHO_HIGH:") then
		read(dumString2,*,iostat=ioErr)highestRho
	elseif(String.eq."MW:") then
		read(dumString2,*,iostat=ioErr) MW	
	elseif(String.eq."TC:") then
		read(dumString2,*,iostat=ioErr) TC
	elseif(String.eq."T_IT:") then
		backspace 2231
		read(2231,*,ioStat=ioErr) dumString1,dumString2, dumString3
		read(dumString3,*,iostat=ioErr) T_IT_subcritical
	elseif(String.eq."T_HIGH:") then
		read(dumString2,*,iostat=ioErr)HighestT
	elseif(String.eq."T_IC1:") then
		write(*,*)
		iFirstIcPt=nItPts-nICs+1
		read(dumString2,*,iostat=ioErr)T_IC_calc(iFirstIcPt,1)
		do i=iFirstIcPt+1,nItPts
			read(2231,*,ioStat=ioErr) dumString1,dumString2
			read(dumString2,*,iostat=ioErr)T_IC_calc(i,1)
		enddo
	endif
end do
close(2231)	

SubCritReducedTemp = T_IT_subcritical / TC


rho_IT_calc(9)=highestRho
rho_increment=highestRho/7.0
rho_IT_calc(1)=rho_increment
rho_IT_calc(2)=rho_IT_calc(1)+rho_increment
rho_IT_calc(3)=rho_IT_calc(2)+rho_increment
rho_IT_calc(4)=rho_IT_calc(3)+rho_increment
rho_IT_calc(5)=rho_IT_calc(4)+rho_increment
rho_IT_calc(6)=rho_IT_calc(5)+rho_increment*0.5
rho_IT_calc(7)=rho_IT_calc(6)+rho_increment*0.5
rho_IT_calc(8)=rho_IT_calc(7)+rho_increment*0.5

rho_IT_calc_vr(4)=highestRho/7.0
rho_IT_calc_vr(3)=highestRho/14.0
rho_IT_calc_vr(2)=highestRho/21.0
rho_IT_calc_vr(1)=highestRho/28.0

do i=iFirstIcPt,nItPts
	rec_T_increment=(1000.0/T_IC_calc(i,1)-1000.0/HighestT)/(nPointsOnIC-1)
	do j=2,nPointsOnIC
		T_IC_calc(i,j)=1000.0/(1000.0/T_IC_calc(i,1)-(j-1)*rec_T_increment)
	enddo
enddo

do i=iFirstIcPt,nItPts
	rho_IC_calc(i)=rho_IT_calc(i)
enddo

do i=1,nItPts
	T_IT_calc(i)=HighestT
enddo



write(*,*) 
write(*,*)'===============================Reading From Simulator Output====================================='
write(*,*)
write(*,'(A2,1x,A9,1x,A9,1x,2(A9,1x),2(A15,1x),A6)')'i','T(K)','rho(g/ml)','Z','Zstd','Ures', 'UresStd','nMolec'

open(4321,file=zures_datafile)
read(4321,*)
do i=1,nMaxData	
	read(4321,*,ioStat=ioErr) Tfile(i),rhofile(i),Zfile(i), Uresfile(i)
	if(ioErr == -1)then !End of file reached
		nData=i-1
		exit
	endif
	write(*,'(I2,1x,f9.2,1x,f9.5,1x,2(f9.3,1x),2(f15.3,1x),f6.1)') &
		i, Tfile(i),rhofile(i),Zfile(i), Uresfile(i)
enddo
close(4321)

write(*,*)
write(*,*)'===============================Splitting Data into IT and IC Arrays====================================='
write(*,*)
do i=1,nData
	do j=1,nItPts
		tolerance=(abs(Tfile(i)-T_IT_calc(j))/Tfile(i)*100.0+abs(rhofile(i)-rho_IT_calc(j))/rhofile(i)*100.0)
		if(tolerance .lt. maxTolerance) then
			T_IT(j)=Tfile(i)
			rho_IT(j)=rhofile(i)
			Z_IT(j)=Zfile(i)
			Ures_IT(j) = Uresfile(i)
		exit
		endif
	enddo
enddo

do i=1,nData
	do j=iFirstIcPt,nItPts
		do k=1,nPointsOnIC
			tolerance=(abs(Tfile(i)-T_IC_calc(j,k))/Tfile(i)*100.0+abs(rhofile(i)-rho_IC_calc(j))/rhofile(i)*100.0)
			if(tolerance .lt. maxTolerance) then
				T_IC(j,k)=Tfile(i)
				rho_IC(j)=rhofile(i)
				Z_IC(j,k)=Zfile(i)
				Ures_IC(j,k) = Uresfile(i)
			exit
			endif
		enddo
	enddo
enddo

do i=1,nData
	do j=1,4
		tolerance=(abs(Tfile(i)-HighestT)/Tfile(i)*100.0+abs(rhofile(i)-rho_IT_calc_vr(j))/rhofile(i)*100.0)
		if(tolerance .lt. maxTolerance) then
			T_IT_vr(j)=Tfile(i)
			rho_IT_vr(j)=rhofile(i)
			Z_IT_vr(j)=Zfile(i)
			Ures_IT_vr(j) = Uresfile(i)
		exit
		endif
	
	enddo
enddo

do i=1,nData
	do j=1,4
		tolerance=(abs(Tfile(i)-TC*SubCritReducedTemp)/Tfile(i)*100.0+abs(rhofile(i)-rho_IT_calc_vr(j))/rhofile(i)*100.0)
		if(tolerance .lt. maxTolerance) then
			T_IT2_vr(j)=Tfile(i)
			rho_IT2_vr(j)=rhofile(i)
			Z_IT2_vr(j)=Zfile(i)
			Ures_IT_vr2(j) = Uresfile(i)
		exit
		endif
	
	enddo
enddo

rho1dif = abs(rho_IT_vr(1)-rho_IT_calc_vr(1))/rho_IT_calc_vr(1)*100
rho2dif = abs(rho_IT_vr(2)-rho_IT_calc_vr(2))/rho_IT_calc_vr(2)*100
rho3dif = abs(rho_IT_vr(3)-rho_IT_calc_vr(3))/rho_IT_calc_vr(3)*100
if ( rho1dif .le. 1.0 .AND. rho2dif .le. 1.0 .AND. rho3dif .le. 1.0) then
	calculateVirialSupercritical = .true. 
else
	calculateVirialSupercritical = .false. 
endif

rho1dif2 = abs(rho_IT2_vr(1)-rho_IT_calc_vr(1))/rho_IT_calc_vr(1)*100
rho2dif2 = abs(rho_IT2_vr(2)-rho_IT_calc_vr(2))/rho_IT_calc_vr(2)*100
rho3dif2 = abs(rho_IT2_vr(3)-rho_IT_calc_vr(3))/rho_IT_calc_vr(3)*100
if ( rho1dif2 .le. 1.0 .AND. rho2dif2 .le. 1.0 .AND. rho3dif2 .le. 1.0) then
	calculateVirialSubcritical = .true. 
else
	calculateVirialSubcritical = .false. 
endif


write(*,*) "Isothermic Points"
write(*,'(A2,1x,A9,1x,A9,1x,2(A15,1x),A6)') 'i','T(K)','rho(g/ml)','Z','Ures','nMolec'				
do i=1,nItPts
	multi_IT(i)=T_IT(i)*rho_IT(i)
	if(multi_IT(i) < 1e-11) then
		write(*,'(A31,1x,I2,A6,1x,2F10.5)') "Warning: Data lacks row number",i,'T,rho:',T_IT_calc(i),rho_IT_calc(i)
		cycle
	endif
	write(*,'(I2,1x,f9.2,1x,f9.5,1x,2(f15.6,1x),f6.1)') i,T_IT(i),rho_IT(i),Z_IT(i),Ures_IT(i)
enddo

do j=iFirstIcPt,nItPts
	write(*,*)
	write(*,'(A28,2x,I2)') "Isochoric Points on Isochore",j
	write(*,'(A2,1x,A9,1x,A9,1x,2(A15,1x),A6)')'i','T(K)','rho(g/ml)','Z','Ures', 'nMolec'
	do k=1,nPointsOnIC
		multi_IC(j,k)=T_IC(j,k)*rho_IC(j)
		if(multi_IC(j,k) < 1e-11) then
			write(*,'(A31,1x,2I2,A6,1x,2F10.5)') "Warning: Data lacks row number",j,k,'T,rho:',T_IC_calc(j,k),rho_IC_calc(j)
			cycle
		endif
		write(*,'(I2,1x,f9.2,1x,f9.5,1x,2(f15.6,1x),f6.1)') k,T_IC(j,k),rho_IC(j),Z_IC(j,k),Ures_IC(j,k)
	enddo
enddo

if(calculateVirialSupercritical) then
write(*,*)
write(*,*) "Supercritical Virial Calculation Points"
write(*,'(A2,1x,A9,1x,A9,1x,2(A15,1x),A6)')'i','T(K)','rho(g/ml)','Z','Ures','nMolec'				
do i=1,4
	multi_IT(i)=T_IT_vr(i)*rho_IT_vr(i)
	if(multi_IT(i) < 1e-11) then
		write(*,'(A31,1x,I2,A6,1x,2F10.5)') "Warning: Data lacks row number",i,'T,rho:',HighestT,rho_IT_calc_vr(i)
		cycle
	endif
	write(*,'(I2,1x,f9.2,1x,f9.5,1x,2(f15.7,1x),f6.1)')&
	i,T_IT_vr(i),rho_IT_vr(i),Z_IT_vr(i),Ures_IT_vr(i)

enddo
endif

if(calculateVirialSubcritical) then
write(*,*)
write(*,*) "Subcritical Virial Calculation  Points"
write(*,'(A2,1x,A9,1x,A9,1x,2(A15,1x),A6)')'i','T(K)','rho(g/ml)','Z','Ures','nMolec'				

do i=1,4
	multi_IT(i)=T_IT2_vr(i)*rho_IT2_vr(i)
	if(multi_IT(i) < 1e-11) then
		write(*,'(A31,1x,I2,A6,1x,2F10.5)') "Warning: Data lacks row number",i,'T,rho:',T_IT2_vr(i),rho_IT_calc_vr(i)
		cycle
	endif
	write(*,'(I2,1x,f9.2,1x,f9.5,1x,2(f15.7,1x),f6.1)')&
	i,T_IT2_vr(i),rho_IT2_vr(i),Z_IT2_vr(i),Ures_IT_vr2(i)

enddo
endif

do i=1,nItPts
	Zmin1OverRho_IT(i)=(Z_IT(i)-1.0)/rho_IT(i)
	uDepT_IT(i) = Ures_IT(i) * T_IT(i)
	ThousandOverT_IT(i)=1000.d0/T_IT(i)
enddo

do j=iFirstIcPt,nItPts
	do k=1,nPointsOnIC
		Zmin1OverRho_IC(j,k)=(Z_IC(j,k)-1.0)/rho_IC(j)
		uDepT_IC(j,k)= Ures_IC(j,k) * T_IC(j,k)
		ThousandOverT_IC(j,k)=1000.d0/T_IC(j,k)
	enddo
enddo	


!========================================Determining B2 Correlation=====================================

write(*,*)
write(*,*)"calculateVirialSupercritical?",calculateVirialSupercritical
write(*,*)"calculateVirialSubcritical?",calculateVirialSubcritical

if(calculateVirialSupercritical) then
	Y1 = (Z_IT_vr(1) - 1.0)/rho_IT_vr(1)
	Y2 = (Z_IT_vr(2) - 1.0)/rho_IT_vr(2)
	Y3 = (Z_IT_vr(3) - 1.0)/rho_IT_vr(3)
	Y4 = (Z_IT_vr(4) - 1.0)/rho_IT_vr(4)

	X1 = rho_IT_vr(1)
	X2 = rho_IT_vr(2)
	X3 = rho_IT_vr(3)
	X4 = rho_IT_vr(4)

	call getInterceptSlope(X1,X2,X3,X4,Y1,Y2,Y3,Y4,B2atHighestT,B3atHighestT)
endif

if(calculateVirialSupercritical .AND. calculateVirialSubcritical) then

	Y1 = (Z_IT2_vr(1) - 1.0)/rho_IT2_vr(1)
	Y2 = (Z_IT2_vr(2) - 1.0)/rho_IT2_vr(2)
	Y3 = (Z_IT2_vr(3) - 1.0)/rho_IT2_vr(3)
	Y4 = (Z_IT2_vr(4) - 1.0)/rho_IT2_vr(4)

	X1 = rho_IT2_vr(1)
	X2 = rho_IT2_vr(2)
	X3 = rho_IT2_vr(3)
	X4 = rho_IT2_vr(4)

	call getInterceptSlope(X1,X2,X3,X4,Y1,Y2,Y3,Y4,B2atSubcriticalIT,B3atSubcriticalIT)

	do i=1,4
		uDep_over_rho_vr(i)=Ures_IT_vr2(i)/rho_IT2_vr(i)
	enddo

	Y1 = uDep_over_rho_vr(1)
	Y2 = uDep_over_rho_vr(2)
	Y3 = uDep_over_rho_vr(3)
	Y4 = uDep_over_rho_vr(4)

	X1 = rho_IT2_vr(1)
	X2 = rho_IT2_vr(2)
	X3 = rho_IT2_vr(3)
	X4 = rho_IT2_vr(4)

	call getInterceptSlope(X1,X2,X3,X4,Y1,Y2,Y3,Y4,uDepOverRhoVsRhoIntercept,uDepOverRhoVsRhoSlope)
	write(*,*)uDepOverRhoVsRhoIntercept,uDepOverRhoVsRhoSlope
	TofSubCriticalIT=TC*SubCritReducedTemp

	amatrix(1,1)=0.0
	amatrix(1,2)=1.0/HighestT - 1.0/TofSubCriticalIT
	amatrix(1,3)=1.0/HighestT**3 - 1.0/TofSubCriticalIT**3

	amatrix(2,1)=0.0
	amatrix(2,2)=1.0/TofSubCriticalIT
	amatrix(2,3)=3.0/TofSubCriticalIT**3

	amatrix(3,1)=1.0
	amatrix(3,2)=1.0/HighestT
	amatrix(3,3)=1.0/HighestT**3

	bmatrix(1)=B2atHighestT - B2atSubcriticalIT
	bmatrix(2)=uDepOverRhoVsRhoIntercept
	bmatrix(3)=B2atHighestT

	xmatrix(3)=(bmatrix(2)-amatrix(2,2)*bmatrix(1)/amatrix(1,2))/(amatrix(2,3)-amatrix(1,3)*amatrix(2,2)/amatrix(1,2))
	xmatrix(2)=bmatrix(1)/amatrix(1,2)-amatrix(1,3)/amatrix(1,2)*xmatrix(3)
	xmatrix(1)=bmatrix(3)-amatrix(3,2)*xmatrix(2)-amatrix(3,3)*xmatrix(3)

	write(*,*)
	write(*,*) "Virial A coef=",xmatrix(1)
	write(*,*) "Virial B coef=",xmatrix(2)
	write(*,*) "Virial C coef=",xmatrix(3)
	Acoef = xmatrix(1)
	Bcoef = xmatrix(2)
	Ccoef = xmatrix(3)
	B2atHighestT = Acoef + Bcoef/HighestT + Ccoef / HighestT**3
endif


write(*,'(A9,1x,2(F10.4,3x))')"B2_IT:",B2atHighestT,HighestT
write(*,'(A9,1x,2(F10.4,3x))')"B2_IT2:",B2atSubcriticalIT,TofSubCriticalIT

Zmin1OverRho_IT(0)=B2atHighestT

aDep_IT(0)= 0.0
aDep_IT(2)= (rho_IT(2)-rho_IT(0))/6.0*(Zmin1OverRho_IT(0)+4.0*Zmin1OverRho_IT(1)+Zmin1OverRho_IT(2))
aDep_IT(3)= 3.0/8.0*(rho_IT(3))/3.0*(Zmin1OverRho_IT(0)+3.0*Zmin1OverRho_IT(1)+3.0*Zmin1OverRho_IT(2)+Zmin1OverRho_IT(3))
aDep_IT(1)= -1.0*(rho_IT(3)-rho_IT(1))/6.0*(Zmin1OverRho_IT(1)+4.0*Zmin1OverRho_IT(2)+Zmin1OverRho_IT(3))+aDep_IT(3)
if(nItPts .eq. 9) then
	aDep_IT(4)= (rho_IT(4)-rho_IT(2))/6.0*(Zmin1OverRho_IT(2)+4.0*Zmin1OverRho_IT(3)+Zmin1OverRho_IT(4))+aDep_IT(2)
	aDep_IT(5)= 3.0/8.0*(rho_IT(5)-rho_IT(2))/3.0*(Zmin1OverRho_IT(2)+3.0*Zmin1OverRho_IT(3)+3.0*Zmin1OverRho_IT(4)+&
			Zmin1OverRho_IT(5))+aDep_IT(2)
	aDep_IT(7)= (rho_IT(7)-rho_IT(4))/6.0*(Zmin1OverRho_IT(4)+4.0*Zmin1OverRho_IT(5)+Zmin1OverRho_IT(7))+aDep_IT(4)
	aDep_IT(8)= 3.0/8.0*(rho_IT(8)-rho_IT(5))/3.0*(Zmin1OverRho_IT(5)+3.0*Zmin1OverRho_IT(6)+3.0*Zmin1OverRho_IT(7)+&
			Zmin1OverRho_IT(8))+aDep_IT(5)
	aDep_IT(6)= -1.0*(rho_IT(8)-rho_IT(6))/6.0*(Zmin1OverRho_IT(6)+4.0*Zmin1OverRho_IT(7)+Zmin1OverRho_IT(8))+aDep_IT(8)
	aDep_IT(9)= (rho_IT(9)-rho_IT(7))/6.0*(Zmin1OverRho_IT(7)+4.0*Zmin1OverRho_IT(8)+Zmin1OverRho_IT(9))+aDep_IT(7)
elseif(nItPts .eq. 11) then
	aDep_IT(5)= (rho_IT(5)-rho_IT(2))/6.0*(Zmin1OverRho_IT(2)+4.0*Zmin1OverRho_IT(3)+Zmin1OverRho_IT(5))+aDep_IT(2)
	aDep_IT(7)= 3.0/8.0*(rho_IT(7)-rho_IT(2))/3.0*(Zmin1OverRho_IT(2)+3.0*Zmin1OverRho_IT(3)+3.0*Zmin1OverRho_IT(5)+&
			Zmin1OverRho_IT(7))+aDep_IT(2)
	aDep_IT(9)= (rho_IT(9)-rho_IT(5))/6.0*(Zmin1OverRho_IT(5)+4.0*Zmin1OverRho_IT(7)+Zmin1OverRho_IT(9))+aDep_IT(5)
	aDep_IT(10)= 3.0/8.0*(rho_IT(10)-rho_IT(7))/3.0*(Zmin1OverRho_IT(7)+3.0*Zmin1OverRho_IT(8)+3.0*Zmin1OverRho_IT(9)+&
			Zmin1OverRho_IT(10))+aDep_IT(7)
	aDep_IT(8)= -1.0*(rho_IT(10)-rho_IT(8))/6.0*(Zmin1OverRho_IT(8)+4.0*Zmin1OverRho_IT(9)+Zmin1OverRho_IT(10))+aDep_IT(10)
	aDep_IT(11)= (rho_IT(11)-rho_IT(9))/6.0*(Zmin1OverRho_IT(9)+4.0*Zmin1OverRho_IT(10)+Zmin1OverRho_IT(11))+aDep_IT(9)
	aDep_IT(6)= -1.0*(rho_IT(8)-rho_IT(6))/6.0*(Zmin1OverRho_IT(6)+4.0*Zmin1OverRho_IT(7)+Zmin1OverRho_IT(8))+aDep_IT(8)
	aDep_IT(4)= -1.0*(rho_IT(6)-rho_IT(4))/6.0*(Zmin1OverRho_IT(4)+4.0*Zmin1OverRho_IT(5)+Zmin1OverRho_IT(6))+aDep_IT(6)
endif

do j=iFirstIcPt,nItPts
	aDep_IC(j,1)= ( ThousandOverT_IC(j,1)-ThousandOverT_IC(j,3) )/6.d0* &
			( uDepT_IC(j,3)+4.d0*uDepT_IC(j,2)+ uDepT_IC(j,1) )/1000.d0+aDep_IT(j)
	if(nPointsOnIC==4) then
	aDep_IC(j,1)= ( ThousandOverT_IC(j,1)-ThousandOverT_IC(j,4) )/8.d0* &
			( uDepT_IC(j,4)+3.d0*uDepT_IC(j,3)+3.d0*uDepT_IC(j,2)+ uDepT_IC(j,1) )/1000.d0+aDep_IT(j)
	endif
enddo

write(*,*)
write(*,*)'===============================Psat Calculation====================================='
write(*,*)

iExpoErr=0
maxIter=500
convergeStatus=-2
do i=iFirstIcPt,nItPts

	tempPsat=Psat(i) !tempPsat is used for Converge/Diverge check
	tempRhov=rhoV(i) !tempRhoV is used for Converge/Diverge check

	WRITE(*,*)
	Write(*,'(A,I2,A,F8.4,A)') "Iterations for Isochore",i,":",rho_IC(i)," (g/ml)"
	write(*,'(2A3,50A12)')"j","i","Tsat","Psat","rhoV","B2sat",&
			"aDepSat","uDep","zLiq","zVap","Hvap","drhoV","d2rhoV"	
	do j=1,maxIter
		if(j==1)then
			zLiq(i)=0.001
		else
			zLiq(i)=Psat(i)/( rho_IC(i)/MW*8.3144598d0*Tsat(i) )
		endif

		Tsat(i)=TsatFinder(zLiq(i),1000.0/T_IC(i,3),Z_IC(i,3),1000.0/T_IC(i,2),Z_IC(i,2) &
			,1000.0/T_IC(i,1),Z_IC(i,1))


		X1 = T_IC(i,1)
		X2 = T_IC(i,2)
		X3 = T_IC(i,3)
		X4 = T_IC(i,3)
		Y1 = uDepT_IC(i,1)
		Y2 = uDepT_IC(i,2)
		Y3 = uDepT_IC(i,3)
		Y4 = uDepT_IC(i,3)
		call getInterceptSlope(X1,X2,X3,X4,Y1,Y2,Y3,Y4,YINTERCEPT,SLOPE)
		uDepTsat(i) = (SLOPE * Tsat(i) + YINTERCEPT)

		if(j==1)then
			aDepSat(i)= ( 1000.d0/Tsat(i)-ThousandOverT_IC(i,1) )*( uDepTsat(i)+uDepT_IC(i,1) )/2000.d0 + aDep_IC(i,1)
		else
			aDepSat(i)= ( 1000.d0/Tsat(i)-1000.d0/TsatOld )*( uDepTsat(i)+uDepSatOld )/2000.d0 + aDepSat(i)
		endif
		uDepSatOld = uDepTsat(i)
		TsatOld = Tsat(i)
		

		if(calculateVirialSupercritical .AND. calculateVirialSubcritical) then
			B2sat(i) = Acoef + Bcoef/Tsat(i) + Ccoef / Tsat(i)**3
		endif

		expArgAZ1=aDepSat(i)+zLiq(i)-1.d0
		if(expArgAZ1 > 33)iExpoErr=1
		if(j==1)then
			expArgB2=-2.d0*B2sat(j)*rho_IC(j)*exp(expArgAZ1) 
			expArgB3=-1.5d0*B3sat(j)*(rho_IC(j)*exp(expArgAZ1))**2
		else
			expArgB2=-2.d0*B2sat(i)*rhoV(i)
			expArgB3=-1.5d0*B3sat(i)*rhoV(i)**2
		endif

		rhoV(i) = rho_IC(i)*exp(expArgAZ1)*exp(expArgB2)*exp(expArgB3)*exp(expArgB4)*exp(expArgB5)*exp(expArgB6)
		rhoVp(i) = rhoV(i)*(-2.d0*B2sat(i) - 3.d0*B3sat(i)*rhoV(i))
		rhoVpp(i) = rhoV(i)*( -3.d0*B3sat(i) + rhoVp(i)*( -2.d0*B2sat(i) - 3.d0*B3sat(i)*rhoV(i)))
		zVsat = 1.d0+B2sat(i)*rhoV(i)+B3sat(i)*rhoV(i)**2
		Psat(i) = ( zVsat )*rhoV(i)/MW*8.3144598d0*Tsat(i)
		denom = ( 1.0/Tsat(i) - 1.0/HighestT ) 
		dB2_dBeta(i) = (B2sat(i) - B2atHighestT)/denom
		hSatVap = dB2_dBeta(i)*rhoV(i)/Tsat(i)+zVsat-1	
		hSatLiq = uDepTSat(i)/Tsat(i)+zLiq(i)-1
		Hvap(i) = ( hSatVap - hSatLiq)*8.3144598d0*Tsat(i)/1000.0

		write(*,'(2I3,F12.5,F15.9,50F12.7)')j,i,Tsat(i),Psat(i),rhoV(i),B2sat(i),&
			aDepSat(i),uDepTsat(i)/Tsat(i),zLiq(i),zVsat,Hvap(i),rhoVp(i),rhoVpp(i) &
			,aDepSat(i) - aDep_IC(i,1)

		!<====Convergence Check==========================================
		PsatIncrement=abs((tempPsat-Psat(i))/Psat(i)*100.0)
		rhoVIncrement=abs((tempRhoV-rhoV(i))/rhoV(i)*100.0)
		if(PsatIncrement .le. 1e-3 .and. rhoVIncrement .le. 1e-3) then
			if(rhoVp(i) .lt. 1.d0 .AND. rhoVpp(i) .gt. 0.d0) then
				convergeStatus(i)=1
				Write(*,'(A,I2,A,F8.4,A)') "Isochore",i,":",rho_IC(i)," (g/ml):"
				write(*,'(A19,I3,A)') "Convergence Status= ",convergeStatus(i)," (Converged to first root)"				
				exit
			elseif(rhoVp(i) .gt. 1.d0 .AND. rhoVpp(i) .gt. 0.d0)then
				convergeStatus(i)=2
				Write(*,'(A,I2,A,F8.4,A)') "Isochore",i,":",rho_IC(i)," (g/ml):"
				write(*,'(A19,I3,A)') "Convergence Status= ",convergeStatus(i)," (Converged to second root)"				
				exit
			elseif(rhoVpp(i) .lt. 0.d0) then
				convergeStatus(i)=3
				Write(*,'(A,I2,A,F8.4,A)') "Isochore",i,":",rho_IC(i)," (g/ml):"
				write(*,'(A19,I3,A)') "Convergence Status= ",convergeStatus(i)," (Converged to third root)"				
				exit
			endif
		endif
		if(j .eq. maxIter)then
			convergeStatus(i)=1
			Write(*,'(A,I2,A,F8.4,A)') "Isochore",i,":",rho_IC(i)," (g/ml):"
			write(*,'(A19,I3,A)') "Convergence Status= ",convergeStatus(i)," (Maximum Iteration s reached)"				
		endif

		tempPsat = Psat(i)
		tempRhoV = rhoV(i)

		if(Psat(i) .lt. 0.0d0 .OR. rhoV(i) .lt. 0.0d0 .OR. Psat(i) .gt. 10.0) then 
			convergeStatus(i)=0
			Write(*,'(A,I2,A,F8.4,A)') "Isochore",i,":",rho_IC(i)," (g/ml):"
			write(*,'(A19,I3,A)') "Convergence Status= ",convergeStatus(i)," (Diverged)"				
		endif
		!>====================================================================

	enddo
		close(53)
		close(7845)
enddo


do i=iFirstIcPt,nItPts
	TsatLnDev(i)=LOG(Tsat(i)/T_IC(i,1))*100.0 
enddo



write(*,*)
write(*,*)'===============================IC Info====================================='
write(*,*)

write(*,'(11(A12,1x))')"(K)","(g/ml)","","",""
write(*,'(11(A12,1x))')"T","rho","Z","aDep","uDep"
do i=iFirstIcPt,nItPts
	do k=1,3

		write(*,'(F12.2,1x,10(F12.6,1x))') T_IC(i,k),rho_IC(i),Z_IC(i,k),&
				aDep_IC(i,k),uDepT_IC(i,k)/T_IC(i,k)
	enddo
enddo
write(*,*)
write(*,*)'===============================IT Info====================================='
write(*,*)
write(*,'(11(A12,1x))')"(K)","(g/ml)","","",""
write(*,'(11(A12,1x))')"T","rho","Z","aDep","uDep","rho"
if(calculateVirialSubcritical) then
	do i=1,4
		write(*,'(F12.2,1x,10(F12.6,1x))') T_IT2_vr(i),rho_IT2_vr(i),Z_IT2_vr(i),0.d0,0.d0,0.d0
	enddo
	write(*,*)
endif

if(calculateVirialSupercritical) then
	do i=1,3
		write(*,'(F12.2,1x,10(F12.6,1x))') T_IT_vr(i),rho_IT_vr(i),Z_IT_vr(i),0.d0,0.d0,0.d0
	enddo
endif

do i=1,nItPts
	write(*,'(F12.2,1x,10(F12.6,1x))') T_IT(i),rho_IT(i),Z_IT(i),&
				aDep_IT(i),uDepT_IT(i)/T_IT(i),rho_IT(i)
enddo

write(*,*)
write(*,*)'===============================Saturation Info====================================='
write(*,*)

write(*,*)'                  (K)      %             (MPa)    (g/ml)      (g/ml)  (KJ/mol)    (ml/g)'
write(*,'(A5,A7,A10,A7,A18,A10,A12,A10,A10,A10,A10)')&
		"Conv","Tr","Tsat","Dev","Psat","rhoL","rhoV","Hvap","B2sat","aDepSat","uDepSat"
do i=iFirstIcPt,nItPts
	write(*,'(I5,F7.3,F10.2,F7.2,F18.12,F10.4,F12.8,F10.4,F10.4,F10.4,F10.4)') &
	convergeStatus(i),Tsat(i)/TC,Tsat(i),TsatLnDev(i),Psat(i),rho_IC(i),rhoV(i),Hvap(i),&
	B2sat(i),aDepSat(i),uDepTSat(i)/Tsat(i)
enddo



end program main
	




subroutine getInterceptSlope(X1,X2,X3,X4,Y1,Y2,Y3,Y4,YINTERCEPT,SLOPE)
	implicit double precision(A-H,O-Z)

	double precision :: xarray(4)
	double precision :: yarray(4)
	
	iCOUNT = 0
	SUMX = 0
	SUMX2 = 0
	SUMY = 0
	SUMXY = 0

	xarray(1) = X1
	xarray(2) = X2
	xarray(3) = X3
	xarray(4) = X4

	yarray(1) = Y1
	yarray(2) = Y2
	yarray(3) = Y3
	yarray(4) = Y4
	
	do i=1,4
		iCOUNT = iCOUNT + 1
		SUMX = SUMX + xarray(i)
		SUMX2 = SUMX2 + xarray(i) ** 2
		SUMY = SUMY + yarray(i)
		SUMXY = SUMXY + xarray(i) * yarray(i)
	enddo

 
	XMEAN = SUMX / iCOUNT
	YMEAN = SUMY / iCOUNT
	SLOPE = (SUMXY - SUMX * YMEAN) / (SUMX2 - SUMX * XMEAN)
	YINTERCEPT = YMEAN - SLOPE * XMEAN

	!B2 = YINTERCEPT
	!B3 = SLOPE
end subroutine getInterceptSlope


double precision function TsatFinder(zLiq,x01,y01,x02,y02,x03,y03)
	implicit none
	integer,parameter::n=3 
	double precision:: zliq,x01,y01,x02,y02,x03,y03,xvalue
	!double precision:: yvalue
	real,dimension (n)::b,x 
	real,dimension(n,n)::a,a1 
	integer::i,j,k,l
	real::z
	a1(1,1)=1.0d0
	a1(1,2)=0.0d0
	a1(1,3)=0.0d0
	a1(2,1)=0.0d0
	a1(2,2)=1.0d0
	a1(2,3)=0.0d0
	a1(3,1)=0.0d0
	a1(3,2)=0.0d0
	a1(3,3)=1.0d0

	a(1,1)=x01**2
	a(1,2)=x01
	a(1,3)=1.d0
	a(2,1)=x02**2
	a(2,2)=x02
	a(2,3)=1.d0
	a(3,1)=x03**2
	a(3,2)=x03
	a(3,3)=1.d0
	b(1)=y01
	b(2)=y02
	b(3)=y03

	!divided all elements of a & a1 by a(i,i) 
	do i=1,n
		z=a(i,i) 
		do j=1,n
			a(i,j)=a(i,j)/z
			a1(i,j)=a1(i,j)/z 
		enddo
		!make zero all entries in column a(j,i) & a1(j,i) 
		do j=i+1,n
		z=a(j,i) 
		do k=1,n
			a(j,k)=a(j,k)-z*a(i,k)
			a1(j,k)=a1(j,k)-z*a1(i,k)
		enddo
		enddo
	enddo
	!subtract appropiate multiple of row j from j-1 
	do i=1,n-1
		do j=i+1,n 
		z=a(i,j) 
			do l=1,n
				a(i,l)=a(i,l)-z*a(j,l) 
				a1(i,l)=a1(i,l)-z*a1(j,l) 
			enddo
		enddo 
	enddo
	do i=1,n 
		do j=1,n 
			x(i)=0 
			do k=1,n
				x(i)=x(i)+a1(i,k)*b(k)
			enddo
		enddo
	enddo
	xvalue=(-1.00*x(2)-SQRT(x(2)**2-4.00*x(1)*(x(3)-zLiq)))/(2.00*(x(1))) 
	TsatFinder=1000.d0/xvalue
end function TsatFinder
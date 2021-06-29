	!Because F90 compiles modules in listed order, the modules defining the FF's must come before the module for selecting a particular FF.

	module TrappeFF
	Implicit DoublePrecision(a-h,o-z)
	parameter (nSiteTypesMax=3, nSitesMax=22)
	dimension iSiteType(nSitesMax)
	dimension iRexpo(nSiteTypesMax),eps_kB(nSiteTypesMax),sigmaNm(nSiteTypesMax),preFactor(22)						   
    data preFactor/0,0,0,0,0,0,0,0,6.75,5.379143536,4.55315466,4,3.603023659,3.303852406,3.070026249,2.882048082,2.727495588,2.598076211,2.48804074,2.393276359,2.310762194,2.238228825/
	!                 CH4,CH3,CH2
	data iRexpo/12,12,12/					
	data eps_kB/148,98,46/					 
	data sigmaNm/0.373d0,0.375d0,0.395d0/
    end module TrappeFF
    
	module TraMieFF
	Implicit DoublePrecision(a-h,o-z)
	parameter (nSiteTypesMax=3, nSitesMax=22)
	dimension iSiteType(nSitesMax)
	dimension iRexpo(nSiteTypesMax),eps_kB(nSiteTypesMax),sigmaNm(nSiteTypesMax),preFactor(22)						  
    data preFactor/0,0,0,0,0,0,0,0,6.75,5.379143536,4.55315466,4,3.603023659,3.303852406,3.070026249,2.882048082,2.727495588,2.598076211,2.48804074,2.393276359,2.310762194,2.238228825/
	!                 CH4,CH3,CH2
	data iRexpo/14,16,16/					 
	data eps_kB/161,121.25d0,61/					 ! 
	data sigmaNm/0.374d0,0.3783d0,0.399d0/
	end module TraMieFF

	module SelectFF
	! Select the desired potential model by uncommenting it, and commenting the other.
	USE TrappeFF !  iSiteType. Set eps,sigma,n__Max, 
	!USE TraMieFF !  iSiteType. Set eps,sigma,n__Max, 
	end module SelectFF

	!Begin main program.
	USE SelectFF !  iSiteType. Set eps,sigma,n__Max, prefactor, ...
	implicit DoublePrecision(a-h,o-z)
	parameter (maxBins=151,maxNgrid=222,zeroTol=1.D-13)
	character*4 siteName,siteNameRd
	character*123 dumString
	LOGICAL LOUD
	dimension rGrid(maxNgrid),rMoleci(nSitesMax,3),rMolecj(nSitesMax,3)
	dimension fMayerAvg(maxNgrid),fMayerAvgSq(maxNgrid),expRef(maxNgrid),expRefSq(maxNgrid),uMatt(maxNgrid),uMattSq(maxNgrid)
	dimension B2(maxBins),B20(maxBins) 
	dimension fMayer(maxBins,maxNgrid),uMapBin(maxBins,maxNgrid),expRefBin(maxBins,maxNgrid),expRefBinSq(maxBins,maxNgrid)
	dimension uMattBin(maxBins,maxNgrid),uMatt2Bin(maxBins,maxNgrid),uMatt3Bin(maxBins,maxNgrid),uMatt4Bin(maxBins,maxNgrid),uMatt5Bin(maxBins,maxNgrid)
	dimension fMayerSq(maxBins,maxNgrid),uMap(maxNgrid),uMapSq(maxNgrid),uMapTmp(maxNgrid),expRefTmp(maxNgrid),uMattOrd(5,maxNgrid) 
	dimension siteName(nSiteTypesMax), B2tmp(5,maxBins),stDevB2tpt(5),B2tpt(9),fMayerCk(maxNgrid),uRefPseudo(maxNgrid)
	common/positions/rMolec(10000,nSitesMax,3)
	data siteName/'CH4','CH3', 'CH2'/
	LOUD=.TRUE.
	LOUD=.FALSE.
	open(51,file='config.txt')
    open(61,file='OFSVC.txt')
    open(62,file='OFSVCsum.txt')
	open(63,file='OneLess.txt')
	read(51,*)nSites,nMolecules,densityG_cc,tKelvin
	write(6  ,'(2i5,F10.5,F10.2,f10.6)')nSites,nMolecules,densityG_cc,tKelvin
	tInv=1/tKelvin  ! this is useful in a lot of the virial coefficient summations.
	nBins=1
	if(nMolecules > 199)nBins=nMolecules*(nMolecules-1)/2/4950	 ! 100 molecules => 4950 interactions is enough to fill one bin.
	if(nBins > maxBins)nBins=maxBins
	Nij_Bin=nMolecules*(nMolecules-1)/2/nBins
	print*,'nBins,Nij/bin=',nBins,Nij_Bin  
	do i=1,nMolecules
		do j=1,nSites
			read(51,'(a123)')dumString
			read(dumString,*,ioStat=ioErr)molec,siteNameRd,( rMolec(i,j,k),k=1,3 )
			if( j  <  nSites-1)write(63,'(i5,a4,1x,3f13.8)')molec,TRIM(siteNameRd),( rMolec(i,j,k),k=1,3 ) 
			if( j == nSites-1)write(63,'(i5,a4,1x,3f13.8)')molec,TRIM( siteName(1) ),( rMolec(i,j,k),k=1,3 ) 
			if(ioErr)write(*,*)'ioErr,dumString=',ioErr,' ',TRIM(dumString)
			do k=1,nSiteTypesMax
				if( TRIM(siteName(k))==TRIM(siteNameRd) )iSiteType(j)=k
			enddo
			do k=1,3
				rMolec(i,j,k)=rMolec(i,j,k)/10 ! convert from Angstroms to nm. 	pdb files use Angstroms for some stupid reason.
			enddo
		enddo
	enddo
	close(51)
	close(63)
	print*,'siteName: ',(siteName( iSiteType(j) ),j=1,nSites)
	write(61,*)'siteName:  ',(   TRIM(  siteName( iSiteType(j) )  ),j=1,nSites   )
	write(61,*)'sigma(nm): ',(       (  sigmaNm( iSiteType(j) )  ),j=1,nSites   )
	write(61,*)'eps/kB(K): ',(       (    eps_kB( iSiteType(j) )  ),j=1,nSites   )
	radGyrationNm=0.373D0+0.154d0*SQRT( DFLOAT(nSites-1) ) ! this is just rough approximation for scale of tabulation.   
	print*,'Radius of Gyration (nm): ',radGyrationNm 
	nShort=10
	do i=1,nShort
		rGrid(i)=i*radGyrationNm/20	   ! r/Rg = [.05,0.5] take big steps expecting <fMayerAvg> ~ -1.  	  
	enddo
	nMiddle=150
	do i=1,nMiddle
		rGrid(i+nShort)=rGrid(nShort)+i*0.005D0*radGyrationNm	! r/Rg = [0.505,1.25]  (0.5+150*0.005 = 1.25)
	enddo
	nLong=15										    ! don't count r = infinity here, so nLong really=16 => even.
	do i=1,nLong
		reciprocal=0.8D0-i*0.05D0				 ! Rg/r = [0.8,0.05] => r/Rg = [1.25,20), but incr=0 at infinity, so we don't need to compute it.
		rGrid(i+nShort+nMiddle)=radGyrationNm/reciprocal
	enddo
	nGrid=nShort+nMiddle+nLong
!	do k=1,nGrid
!		print*,'Main: k,rGrid=',k,rGrid(k)
!	enddo
!	pause 'check rGrid'
    tSim=tKelvin
    tKelvin=INT(tSim*0.55/100)*100
    delTK=.25d0*tSim
    delTK=100
	do iTemperature=1,25
        if(tKelvin >  995)delTK=200
        if(tKelvin > 1995)delTK=2000
        tKelvin=tKelvin+delTK
        if(tKelvin > 10000)exit
	    tInv=1/tKelvin  ! this is useful in a lot of the virial coefficient summations.
	    if(LOUD)print*,'starting Mayer moves. nGrid,rGrid(1)=',nGrid,rGrid(1)
	    fMayerAvg=0										  ! vector init,  Technically, 1+ Mayer function, < 1+f > = < exp(-beta*u) >. This form is more convenient  
	    fMayerAvgSq=0									! vector init
	    expRef=0										    ! < exp(-beta*u0) >, reference averaged Mayer function, W.R. Smith, Can J. Phys. (1974).
	    expRefSq=0								       	  ! vector init
	    uMattOrd=0										  ! vector init for all orders and all k, ~ < -uAtt/kB^ord >
	    nij=0
	    iBin=1
	    do i=1,nMolecules-1
		    do j=i+1,nMolecules
			    nij=nij+1
			    call TranslateToOrigin( i,rMoleci,nSites) ! center to origin before each round of movement.
			    call TranslateToOrigin( j,rMolecj,nSites)
			    call MoveThruRgrid( rMoleci,rMolecj,nSites,rGrid,nGrid,tKelvin,fMayerAvg,expRef,uMattOrd,fMayerAvgSq,expRefSq )
			    call MoveThruRgrid( rMolecj,rMoleci,nSites,rGrid,nGrid,tKelvin,fMayerAvg,expRef,uMattOrd,fMayerAvgSq,expRefSq )
			    !if(LOUD)pause 'check fMayer iMolec,jMolec'
			    if(mod(nij,Nij_bin)==0)then	!start new bin
				    do k=1,nGrid
					    fMayer(iBin,k)=fMayerAvg(k)/(2*Nij_bin)
					    fMayerSq(iBin,k)=fMayerAvgSq(k)/(2*Nij_bin)
					    fMayerAvg(k)=0
					    fMayerAvgSq(k)=0
					    expRefBin(iBin,k)=expRef(k)/(2*Nij_bin)
					    expRefBinSq(iBin,k)=expRefSq(k)/(2*Nij_bin)
					    expRef(k)=0
					    expRefSq(k)=0
					    uMattBin(iBin,k)  =uMattOrd(1,k)/(2*Nij_bin)
					    uMatt2Bin(iBin,k)=uMattOrd(2,k)/(2*Nij_bin)
					    uMatt3Bin(iBin,k)=uMattOrd(3,k)/(2*Nij_bin)
					    uMatt4Bin(iBin,k)=uMattOrd(4,k)/(2*Nij_bin)
					    uMatt5Bin(iBin,k)=uMattOrd(5,k)/(2*Nij_bin)
					    fMayerCk(k) = expRefBin(iBin,k)+tinv*(uMattBin(iBin,k)+tinv/2*(uMatt2Bin(iBin,k)+tInv/3*(uMatt3Bin(iBin,k)+tInv/4*(uMatt4Bin(iBin,k)+tInv/5*uMatt5Bin(iBin,k)))))
					    uRef=0
					    if(expRefBin(iBin,k) > 0)uRef=LOG( expRefBin(iBin,k) )
					    !print*,' nij,r,fMayer=',nij,rGrid(k),fMayer(iBin,k),uMattBin(iBin,k),uRef
					    !write(*,'( a,i4,f8.4,3(1x,F12.2) )' )  ' nij,r,fMayer=',nij,rGrid(k),fMayer(iBin,k),uMattBin(iBin,k),uRef
				    enddo
				    uMattOrd=0 ! reset all orders to zero for all k.  
				    nij=0 ! restart counting
				    print*,'iBin=',iBin
				    !pause 'Check uMatt'
				    iBin=iBin+1
				    if(iBin > nBins) exit ! exits just this do loop.
			    endif ! new bin started
		    enddo ! j = molecules
		    if(iBin > nBins) exit !  to exit both do loops.
	    enddo  ! i = molecules
	    if(LOUD)print*,'r/Rg, fM,fMck,eRef,<uP/T>,<uP^2>,<uP^3>,<uP^4>'
	    do k=1,nGrid
		    if(LOUD)write(*,'(f7.3,4f7.4,3f7.2)')rGrid(k),fMayer(1,k),fMayerCk(k),expRefBin(1,k),uMattBin(1,k)/tKelvin,uMatt2Bin(1,k),uMatt3Bin(1,k),uMatt4Bin(1,k)
	    enddo
	    print*,'nBins,nij leftover=',nBins,nij
	    !pause 'check fMayers, nBins'

	    write(6  ,'(2i5,F10.5,F10.2,f10.6)')nSites,nMolecules,densityG_cc,tKelvin,radGyrationNm
	    write(61,'(2i5,F10.5,F10.2,f10.6)')nSites,nMolecules,densityG_cc,tKelvin,radGyrationNm
	    fMayerAvg=0
	    fMayerAvgSq=0
	    uMap=0					! uMap(r) [=] K = -tKelvin*log( < exp(-uism(r)/kT) -1 > );  aka. "Mayer averaged potential," equivalent to low density pmf.
	    uMapSq=0
	    expRef=0
	    expRefSq=0
	    uMatt=0
	    uMattSq=0
	    do iBin=1,nBins
		    do k=1,nGrid
			    !print*,'nBins,iBin,k=',nBins,iBin,k
                if(fMayer(iBin,k) < 1.D-175)fMayer(iBin,k)=1.D-175
			    uMapNow= -tKelvin*LOG( fMayer(iBin,k) )
			    uMapBin(iBin,k)=uMapNow	
			    uMap(k)=uMap(k)+uMapNow/nBins								!sum(xi-mean)^2 = sum( xi^2 - 2xi*mean + mean^2) = sum(xi^2) -2*mean*sum(xi) + mean^2 = sum(xi^2) -2*N*mean^2 + mean^2
			    uMapSq(k)=uMapSq(k)+uMapNow*uMapNow/nBins
			    !print*,'iBin,k=',iBin,k
			    fMayer(iBin,k)=fMayer(iBin,k)-1 ! preparing for B2 calculation
			    fMayerAvg(k)=fMayerAvg(k)+fMayer(iBin,k)/nBins									!sum(xi-mean)^2 = sum( xi^2 - 2xi*mean + mean^2) = sum(xi^2) -2*mean*sum(xi) + mean^2 = sum(xi^2) -2*N*mean^2 + mean^2
			    fMayerAvgSq(k)=fMayerAvgSq(k)+fMayerSq(iBin,k)/nBins
			    expRef(k)=expRef(k)+expRefBin(iBin,k)/nBins									!sum(xi-mean)^2 = sum( xi^2 - 2xi*mean + mean^2) = sum(xi^2) -2*mean*sum(xi) + mean^2 = sum(xi^2) -2*N*mean^2 + mean^2
			    expRefSq(k)=expRefSq(k)+expRefBinSq(iBin,k)/nBins
			    uMatt(k)=uMatt(k)+uMattBin(iBin,k)/nBins									!sum(xi-mean)^2 = sum( xi^2 - 2xi*mean + mean^2) = sum(xi^2) -2*mean*sum(xi) + mean^2 = sum(xi^2) -2*N*mean^2 + mean^2
			    uMattSq(k)=uMattSq(k)+uMatt2Bin(iBin,k)/nBins
		    enddo ! nGrid
	    enddo  ! nBins
	    uMin=1234		! We must find the best average before seeking the best uMin.
	    do k=1,nGrid
    !		uMatt(k)= uMap(k)	 ! for now, just store everything and find the minimum.
    !		uMattSq(k)=uMapSq(k)
		    if(uMap(k) < uMin)then
			    uMin=uMap(k)
			    kMin=k
			    rMin=rGrid(k)
		    endif
	    enddo ! nGrid
	    uRefPseudo=uMap(1)-uMin	 ! Set first value here 
	    dEhs = 0
	    sigma=0
	    do k=2,kMin
    !		uMatt(k) = uMin
		    if(uMap(k)*uMap(k-1) < 0)sigma=rGrid(k-1)+( 0-uMap(k-1) )/( uMap(k)-uMap(k-1) )*( rGrid(k)-rGrid(k-1) )
		    uRefPseudo(k) = uMap(k)-uMin	   ! This result, combined with uMap, is ready for application as a pseudopotential
		    dEhs=dEhs+( exp(-uRefPseudo(k)/tKelvin)+exp(-uRefPseudo(k-1)/tKelvin) )/2*( rGrid(k)-rGrid(k-1) )																		
    !		uMattSq(k)=uMin*uMin
	    enddo
	    dEhs=rMin-dEhs
	    write(6  ,'(a, i3, 4f10.4,4f8.3)') ' kMin,rMin,uMin,sigmaMap,dEhs:',kMin,rMin,uMin,sigma,dEhs
	    write(61,'(a, i3, 4f10.4,4f8.3)') ' kMin,rMin,uMin,sigma,dEhs:',kMin,rMin,uMin,sigma,dEhs
	    !pause 'check uMin etc'

	    write(*  ,'(a,i3,f7.4,6f10.5)')' k,rGrid,fMayerAvg,uMap,expRef,uMatt:'
	    write(61,'(a,i3,f7.4,6f10.5)')'  k    rGrid fMayerAvg      +/-      uMap        +/-         expRef       +/-     e0*uMatt      +/-:'
	    do k=1,nGrid
		    variance=( fMayerAvgSq(k)-fMayerAvg(k)*fMayerAvg(k) )
		    if( variance < 0) variance=1.D-11
		    stDevF=SQRT( variance )
		    variance=( uMapSq(k)-uMap(k)*uMap(k) )
		    if( variance < 0) variance=1.D-11
		    stDevU=SQRT( variance )
		    variance=( uMattSq(k)-uMatt(k)*uMatt(k) )
		    if( variance < 0) variance=1.D-11
		    stDevA=SQRT( variance )
		    variance=( expRefSq(k)-expRef(k)*expRef(k) )
		    if( variance < 0) variance=1.D-11
		    stDevE=SQRT( variance )
		    write(*  ,'( i5,f7.3,f10.5,6(1PE12.4) )')k,rGrid(k),fMayerAvg(k),uMap(k),expRef(k),uMatt(k)
		    write(61,'(i5,f7.3,2f10.4,6(1PE12.4) )')k,rGrid(k),fMayerAvg(k),stDevF,uMap(k),stDevU,expRef(k),stDevE,uMatt(k),stDevA
	    enddo
	    if(LOUD)pause 'Check averaged fMayer & uMap'
	    !Perform integration: Simpson's rule for now.
	    write(61,'(a )') ' k,rGrid(k),fMayerAvg(k),stDevF,uMap(k),stDevU,expRefTmp(k),stDevE,uMattTmp(k),stDevA	 '
	    do iBin=1,nBins
		    do k=1,nGrid
			    uMapTmp(k)=uMapBin(iBin,k)
			    expRefTmp(k)=expRefBin(iBin,k)
			    uMattOrd(1,k)=uMattBin(iBin,k)		 ! re-use uMatts
			    uMattOrd(2,k)=uMatt2Bin(iBin,k)
			    uMattOrd(3,k)=uMatt3Bin(iBin,k)
			    uMattOrd(4,k)=uMatt4Bin(iBin,k)
			    uMattOrd(5,k)=uMatt5Bin(iBin,k)
		    enddo
		    call integrateB2(rGrid,nShort,nMiddle,nLong,uMapTmp,expRefTmp,uMattOrd,tKelvin,B2(iBin),B20(iBin),B2tmp(1,iBin) ) ! B2tmp(1,iBin) is a pointer to the first order at iBin.  iBin must be the 2nd subscript for this to work.
		    !write(6,'(a,i3,3F10.7)')' iBin,B2,B2sq ',iBin,B2(iBin),B2sq(iBin)
	    enddo ! iBin
	    call B2average(nBins,B2  ,B2avg  ,stDevB2)
	    call B2average(nBins,B20,B20avg,stDevB20)
	    do j=1,5
		    B2(1:nBins)=B2tmp(j,1:nBins)  ! reuse B2() vector
		    call B2average(nBins,B2,B2tpt(j),stDevB2tpt(j) )  ! reuse B2tpt() for final results.
	    enddo
	    write(6 ,'(a,4F10.3)') ' T(K)     B2(cc/mol)    B20     B21/T    B22/T^2     B23/T^3      B24/T^4      B25/T^5'  
	    write(6 ,'( f7.2,2F10.3,5(f12.4) )') tKelvin, B2avg, B20avg,(B2tpt(j),j=1,5)
	    write(6 ,'( f7.2,2F10.3,5(f12.4) )') tKelvin, stDevB2, stDevB20,(stDevB2tpt(j),  j=1,5)
	    write(61,'(a,4F10.3)') ' T(K)     B2(cc/mol)    B20     B21/T    B22/T^2     B23/T^3      B24/T^4      B25/T^5'  
	    write(61,'( f7.2,2F10.3,5(f12.4) )') tKelvin, B2avg, B20avg,(B2tpt(j),j=1,5)
	    write(61,'( f7.2,2F10.3,5(f12.4) )') tKelvin, stDevB2, stDevB20,(stDevB2tpt(j),  j=1,5)
	    write(62,'( f7.2,2F10.3,5(f12.4) )') tKelvin, B2avg, B20avg,(B2tpt(j),j=1,5)
	    write(62,'( f7.2,2F10.3,5(f12.4) )') tKelvin, stDevB2, stDevB20,(stDevB2tpt(j),  j=1,5)
	    !write(6 ,'(a,4F10.3)') ' T(K),B2(cc/mol),+/-= ',tKelvin,B2avg,stDevB2
	    !write(61,'(a,4F10.3)') ' T(K),B2(cc/mol),+/-= ',tKelvin,B2avg,stDevB2
	    !pause 'check B2'
        iniTemp=200
        if(nSites==2)iniTemp=100 ! go low for ethane.
	    do kelvins=iniTemp,1500,50
		    beta=tKelvin/kelvins !relative to refTemp
		    B2tpt(6)=B2tpt(5)*(  B2tpt(5)/( B2tpt(4) )  )
		    B2tpt(7)=B2tpt(6)*B2tpt(5)/B2tpt(4)
		    B2tot=B20avg+beta*(B2tpt(1)+beta*(B2tpt(2)+beta*(B2tpt(3)+beta*(B2tpt(4)+beta*(B2tpt(5)+beta*(B2tpt(6)+beta*B2tpt(7)))))))
		    write(6  ,'(a,i5,f11.3)' ) '  T(K), B2tot(cc/mol)= ' , kelvins, B2tot
		    write(61,'(a,i5,f11.3)' ) '  T(K), B2tot(cc/mol)= '	, kelvins, B2tot
	    enddo
	    write(61,*)'       rGrid   exp(-beta*uRefPseudo)  uRefPseudo '
    enddo ! iTemperature    
	close(62)
    open(62,file='OFSVCsum.txt') !read/copy summary to end of long file.
    do i=1,111
        read(62,*,ioStat=ioErr)tKelvin, B2avg, B20avg,(B2tpt(j),j=1,5)
        if(ioErr)exit
        read(62,*,ioStat=ioErr)tKelvin, stDevB2, stDevB20,(stDevB2tpt(j),  j=1,5)
	    write(61,'( f7.2,2F10.3,5(f12.4) )') tKelvin, B2avg, B20avg,(B2tpt(j),j=1,5)
	    write(61,'( f7.2,2F10.3,5(f12.4) )') tKelvin, stDevB2, stDevB20,(stDevB2tpt(j),  j=1,5)
    enddo
    close(62)
	close(61)
	stop
	end

	Subroutine MoveThruRgrid( rMoleci,rMolecj,nSites,rGrid,nGrid,tKelvin,fMayerAvg,expRef,uMattOrd,fMayerAvgSq,expRefSq )
	USE SelectFF !  iSiteType. Set eps,sigma,n__Max, 
	implicit DoublePrecision(a-h,o-z)
	parameter (maxNgrid=222)
	!LOGICAL LOUD
	dimension rGrid(maxNgrid),unitVector(3),rMove(nSitesMax,3),rMoleci(nSitesMax,3),rMolecj(nSitesMax,3)
	dimension fMayerAvg(maxNgrid),fMayerAvgSq(maxNgrid),expRef(maxNgrid),expRefSq(maxNgrid)
	dimension uMattOrd(5,maxNgrid) ! uMatt(maxNgrid),uMattSq(maxNgrid)
	!dimension uMattTmp(maxNgrid)
	do m=1,3  ! move moleci	in 3 dimensions, holding molecj at the origin.
		unitVector=0
		unitVector(m)=1
		do k=1,nGrid  
			rMove=rMoleci  ! set rMove initially. This should be centered to origin. Then move the mth direction.  
			do iSite=1,nSites
				rMove(iSite,m)=rMoleci(iSite,m)+unitVector(m)*rGrid(k) !just move in the mth direction.
				!if(LOUD)print*,'rMove=',(rMove(iSite,n),n=1,3)
			enddo
			!print*,'calling Mayer. nSites,rGrid=',nSites,rGrid(k)
			call ismPot( rMove,rMolecj,nSites,uism_kB,uRef_kB,uAtt_kB, iErr )
			!write(*,'(a,f8.4,3(1x,F12.2) )')' Move...: r,u,uRef,uAtt:',rGrid(k),uism_kB,uRef_kB,uAtt_kB
			if( iErr .ne. 0)pause 'Main: iErr.ne.0 from ismPot'
			betaU=uism_kB/tKelvin
			fMayerNow=EXP( -uism_kB/tKelvin )
			if( fMayerNow < 1.D-33 ) fMayNow=1.D-33  ! dunno how this could happen, but we take the log later so we need to be sure. 
			fMayerAvg(k)=fMayerAvg(k)+fMayerNow/3
			fMayerAvgSq(k)=fMayerAvgSq(k)+fMayerNow*fMayerNow/3
			eRef=EXP( -uRef_kB/tKelvin )
            if(eRef > 1.D33)eRef=1.D33
			expRef(k)=expRef(k)+eRef/3
			expRefSq(k)=expRefSq(k)+eRef*eRef/3
			uMatt_kB=    ( -uAtt_kB ) 				 ! fMayer ~ exp( -u ) so uMatt ~ -uAtt.  Note: uMatt > 0 always. 
			uMattOrd(1,k)=uMattOrd(1,k)+eRef*uMatt_kB/3
			uMattOrd(2,k)=uMattOrd(2,k)+eRef*uMatt_kB/3*uMatt_kB
			uMattOrd(3,k)=uMattOrd(3,k)+eRef*uMatt_kB/3*uMatt_kB*uMatt_kB
			uMattOrd(4,k)=uMattOrd(4,k)+eRef*uMatt_kB/3*uMatt_kB*uMatt_kB*uMatt_kB
			uMattOrd(5,k)=uMattOrd(5,k)+eRef*uMatt_kB/3*uMatt_kB*uMatt_kB*uMatt_kB*uMatt_kB
			!if(LOUD.and.rGrid(k) < radGyrationNm*1.5)write(*,'(a,2f10.5,E11.4)')' Main: rGrid,fMayer=',rGrid(k),fMayerNow,uism_kB
        enddo
        continue
		!if(LOUD)pause 'check betaU'
		!pause 'check uAtt'
	enddo ! move moleci done
	return
	end

	Subroutine B2average(nBins,B2,B2avg,stDev)
	implicit DoublePrecision(A-H,O-Z)
	!logical LOUD
	dimension B2(nBins) !uMap(maxBins,maxNgrid) !use pointer here so we can receive any particular iBin.
	B2avg=0
	B2avgSq=0
	do i=1,nBins
		B2avg=B2avg+B2(i)
		B2avgSq=B2avgSq+B2(i)*B2(i)
	enddo
	B2avg  =  B2avg/nBins
	B2avgSq=B2avgSq/nBins
	variance=B2avgSq-B2avg*B2avg
	stDev = 0	
	if(variance > 0)stDev=SQRT( variance )
	return
	end

	Subroutine integrateB2(rGrid,nShort,nMiddle,nLong,uMap,expRef,uMattOrd,tKelvin,B2,B20,B2tpt  )
	implicit DoublePrecision(A-H,O-Z)
	logical LOUD
	dimension rGrid(nshort+nMiddle+nLong),uMap(nshort+nMiddle+nLong),expRef(nshort+nMiddle+nLong) !,uMatt(nshort+nMiddle+nLong) !uMap(maxBins,maxNgrid) !use pointer here so we can receive any particular iBin.
	dimension fMayer(nshort+nMiddle+nLong)
	dimension B2tpt(5),B2tmp(5),uMattOrd(5,nshort+nMiddle+nLong)
	LOUD=.FALSE.
	!LOUD=.TRUE.
	pi=3.14159265D0
	avoNum=602.22D0
	nGrid=nShort+nMiddle+nLong
	if(LOUD)write(*,'(11f7.4)')(rGrid(i),i=1,11)
	if(LOUD)write(*,'(7E11.4)')(uMap(i)/tKelvin,i=1,7)
	if(LOUD)print*,'integrateB2: nShort,nMiddle,nLong=',nShort,nMiddle,nLong
	B2short=0 !first step,  fMayer*r^2 = 0 at rGrid = 0
	B20 = 0.
	B2tpt=0	  !vector init
	do i=1,nShort-1
		fMayer(i)=EXP( -uMap(i)/tKelvin ) -1
		iFactor=2+2*mod(i,2)  ! = 4 if i=odd, 2 if i=even.
		B2short=B2short+iFactor*rGrid(i)*rGrid(i)*fMayer(i) 
		B20     =B20      +iFactor*rGrid(i)*rGrid(i)*(expRef(i)-1)
		do j=1,5
			uTptj=uMattOrd(j,i)
			B2tpt(j)=B2tpt(j)+iFactor*rGrid(i)*rGrid(i)*uTptj ! NOTE: factorial in denom is accounted for at the end.
		enddo
		!write(6  ,'(a,2i4, 2f10.5)')' iBin,k,fMayer,incr',iBin,i,fMayer ,B2increment
		!if(LOUD)write(6,'(a,i4, 2f10.5)')' k,fMayer,incr',i,fMayer(i) 
	enddo  ! nShort must be even so we start even (zero=even) and end even and odds will be iFactor=4.  e.g. 0,1,2 => f(x0)+4*f(x1)+f(x2)
	B2short=( B2short+( EXP( -uMap(nShort)/tKelvin )-1 )*rGrid(nShort)**2 )*rGrid(1)/3	!rGrid(1) = dr
	B20stor     =( B20      + rGrid(nShort)*rGrid(nShort)*(expRef(nShort)-1) )*rGrid(1)/3
	do j=1,5
		uTptj= uMattOrd(j,nShort)	  ! NOTE: fMayer=exp( -u ) so B2tpt ~ -u
		stepIntegrand=rGrid(nShort)*rGrid(nShort)*uTptj
		B2tmp(j)=( B2tpt(j)+stepIntegrand )*rGrid(1)/3
		B2tpt(j) = stepIntegrand	! initiating the middle integration.
	enddo
	if(LOUD)print*,'B2short,B20short=',B2short,B20stor


	B2middle= ( EXP( -uMap(nShort)/tKelvin )-1 )*rGrid(nShort)**2
	B20        = ( expRef(nShort)-1 )*rGrid(nShort)**2
	dRmiddle= rGrid(nShort+2)-rGrid(nShort+1)
	do i=nShort+1,nShort+nMiddle-1	 ! NOTE: stopping one short of nMiddle
		fMayer(i)=EXP( -uMap(i)/tKelvin ) -1
		iFactor=2+2*mod(i,2)  ! = 4 if i=odd, 2 if i=even.
		rGrid2 = rGrid(i)*rGrid(i)
		B2middle=B2middle+iFactor*( fMayer(i) *rGrid2 )
		B20=B20+iFactor*( expRef(i)-1 ) *rGrid2 
		uTptj=1 
		do j=1,5
			uTptj= uMattOrd(j,i)	  ! NOTE: fMayer=exp( -u ) so B2tpt ~ -u
			B2tpt(j)=B2tpt(j)+iFactor*uTptj *rGrid2 
		enddo
		!write(*,'(a,4f15.5)') ' r,B21intgrd',rGrid(i),expRef(i)*uMatt(i)*rGrid2,B21,B2tpt(1)
		!if(LOUD)write(6,'(a,2i4, 2f10.5)')' k,iFactor,fMayer ',i,iFactor,fMayer(i) 
	enddo  ! nMiddle must be even so we start and end even and odds will be iFactor=4
	rGrid2=rGrid(nShort+nMiddle)*rGrid(nShort+nMiddle)
	B2middle=( B2middle+( EXP( -uMap(nShort+nMiddle)/tKelvin )-1 )*rGrid2 )*( dRmiddle )/3
	B20        =( B20        +( expRef(nShort+nMiddle)-1 )*rGrid2 )*( dRmiddle )/3
	B20stor=B20stor+B20
	uTptj=1 
	do j=1,5
		uTptj= uMattOrd(j,nShort+nMiddle)	  ! NOTE: fMayer=exp( -u ) so B2tpt ~ -u
        stepIntegrand=uTptj*rGrid2
		B2tmp(j)=B2tmp(j)+( B2tpt(j)+stepIntegrand )*( dRmiddle )/ 3
		B2tpt(j) = stepIntegrand*rGrid2	! initiating the long integration.  Note: rGrid^2 already included in stepIntegrand so just factor rGrid^2 more.
	enddo
	if(LOUD)print*,'B2middle,B20mid=',B2middle,B20


	B2Long= ( EXP( -uMap(nShort+nMiddle)/tKelvin )-1 )*rGrid(nShort+nMiddle)**4
	B20     = ( expRef(nShort+nMiddle)-1 )*rGrid(nShort+nMiddle)**4
	do i=nShort+nMiddle+1,nGrid
		fMayer(i)=EXP( -uMap(i)/tKelvin ) -1
		iFactor=2+2*mod(i,2)  ! = 4 if i=odd, 2 if i=even.
		rGrid4 = rGrid(i)*rGrid(i)*rGrid(i)*rGrid(i)
		B2Long=B2Long+iFactor*( fMayer(i) *rGrid4 )
		B20     =B20      +iFactor*( expRef(i)-1 ) *rGrid4 
		uTptj=1 
		do j=1,5
			uTptj= uMattOrd(j,i)
			B2tpt(j)=B2tpt(j)+iFactor*uTptj*rGrid4
		enddo
		!if(LOUD)write(6,'(a,2i4, 2f10.5)')' k,iFactor,fMayer ',i,iFactor,fMayer(i) 
	enddo
	reciprocalDelta= -( 1/rGrid(nGrid)-1/rGrid(nGrid-1) )
	B2Long= ( B2Long+ 0 )*reciprocalDelta/3
	B20= ( B20+ 0 )*reciprocalDelta/3
	if(LOUD)print*,'B2Long=',B2Long


	B2  = -2*pi*avoNum*(B2short+B2middle+B2Long)
	B20= -2*pi*avoNum*(B20stor+B20)
	factorial=1
	tFactor=1
	do j=1,5
		factorial=factorial*j
		tFactor=tFactor*tKelvin
		B2tpt(j)= -2*pi*avoNum*( B2tmp(j)+B2tpt(j)*reciprocalDelta/3 ) / ( factorial*tFactor )	! combining 4 steps here to minimize looping: reciprocal data, adding , factorial, Tscale
	enddo
	!print*,'B21=',B2tpt(1)
	!pause 'integrateB2: check B21'
	return
	end

	Subroutine ismPot( rMoleci,rMolecj,nSites, UijMolec, uRef_kB,uAtt_kB, iErr )
	USE SelectFF !  iSiteType. Set eps,sigma,n__Max, 
	implicit DoublePrecision(a-h,o-z)
	LOGICAL LOUD
	dimension rMoleci(nSitesMax,3),rMolecj(nSitesMax,3)
	!data    rMin/1.122462048D0/
	data rMinSq/1.25992105D0/
	!dimension centeri(3),centerj(3),rTranslate(3),unitVector(3)  ! work arrays.
	LOUD=.FALSE.
	!LOUD=.TRUE.
	iErr=0
	UijMolec=0
	uRef_kB=0
	uAtt_kB=0
	if(LOUD)print*,'ismPot: '
	do i=1,nSites
		do j=1,nSites
			epsij=SQRT(  eps_kB( iSiteType(i) )*eps_kB( iSiteType( j ) )  )
			sigij = (  sigmaNm( iSiteType(i) )+sigmaNm( iSiteType( j ) )  )/2
            iRexpoij=(  iRexpo( iSiteType(i) )+iRexpo( iSiteType( j ) )  )/2
			if(LOUD)print*,'ismPot: epsij,sigij',epsij,sigij
			rijSq=0
			do k=1,3
				rijSq=rijSq+( rMoleci(i,k)-rMolecj( j,k) )*( rMoleci(i,k)-rMolecj( j,k) ) !
			enddo
			sij_rijSq=sigij*sigij/rijSq
			if(LOUD)print*,'rijSq,rijSq/sigSq=',rijSq,1/sij_rijSq
			rInv=SQRT(sij_rijSq)
            rInv2=sij_rijSq
            rInv6=rInv2*rInv2*rInv2
			sqrtUij= 2*rInv6 - 1	 ! = 2(sigma/rij)^6 - 1
			if(sqrtUij > 3E25)then ! 1E307 ~ 3E25^12 is largest number in Excel.
				Uij_kB= 1E33	 ! exp(-inf) = 0.
				cycle
			elseif(iRexpoij==12)then
				Uij_kB=( sqrtUij*sqrtUij - 1 )
            elseif(iRexpoij==13)then
                rInv13=rInv*rInv6*rInv6
                Uij_kB=preFactor(13)*( rInv13 - rInv6 ) !
            elseif(iRexpoij==14)then
                rInv14=rInv2*rInv6*rInv6
			    Uij_kB=preFactor(14)*( rInv14 - rInv6 ) !
            elseif(iRexpoij==15)then
                rInv15=rInv*rInv2*rInv6*rInv6
			    Uij_kB=preFactor(15)*( rInv15 - rInv6 ) !
            elseif(iRexpoij==16)then
                rInv16=rInv2*rInv2*rInv6*rInv6
			    Uij_kB=preFactor(16)*( rInv16 - rInv6 ) !
            else
                pause 'uism: unsupported iRexpo'
			endif
			UijMolec=UijMolec+Uij_kB*epsij
			if(1/sij_rijSq .le. rMinSq)then
				uRef_kB=uRef_kB+(Uij_kB+1)*epsij
				uAtt_kB=uAtt_kB - epsij
			else
				uAtt_kB=uAtt_kB + Uij_kB*epsij	
			endif
			if(LOUD)print*,'ismPot: UijMolec,uRef,uAtt',Uij_kB*epsij,(Uij_kB+1)*epsij, -epsij
		enddo !j=1,nSites
	enddo ! i=1,nSites
	!uAtt_kB = UijMolec-uRef_kB
	if(LOUD)print*,'UijMolec=',UijMolec
	if(LOUD)pause 'check Uij'
	return
	end 

	Subroutine TranslateToOrigin( jMolec,rMolecj, nSites )
	USE SelectFF !  iSiteType. Set eps,sigma,n__Max, 
	implicit DoublePrecision(a-h,o-z)
	LOGICAL LOUD
	dimension centerj(3),rMolecj(nSitesMax,3)  ! work arrays.
	common/positions/rMolec(10000,nSitesMax,3)
	LOUD=.FALSE.
	!LOUD=.TRUE.
	centerj=0
	do i=1,nSites
		do k=1,3
			centerj(k)=centerj(k)+rMolec(jMolec,i,k)/nSites
		enddo
	enddo
	if(LOUD)write(*,'(a,3f10.4)')' centerj=',(centerj(k),k=1,3)
	do i=1,nSites
		do k=1,3
			rMolecj(i,k)=rMolec( jMolec,i,k )-centerj(k)
		enddo
		if(LOUD)write(*,'(a,i3,3f10.4)')' iSite,rMolecj=',i,(rMolecj(i,k),k=1,3)
	enddo
	return
	end
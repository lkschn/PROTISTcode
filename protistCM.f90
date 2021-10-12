	! 6 char name for process mathc with second line of PDF  
subroutine PROTCM     ( pmsa   , fl     , ipoint , increm, noseg , &                            
							noflux , iexpnt , iknmrk , noq1  , noq2  , &                            
							noq3   , noq4   )     
! dont touch next line, replace name though
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS: 'PROTCM' :: PROTCM                                     
!                                                                                                     
!*******************************************************************************                      
!  
use protist_math_functions
use protist_cell_functions
use protist_uptake_functions
use protist_photosynthesis_functions
use protist_phagotrophy_functions

	IMPLICIT NONE                                                                                   
!                                                                                                     
!     Type    Name         I/O Description                                                            
!          
    integer, parameter :: plen = 216 ! total length of the PMSA input and output array
	real(4) pmsa(*)      ! I/O Process Manager System Array, window of routine to process library     
	real(4) fl(*)        ! O  Array of fluxes made by this process in mass/volume/time               
	integer ipoint(plen) ! I  Array of pointers in pmsa to get and store the data                    
	integer increm(plen) ! I  Increments in ipoint for segment loop, 0=constant, 1=spatially varying 
	integer noseg        ! I  Number of computational elements in the whole model schematisation     
	integer noflux       ! I  Number of fluxes, increment in the fl array                            
	integer iexpnt(4,*)  ! I  From, To, From-1 and To+1 segment numbers of the exchange surfaces     
	integer iknmrk(*)    ! I  Active-Inactive, Surface-water-bottom, see manual for use              
	integer noq1         ! I  Nr of exchanges in 1st direction (the horizontal dir if irregular mesh)
	integer noq2         ! I  Nr of exchanges in 2nd direction, noq1+noq2 gives hor. dir. reg. grid  
	integer noq3         ! I  Nr of exchanges in 3rd direction, vertical direction, pos. downward    
	integer noq4         ! I  Nr of exchanges in the bottom (bottom layers, specialist use only)     
!                                                                                                     
!*******************************************************************************                      
!                                                                                                     
!     Type    Name         I/O Description                                        Unit                
!                                                                                                     
!     support variables
	integer ipnt(plen)    ! Local work array for the pointering                                    
	integer iseg          ! Local loop counter for computational element loop                      
	integer ioq
	integer iflux   
	integer ikmrk1        ! first segment attribute
    	  
	integer ispec         ! local species number counter
	integer spInc         ! local species PMSA number counter
	integer inpItems      ! nr of input items need for output PMSA
	
	 !input parameters
     integer    maxNrSp, nrSp, nrSpCon, nrSpInd     ! constant and species numbers   
     integer    maxNrPr                             ! maxNrPrey
     real    relPhag                                     ! relative phagotrophy night:day
     real    UmRT, Q10, RT, CR                           ! growth and respiration rate calculation 
     real    NCm, NO3Cm, PCm, ChlCm                      ! maximum NC, PC, ChlC quotas
     real    NCo, PCo, ChlCo                             ! minimum NC and PC quotas
     real    NCopt, NO3Copt, PCopt                       ! optimal NC and PC quotas
     real    KtP, KtNH4, KtNO3                           ! half saturation constants
     real    PCoNCopt, PCoNCm                            ! P status influence on optimum NC
     real    ReUmNH4, ReUmNO3, redco, PSDOC, maxPSreq, relPS     ! relative growth rates with specific nutrients  
     real    CcellProt, rProt                            ! parameters for protozooplankton cell
     real    optCR                                       ! parameters for encounter
     real    kAE, AEm, AEo                               ! parameters for assimilation efficiency
     real    SDA                                         ! specific dynamic action
     real    MrtRT, FrAut, FrDet                         ! reference mortality and fractions
     real    alpha                                       ! inital slope
 
     ! input state variables
	 real    protC, protChl, protN, protP                ! protist state variables
	 real    PO4, NH4, NO3                               ! nutrient state variables
	 real    Temp                                        ! physical abiotic variables 
     real    PFD, atten, exat                            ! available light and extinction

     ! prey state variables
     real, dimension(:), allocatable :: vec_preyC, vec_preyChl, vec_preyN, vec_preyP, vec_preySi 
     ! other prey input parameters
     real, dimension(:), allocatable :: vec_CcellPrey, vec_rPrey, vec_motPrey, vec_PR
		 
	 ! auxiliaries
     real    lightInh                                    ! inhibtion of feeding in dark 
     real    NC, PC, ChlC                                ! cell nutrient quotas
     real    UmT, BR                                     ! growth and repsiration rates
     real    NCu, PCu, NPCu                              ! nutrient status within the cell
     real    motProt                                     ! motility 
     real    upP, upNH4, upNO3                          ! nutrient uptake
     real    PSqm, Cfix, synChl, degChl                 ! plateau and Cifx through photosynthesis
     real    CfixPS                                     ! C fix minus phototsynthesis related respiration 
     real    PS                                          ! req for C to come from PS 
     ! food quantity
     real    sumCP        ! total captured prey
     real    ingNC, ingPC ! total ingested N and P 
     real    preyFlag     ! sum of preyFlag (can be 0 = both low, 1 = 1 ok, 2 = both ok)
     real, dimension(:), allocatable :: vec_preyFlag                       ! protection aginst small prey conc. 
     real, dimension(:), allocatable :: vec_nrPrey                         ! prey abundance
     real, dimension(:), allocatable :: vec_smallerVel, vec_largerVel      ! velocities
     real, dimension(:), allocatable :: vec_encPrey                        ! prey encounter
     real, dimension(:), allocatable :: vec_capturedPrey                   ! prey capture
     real, dimension(:), allocatable :: vec_propPrey, vec_ingNC, vec_ingPC ! preyN and preyP proprtion in diet
     ! food quality
     real    stoichP, ppNC, ppPC                    ! stoichiometry comparison                 
     real    opAE                                   ! assimilation efficiency
     ! ingestion and assimilation 
     real    reqPred                                ! required Predation
     real    maxIng, ingSat, ingC, ingN, ingP, KI   ! ingestion
     real    assC, assN, assP                       ! assimilation 
     ! respiration, Cu and mortality
     real    totR, Cu                                    ! respiration and C-growth
     real    mrt, mrtFrAut, mrtFrDet                     ! mortality to detritus and autolysis
     
     ! other parameters
     real, parameter :: wTurb = 0.0 ! this needs to be an input!!!!
     real(8),  parameter :: PI_8  = 4 * atan (1.0_8) 
     real, parameter :: threshNut = 1.0E-10 ! protection against small nut conc.
     
     ! loop counter 
     integer iPrey      ! counter for loops

	 ! Fluxes
     real    dNH4up, dNO3up, dPup                        ! uptake fluxes
     real    dCfix                                       ! photosynthesis flux
     real    dChlsyn, dChldeg                            ! Chl synthesis  and degradation flux
     real    dCresp                                      ! respiration flux
     real    dDOCleak                                    ! C leak through photosynthesis 
     real    dDOCvoid, dNH4out, dPout                    ! voiding fluxes
     real    dAutC, dAutN, dAutP, dAutChl                ! autolysis fluxes                          
     real    dDetC, dDetN, dDetP, dDetChl                ! voiding fluxes
     real    dCeat, dNeat, dPeat                         ! assimilation fluxes
     real    dPOCout, dPONout, dPOPout                   ! voiding fluxes
     real, dimension(:), allocatable :: vec_dPreyC, vec_dPreyChl, vec_dPreyN, vec_dPreyP, vec_dPreySi


!                                                                                                     
!******************************************************************************* 
!                                                                                                     
	ipnt        = ipoint
		   
	iflux = 0
	
	! segment and species independent items
    maxNrSp   = PMSA(ipnt(   1 ))   !   total nr species implemented in process                (dl)
	nrSp      = PMSA(ipnt(   2 ))   !   nr of species to be modelled                           (dl)                
	nrSpCon   = PMSA(ipnt(   3 ))   !   nr of species dependent items                          (dl)                
	nrSpInd   = PMSA(ipnt(   4 ))   !   nr of species independent items                        (dl)  
    maxNrPr   = PMSA(ipnt(   5 ))   !   nr of prey species implemented                         (dl)
    !nrPr      = PMSA(ipnt(   6 ))   !   nr of prey species user wants to model                 (dl)
    
    ! can I create a module to allocate all of these???? can a module use input to allocate??? 
    ! allocation of prey input array
    allocate( vec_preyC(maxNrPr)        )    
    allocate( vec_preyChl(maxNrPr)      )  
    allocate( vec_preyN(maxNrPr)        )    
    allocate( vec_preyP(maxNrPr)        )    
    allocate( vec_preySi(maxNrPr)       )   
    allocate( vec_CcellPrey(maxNrPr)    )
    allocate( vec_rPrey(maxNrPr)        )    
    allocate( vec_motPrey(maxNrPr)      )  
    allocate( vec_PR(maxNrPr)           )       
    
    ! allocation of food quantity and quality arrays
    allocate( vec_nrPrey(maxNrPr)       )
    allocate( vec_preyFlag(maxNrPr)     )
    allocate( vec_smallerVel(maxNrPr)   )
    allocate( vec_largerVel(maxNrPr)    )
    allocate( vec_encPrey(maxNrPr)      )
    allocate( vec_capturedPrey(maxNrPr) )
    allocate( vec_propPrey(maxNrPr)     )
    allocate( vec_ingNC(maxNrPr)        )
    allocate( vec_ingPC(maxNrPr)        )
    allocate( vec_dPreyC(maxNrPr)       )  
    allocate( vec_dPreyChl(maxNrPr)     )
    allocate( vec_dPreyN(maxNrPr)       )  
    allocate( vec_dPreyP(maxNrPr)       )  
    allocate( vec_dPreySi(maxNrPr)      ) 

   	  
	! length of the PMSA input array. 
	inpItems = nrSpInd   + maxNrSp * nrSpCon + maxNrPr * 10
   
	! segment loop
	do iseg = 1 , noseg
		call dhkmrk(1,iknmrk(iseg),ikmrk1)
		if (ikmrk1.eq.1) then
			
		! species independent items
		PO4          = PMSA(ipnt(  6 ))  !    initial external DIP                                   (gP m-3)
		NH4          = PMSA(ipnt(  7 ))  !    initial external NH4                                   (gN m-3)
		NO3          = PMSA(ipnt(  8 ))  !    initial external NO3                                   (gN m-3)
		Temp         = PMSA(ipnt(  9 ))  !    ambient water temperature                              (oC)               
		PFD          = PMSA(ipnt( 10 ))  !    from rad to photon flux density                        (umol photon m-2)           
		atten        = PMSA(ipnt( 11 ))  !    attenuation of light by water + plankton Chl           (dl)                            
		exat         = PMSA(ipnt( 12 ))  !    -ve exponent of attenuation                            (dl)    
        		  	  
		! species loop
		do iSpec = 0, (nrSp-1)

			spInc = nrSpCon * iSpec
			   
			! species dependent items
			! (number of species independent items + location of input item in vector + species loop)
            protC        = PMSA(ipnt( nrSpInd +  1 + spInc ))   !     C-biomass                                              (gC m-3)  
            protChl      = PMSA(ipnt( nrSpInd +  2 + spInc ))   !     Chl-biomass                                            (gChl m-3)   
            protN        = PMSA(ipnt( nrSpInd +  3 + spInc ))   !     N-biomass                                              (gN m-3)   
            protP        = PMSA(ipnt( nrSpInd +  4 + spInc ))   !     P-biomass                                              (gP m-3)
            AEm          = PMSA(ipnt( nrSpInd +  5 + spInc ))   !     maximum assimilation efficiency (AE)                   (dl)
            AEo          = PMSA(ipnt( nrSpInd +  6 + spInc ))   !     minimum AE                                             (dl)
            alpha        = PMSA(ipnt( nrSpInd +  7 + spInc ))   !     alpha for photosynthesis in protist                    (Figure this out!)  
            CcellProt    = PMSA(ipnt( nrSpInd +  8 + spInc ))   !     C content of protist cell                              (pgC cell-1) 
            ChlCm        = PMSA(ipnt( nrSpInd +  9 + spInc ))   !     maximum cellular Chl:C ratio                           (gChl gC-1)
            ChlCo        = PMSA(ipnt( nrSpInd + 10 + spInc ))   !     minimum cellular Chl:C ratio                           (gChl gC-1)
            CR           = PMSA(ipnt( nrSpInd + 11 + spInc ))   !     catabolic respiration quotient                         (dl)
            FrAut        = PMSA(ipnt( nrSpInd + 12 + spInc ))   !     fraction of mortality to autolysis                     (dl)
            FrDet        = PMSA(ipnt( nrSpInd + 13 + spInc ))   !     fraction of mortality to detritus                      (dl)
            kAE          = PMSA(ipnt( nrSpInd + 14 + spInc ))   !     Control of AE in response to prey quality              (dl)
            KtNH4        = PMSA(ipnt( nrSpInd + 15 + spInc ))   !     Kt for NH4 transport                                   (gN m-3)
            KtNO3        = PMSA(ipnt( nrSpInd + 16 + spInc ))   !     Kt for NO3 transport                                   (gN m-3)
            KtP          = PMSA(ipnt( nrSpInd + 17 + spInc ))   !     Kt for DIP transport                                   (gP m-3)
            MrtRT        = PMSA(ipnt( nrSpInd + 18 + spInc ))   !     mortality at reference temperature                     (dl)
            maxPSreq     = PMSA(ipnt( nrSpInd + 19 + spInc ))   !     maximum C to come from PS                              (dl)
            NCm          = PMSA(ipnt( nrSpInd + 20 + spInc ))   !     N:C that totally represses NH4 transport               (gN gC-1)
            NCo          = PMSA(ipnt( nrSpInd + 21 + spInc ))   !     minimum N-quota                                        (gN gC-1)
            NCopt        = PMSA(ipnt( nrSpInd + 22 + spInc ))   !     N:C for growth under optimal conditions                (gN gC-1)
            NO3Cm        = PMSA(ipnt( nrSpInd + 23 + spInc ))   !     N:C that totally represses NO3 transport               (gN gC-1)
            NO3Copt      = PMSA(ipnt( nrSpInd + 24 + spInc ))   !     N:C for growth on NO3 under optimal conditions         (gN gC-1)
            optCR        = PMSA(ipnt( nrSpInd + 25 + spInc ))   !     proportion of prey captured by starved Prot               (dl)
            PCm          = PMSA(ipnt( nrSpInd + 26 + spInc ))   !     PC maximum quota                                       (gP gC-1)
            PCo          = PMSA(ipnt( nrSpInd + 27 + spInc ))   !     PC minimum quota                                       (gP gC-1)
            PCoNCm       = PMSA(ipnt( nrSpInd + 28 + spInc ))   !     maximum NC when PC is minimum (PCu = 0)                (gN gC-1)
            PCoNCopt     = PMSA(ipnt( nrSpInd + 29 + spInc ))   !     optimum NC when PC is minimum (PCu = 0)                (gN gC-1)
            PCopt        = PMSA(ipnt( nrSpInd + 30 + spInc ))   !     PC optimum quota                                       (gP gC-1)
            PSDOC        = PMSA(ipnt( nrSpInd + 31 + spInc ))   !     proportion of current PS being leaked as DOC           (dl)
            Q10          = PMSA(ipnt( nrSpInd + 32 + spInc ))   !     Q10 for UmRT                                           (dl)
            rProt        = PMSA(ipnt( nrSpInd + 33 + spInc ))   !     radius of nutrient repleted protist cell               (um)
            redco        = PMSA(ipnt( nrSpInd + 34 + spInc ))   !     C respired to support nitrate reduction for NH4        (gC gN-1)
            relPhag      = PMSA(ipnt( nrSpInd + 35 + spInc ))   !     rel. phagotrophy in dark : in light                    (dl)
            relPS        = PMSA(ipnt( nrSpInd + 36 + spInc ))   !     relative PSmax:Umax on phototrophy                     (dl)
            ReUmNH4      = PMSA(ipnt( nrSpInd + 37 + spInc ))   !     max. growth rate supported by NH4-N:Umax               (dl)
            ReUmNO3      = PMSA(ipnt( nrSpInd + 38 + spInc ))   !     max. growth rate supported by NO3-N:Umax               (dl)
            RT           = PMSA(ipnt( nrSpInd + 39 + spInc ))   !     reference temperature for UmRT                         (deg C)
            SDA          = PMSA(ipnt( nrSpInd + 40 + spInc ))   !     specific dynamic action                                (dl)
            UmRT         = PMSA(ipnt( nrSpInd + 41 + spInc ))   !     maximum growth rate at reference T                     (d-1) 
            
            if (protC <= 1.0E-9) then 
                cycle
            end if
                        
            ! Calculate the nutrient quota of the cell-------------------------------------------------------------------------------    						 
			! Units: gNut gC-1  
            NC   = quota(protN, protC)
			PC   = quota(protP, protC)
            ChlC = quota(protChl, protC)
                        
            ! Calculate maximum growth and respiration -------------------------------------------------------------------------------    
            ! Units: gC gC-1 d-1
            UmT = Q10rate(UmRT, Q10, Temp, RT)
            BR  = basal_respiration(UmT, CR)  

            !! PHOTOTROPHY -------------------------------------------------------------------------------    
            
            ! Calculate nutrient status within cell compared to ideal status (nutrient status = 1) --------------------------------------- 
            ! Determine minimum of N-P-Si limitation; Liebig-style limitation of growth (NPCu)
            ! Units: dl
            NCu = statusNC(NC, NCo, NCopt)                        
            PCu = statusPC(PC, PCo, PCopt)
            NPCu = min(NCu, PCu)	  
                        
            ! swimming speed -------------------------------------------------------------------------------    
            ! Units: m s-1
            motProt = motility(rProt) 
            
            ! Calculate uptake for the nutrients --------------------------------------- 
            ! Units: gNut gC-1 d-1    
            upP = uptakeP(PC, PCo, PCopt, PCm, UmT, PO4, KtP)
            upNH4 = uptakeNH4(PCoNCopt, PCoNCm, PCu, NCu, NC, NCo, NCopt, NCm, UmT, ReUmNH4, NH4, KtNH4)             
            upNO3 = uptakeNO3(PCoNCm, PCu, NC, NC, NCo, NO3Copt, NO3Cm, UmT, ReUmNO3, NO3, KtNO3) 
                   
            ! Calculate photosynthesis related equation --------------------------------------- 
            ! Units: gC gC-1 d-1
            ! I do not like the variable maxPSreq. No measureable and pretty "strong" influence
            PSqm = plateauPS(UmT, maxPSreq, relPS, NCopt, redco, NPCu, BR, PSDOC)
            PS   = grossPS(ChlC, PFD, exat, atten, PSqm, alpha)
            Cfix = netPS(PS, PSDOC)
            
            ! rate of (positive) net phototrophy
            ! Units: gC gC-1 d-1     
            CfixPS = Cfix - totalRespiration(redco, upNO3, upNH4, 0.0, 0.0, 0.0, 0.0)
                        
            ! Calculate chlorophyll synthesis and degradation --------------------------------------- 
            ! Units: gChl gC-1 d-1          
            synChl = synthesisChl(ChlC, ChlCo, ChlCm, UmT, maxPSreq, NPCu, Cfix, PSqm)
            degChl = degradeChl(ChlC, ChlCm, UmT, NPCu)

            !! PHAGOTRPOHY -------------------------------------------------------------------------------    
            
            ! ideally I would like to have this in four subroutines:
            ! foodQuantity(, sumCP, ingNC, ingPC) 
            ! foodQuality(, ppNC, ppPC, stoichP, opAE) 
            ! ingestion(, maxIng, ingSat, ingC, ingN, ingP)           
            ! assimilation(, assC, assN, assP) 
            ! but it keeps crashing when I change the equations into subroutines 
                        
            do iPrey = 0, (maxNrPr - 1)
                !prey specific input
                ! independentItems + all input items of all zoo species + first prey item + current PreyNumber * total nr of prey specific items
                vec_preyC(iPrey + 1)         = PMSA(ipnt( nrSpInd + maxNrSp * nrSpCon + 1 + iPrey * 10))   !      C-biomass                                              (gC m-3)   
                vec_preyChl(iPrey + 1)       = PMSA(ipnt( nrSpInd + maxNrSp * nrSpCon + 2 + iPrey * 10))   !      Chl-biomass                                            (gC m-3)  
                vec_preyN(iPrey + 1)         = PMSA(ipnt( nrSpInd + maxNrSp * nrSpCon + 3 + iPrey * 10))   !      N-biomass                                              (gN m-3)  
                vec_preyP(iPrey + 1)         = PMSA(ipnt( nrSpInd + maxNrSp * nrSpCon + 4 + iPrey * 10))   !      P-biomass                                              (gP m-3)  
                vec_preySi(iPrey + 1)        = PMSA(ipnt( nrSpInd + maxNrSp * nrSpCon + 5 + iPrey * 10))   !      Si-biomass                                             (gP m-3)  
                vec_CcellPrey(iPrey + 1)     = PMSA(ipnt( nrSpInd + maxNrSp * nrSpCon + 6 + iPrey * 10))   !      C content of protist cell                              (pgC cell-1) 
                vec_rPrey(iPrey + 1)         = PMSA(ipnt( nrSpInd + maxNrSp * nrSpCon + 7 + iPrey * 10))   !      radius of nutrient repleted protist cell               (um)
                vec_motPrey(iPrey + 1)       = PMSA(ipnt( nrSpInd + maxNrSp * nrSpCon + 8 + iPrey * 10))   !      swimming velocity                                      (m s-1)
                vec_PR(iPrey + 1)            = PMSA(ipnt( nrSpInd + maxNrSp * nrSpCon + 9 + iSpec + iPrey * 10))   !      handling index of prey 1 by pred 1                     (dl)  
                
                ! if loop to protect against small preys [dl]
                if (vec_preyC(iPrey + 1) >= 1.0E-5) then 
                    vec_preyFlag(iPrey + 1) = 1.0
                else 
                    vec_preyFlag(iPrey + 1) = 0.0
                end if 
                                
            end do
            
            ! for output [dl]
            preyFlag = sum(vec_preyFlag)
            
            !! FOOD QUANTITY -------------------------------------------------------------------------------                
            ! reduction of phagotrophy during night 
            ! Units: dl
            relPhag = 1.0 - relPhag
            lightInh = lightInhibition(PFD, relPhag)
            
            ! cell abundance of prey per m3
            ! Units: nr cells m-3 (1e12: transform between g (preyC) and pg (CcontentPrey))
            vec_nrPrey = lightInh * vec_preyFlag * 1e12 * vec_preyC / vec_CcellPrey
            
            ! encounter Rate 
            ! Units: prey predator-1 d-1
            vec_smallerVel = min(vec_motPrey, motProt  )
            vec_largerVel = max(vec_motPrey, motProt  )
            vec_encPrey = (24.0 * 60.0 * 60.0) * PI_8 *(vec_rPrey / 1E6 + rProt / 1E6)**2 * vec_nrPrey * ((vec_smallerVel**2 + 3 * vec_largerVel**2 + 4 * wTurb**2) * ((vec_largerVel**2 + wTurb**2)**(-0.5))) * 3.0**(-1.0)
            
            ! potential C-specific capture of prey
            ! Units: gC gC-1 d-1 
            vec_capturedPrey = vec_encPrey * vec_PR * optCR * vec_CcellPrey / CcellProt  
            
            ! sum of all potential C-specific prey captures 
            ! Units: gC gC-1 d-1 (1.0E-20: protection against division by 0 in following line)
            sumCP = sum(vec_capturedPrey) + 1.0E-20
             
            ! proportion iPrey of total prey 
            ! Units: dl
            vec_propPrey = vec_capturedPrey / sumCP
            
            ! total captured prey Nut:C
            ! Units: gNut gC-1 d-1
            vec_ingNC = vec_propPrey * vec_preyN / vec_preyC             
            vec_ingPC = vec_propPrey * vec_preyP / vec_preyC 
            ingNC = sum(vec_ingNC)
            ingPC = sum(vec_ingPC)
     
            ! FOOD QUALITY -------------------------------------------------------------------------------   
            ! quota of captured prey in relation to predator    
            ! Units: dl 
            ppNC = ingNC / NCopt
            ppPC = ingPC / PCopt
            ! determine limiting nutrient in prey or set to 1 if preNut > predNut
            ! Units: dl
            stoichP = min(ppNC, ppPC, 1.0)
                        
            !! assimilation efficiency for prey 
            !! Units: dl  
            opAE = (AEo + (AEm - AEo) * stoichP / (stoichP + kAE) * (1.0 + kAE)) * stoichP + 1.0E-20
             
            ! INGESTION  ------------------------------------------------------------------------------- 
            ! required predation 
            ! can also fall below 0.0 if CfixPS is larger than Umt+BR
            ! of course the question remains what UmT is worth if Cfix > UmT 
            ! Units: gC gC-1 d-1
            ! I removed the CfixPS part for now because it produced jumps which didn't look nice. 
            ! I need to rethink this part
            !reqPred = ((UmT + BR - CfixPS) / (1.0 - SDA)) / opAE
            reqPred = ((UmT + BR - 0.0) / (1.0 - SDA)) / opAE

            ! maximum ingestion if needed
            ! if 0.0 then there is no need for ingestion because of high Cfix
            ! Units: gC gC-1 d-1
            ! with the upper adjustment (CfixPS == 0) this line becomes obsolete e.g. reqPred > 0.0 => reqPred = maxIng
            maxIng = max(0.0, reqPred)            

            ! half saturation constant for satiation feedback (see paper by Flynn and Mitra 2016)
            ! Units: gC gC-1 d-1
            KI = (maxIng / 4.0)
            
            ! ingestion with satiation feedback included
            ! Units: gC gC-1 d-1
            ingSat = maxIng * sumCP / (sumCP + KI)
            
            ! ingestion of C
            ! Units: gC gC d-1
            ingC = min(ingSat, sumCP)	
            ingN = ingC * ingNC            
            ingP = ingC * ingPC    
            
            ! ASSIMILATION ------------------------------------------------------------------------------- 
            ! assimilation of ingested prey
            ! Units: gC gC-1 d-1 / gNut gC-1 d-1
            assC = ingC * opAE
            assN = assC * NCopt            
            assP = assC * PCopt
                        
            ! Calculate respiration and C-growth  --------------------------------------- 
            ! Units: gC gC-1 d-1             
            ! 0.0 because it cannot assimilate prey
            if (protC > 1.0E-5) then 
                totR = totalRespiration(redco, upNO3, upNH4, assC, assN, SDA, BR)
            else
                totR = 0.0
            end if
            !totR = totalRespiration(redco, upNO3, upNH4, assC, assN, SDA, BR)
			Cu	 = Cfix + assC - totR
                        
            ! Calculate mortality  --------------------------------------- 
            ! Units: gC gC-1 d-1 
            if (protC > 1.0E-5) then 
                mrt = Q10rate(MrtRT, Q10, Temp, RT) 
            else 
                mrt = 0.0
            end if
            !mrt = Q10rate(MrtRT, Q10, Temp, RT) 
            mrtFrAut = mortality(mrt, FrAut)           
            mrtFrDet = mortality(mrt, FrDet)      
				   
			! Output -------------------------------------------------------------------
			   
			! (input items + position of specific output item in vector + species loop * total number of output) 
            PMSA(ipnt( inpItems +  1 + iSpec * 41 )) = NC 
            PMSA(ipnt( inpItems +  2 + iSpec * 41 )) = PC 
            PMSA(ipnt( inpItems +  3 + iSpec * 41 )) = ChlC 
            PMSA(ipnt( inpItems +  4 + iSpec * 41 )) = UmT 
            PMSA(ipnt( inpItems +  5 + iSpec * 41 )) = BR
            PMSA(ipnt( inpItems +  6 + iSpec * 41 )) = NCu 
            PMSA(ipnt( inpItems +  7 + iSpec * 41 )) = PCu 
            PMSA(ipnt( inpItems +  8 + iSpec * 41 )) = NPCu
            PMSA(ipnt( inpItems +  9 + iSpec * 41 )) = motProt
            PMSA(ipnt( inpItems + 10 + iSpec * 41 )) = upP 
            PMSA(ipnt( inpItems + 11 + iSpec * 41 )) = upNH4 
            PMSA(ipnt( inpItems + 12 + iSpec * 41 )) = upNO3 
            PMSA(ipnt( inpItems + 13 + iSpec * 41 )) = PSqm 
            PMSA(ipnt( inpItems + 14 + iSpec * 41 )) = PS
            PMSA(ipnt( inpItems + 15 + iSpec * 41 )) = Cfix 
            PMSA(ipnt( inpItems + 16 + iSpec * 41 )) = CfixPS
            PMSA(ipnt( inpItems + 17 + iSpec * 41 )) = synChl
            PMSA(ipnt( inpItems + 18 + iSpec * 41 )) = degChl
            PMSA(ipnt( inpItems + 19 + iSpec * 41 )) = sumCP
            PMSA(ipnt( inpItems + 20 + iSpec * 41 )) = ingNC
            PMSA(ipnt( inpItems + 21 + iSpec * 41 )) = ingPC
            PMSA(ipnt( inpItems + 22 + iSpec * 41 )) = ppNC
            PMSA(ipnt( inpItems + 23 + iSpec * 41 )) = ppPC
            PMSA(ipnt( inpItems + 24 + iSpec * 41 )) = stoichP
            PMSA(ipnt( inpItems + 25 + iSpec * 41 )) = opAE
            PMSA(ipnt( inpItems + 26 + iSpec * 41 )) = reqPred
            PMSA(ipnt( inpItems + 27 + iSpec * 41 )) = maxIng
            PMSA(ipnt( inpItems + 28 + iSpec * 41 )) = ingSat
            PMSA(ipnt( inpItems + 29 + iSpec * 41 )) = ingC  
            PMSA(ipnt( inpItems + 30 + iSpec * 41 )) = assC  
            PMSA(ipnt( inpItems + 31 + iSpec * 41 )) = ingN
            PMSA(ipnt( inpItems + 32 + iSpec * 41 )) = ingP
            PMSA(ipnt( inpItems + 33 + iSpec * 41 )) = assN
            PMSA(ipnt( inpItems + 34 + iSpec * 41 )) = assP
            PMSA(ipnt( inpItems + 35 + iSpec * 41 )) = totR 
            PMSA(ipnt( inpItems + 36 + iSpec * 41 )) = Cu
            PMSA(ipnt( inpItems + 37 + iSpec * 41 )) = mrt 
            PMSA(ipnt( inpItems + 38 + iSpec * 41 )) = mrtFrAut 
            PMSA(ipnt( inpItems + 39 + iSpec * 41 )) = mrtFrDet
            PMSA(ipnt( inpItems + 40 + iSpec * 41 )) = preyFlag
            PMSA(ipnt( inpItems + 41 + iSpec * 41 )) = lightInh

    		! FLUXES -------------------------------------------------------------------   
            ! Protist gains------------------------------------------------------------   
            ! Protist growth through assimilation -----------------------------------------------------
			! gX m-3 d-1 assimilation of X from prey
			dCeat = protC * assC
			dNeat = protC * assN	
			dPeat = protC * assP
            
            ! gNut m-3 d-1   uptake of nutrients into algal biomass
			dNH4up = protC * upNH4	
			dNO3up = protC * upNO3	
			dPup   = protC * upP
            			
			! gC m-3 d-1   total contribution to biomass growth from C-fixation
			dCfix = protC * Cfix
            			
            ! gChl m-3 d-1 Chl synthesis or degradation
			dChlsyn	= protC * synChl	
            dChldeg = protC * degChl
            
            ! Protist losses-----------------------------------------------------------
            ! gC m-3 d-1   total respiration rate
			dCresp = protC * totR	
            			
            ! gC m-3 d-1   release of DOC 
			dDOCleak = protC * (PS - Cfix)
            	           
            ! gC m-3 d-1   voiding of C as DOC if NC falls below NCo
			if (NC < NCo) then 
			    dDOCvoid = protC - protN / NCo
			else 
			    dDOCvoid = 0.0
            end if 			
            
            ! gNut m-3 d-1 voiding of nutrient P and N if interanl maximum is reached
            dNH4out = voiding(protN, protC, NCm)
            dPout   = voiding(protP, protC, PCm)
            
            ! gX m-3 d-1  rate of voiding of X as particulates
			dPOCout = protC * (ingC - assC)
			dPONout = protC * (ingN - assN)	
			dPOPout = protC * (ingP - assP) 
                                    
            ! gNut m-3 d-1 mortality
            dAutC       = protC * mrtFrAut
			dDetC       = protC * mrtFrDet	
            dAutN       = protN * mrtFrAut
			dDetN       = protN * mrtFrDet			
            dAutP       = protP * mrtFrAut
			dDetP       = protP * mrtFrDet			
            dAutChl     = protChl * mrtFrAut
			dDetChl     = protChl * mrtFrDet
			  
			! (1 + SpeciesLoop * (nr of fluxes per individual species) + total number of fluxes) 
            fl (  1 + (25 + maxNrPr * 5) * iSpec + iflux )  = dNH4up    
            fl (  2 + (25 + maxNrPr * 5) * iSpec + iflux )  = dNO3up    
            fl (  3 + (25 + maxNrPr * 5) * iSpec + iflux )  = dPup      
            fl (  4 + (25 + maxNrPr * 5) * iSpec + iflux )  = dCfix     
            fl (  5 + (25 + maxNrPr * 5) * iSpec + iflux )  = dChlsyn   
            fl (  6 + (25 + maxNrPr * 5) * iSpec + iflux )  = dChldeg   
            fl (  7 + (25 + maxNrPr * 5) * iSpec + iflux )  = dCresp    
            fl (  8 + (25 + maxNrPr * 5) * iSpec + iflux )  = dDOCleak    
            fl (  9 + (25 + maxNrPr * 5) * iSpec + iflux )  = dDOCvoid    
            fl ( 10 + (25 + maxNrPr * 5) * iSpec + iflux )  = dNH4out   
            fl ( 11 + (25 + maxNrPr * 5) * iSpec + iflux )  = dPout     
            fl ( 12 + (25 + maxNrPr * 5) * iSpec + iflux )  = dCeat     
            fl ( 13 + (25 + maxNrPr * 5) * iSpec + iflux )  = dNeat     
            fl ( 14 + (25 + maxNrPr * 5) * iSpec + iflux )  = dPeat     
            fl ( 15 + (25 + maxNrPr * 5) * iSpec + iflux )  = dPOCout   
            fl ( 16 + (25 + maxNrPr * 5) * iSpec + iflux )  = dPONout   
            fl ( 17 + (25 + maxNrPr * 5) * iSpec + iflux )  = dPOPout   
            fl ( 18 + (25 + maxNrPr * 5) * iSpec + iflux )  = dAutC     
            fl ( 19 + (25 + maxNrPr * 5) * iSpec + iflux )  = dDetC     
            fl ( 20 + (25 + maxNrPr * 5) * iSpec + iflux )  = dAutN     
            fl ( 21 + (25 + maxNrPr * 5) * iSpec + iflux )  = dDetN     
            fl ( 22 + (25 + maxNrPr * 5) * iSpec + iflux )  = dAutP     
            fl ( 23 + (25 + maxNrPr * 5) * iSpec + iflux )  = dDetP     
            fl ( 24 + (25 + maxNrPr * 5) * iSpec + iflux )  = dAutChl   
            fl ( 25 + (25 + maxNrPr * 5) * iSpec + iflux )  = dDetChl   
            !iSPec * 25
            
            ! Prey losses through pred ing. ----------------------------------------------------  
                            
            ! ingestion of nut of iPrey through iPred gNut m-3 d-1	
            vec_dPreyC    = protC * (ingC * vec_propPrey)	
            vec_dPreyChl  = vec_dPreyC * (vec_preyChl / vec_preyC)            
            vec_dPreyN    = vec_dPreyC * (vec_preyN / vec_preyC) 
            vec_dPreyP    = vec_dPreyC * (vec_preyP / vec_preyC) 
            vec_dPreySi   = vec_dPreyC * (vec_preySi / vec_preyC)   
                        
            ! loop over prey ingestion fluxes
            do iPrey = 0, (maxNrPr - 1)                                
                ! (nr prey independent fluxes + prey Flux # + loop) + (move on to next predator) + total number of fluxes
                fl ( (25 + 1 + iPrey * 5) + (25 + maxNrPr * 5) * iSpec + iflux ) = vec_dPreyC(iPrey + 1)  
                fl ( (25 + 2 + iPrey * 5) + (25 + maxNrPr * 5) * iSpec + iflux ) = vec_dPreyChl(iPrey + 1)
                fl ( (25 + 3 + iPrey * 5) + (25 + maxNrPr * 5) * iSpec + iflux ) = vec_dPreyN(iPrey + 1)  
                fl ( (25 + 4 + iPrey * 5) + (25 + maxNrPr * 5) * iSpec + iflux ) = vec_dPreyP(iPrey + 1)  
                fl ( (25 + 5 + iPrey * 5) + (25 + maxNrPr * 5) * iSpec + iflux ) = vec_dPreySi(iPrey + 1) 
            end do  
            
            !write(82,*) protC, Cfix, ingC, assC, totR, mrt, iseg
            !write(78,*) protC
                        
            if ( isnan(protC) ) write ( 70, '(''ERROR: NaN in protC in segment:'', i10)' )    iseg
            if ( isnan(Cfix) )  write ( 71, '(''ERROR: NaN in Cfix in segment:'', i10)' )    iseg
            if ( isnan(totR) )  write ( 72, '(''ERROR: NaN in totR in segment:'', i10)' )    iseg
            if ( isnan(mrt) )   write ( 73, '(''ERROR: NaN in mrt in segment:'', i10)' )    iseg
            if ( isnan(NC) )    write ( 74, '(''ERROR: NaN in NC in segment:'', i10)' )    iseg
            if ( isnan(PC) )    write ( 75, '(''ERROR: NaN in PC in segment:'', i10)' )    iseg
            if ( isnan(ChlC) )  write ( 76, '(''ERROR: NaN in ChlC in segment:'', i10)' )    iseg
            if ( isnan(ingC) )    write ( 77, '(''ERROR: NaN in ingC in segment:'', i10)' )    iseg
  			   
		enddo ! end loop over species 

		endif ! end if check for dry cell 

		!allocate pointers
		iflux = iflux + noflux
		ipnt = ipnt + increm

	enddo ! end loop over segments
	return
end ! end subroutine 

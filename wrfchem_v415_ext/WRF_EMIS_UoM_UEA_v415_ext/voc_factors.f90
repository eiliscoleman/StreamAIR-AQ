MODULE voc_factors
IMPLICIT NONE
INTEGER, PARAMETER :: nvocs = 47, nvocsrc = 9
REAL :: factor(nvocs,nvocsrc)

! List of the 47 of the top 50 NMVOC species from the NAEI dataset (propylene, dichloromethane and tetrachloroethene
! are excluded as they are not mappable onto the current WRF chemical scheme) and the fraction of the whole NMVOC
! emission which comes from each species/source. The numbers are calculated from Table 2.13 of 
!
! "UK Emissions of Air Pollutants 1970 to 2005", UK Emissions Inventory Team, AEA Energy & Environment: C J Dore, 
! J D Watterson, T P Murrells, N R Passant, M M Hobson, S L Choudrie, G Thistlethwaite, A Wagner, J Jackson, Y Li, 
! T Bush, K R King, J Norris, P J Coleman, C Walker, R A Stewart, J W L Goodwin, I Tsagatakis, C Conolly, M K Downes,
! N Brophy, M R Hann, 2005.
!
! which is available from the NAEI website. These numbers are for 2005 and may need to be revised if speciation for
! other years differs significantly.
!
! Numbers by Øivind Hodnebrog.
! Conversion to module and documentation by T. Pugh 10/09/09.
!
! Key to species:
! 1) Ethanol
! 2) Butane
! 3) Ethane
! 4) Propane
! 5) Toluene
! 6) Methanol
! 7) Ethylene
! 8) 2-Methylbutane
! 9) Pentane
! 10) Acetone
! 11) Hexane
! 12) m-Xylene
! 13) 2-methylpropane
! 14) Formaldehyde
! 15) Trichloroethene
! 16) Benzene
! 17) 2-Butanone
! 18) Butyl acetate
! 19) Decane
! 20) Ethylbenzene
! 21) 1,2,4-Trimethylbenzene
! 22) 2-propanol
! 23) Heptane
! 24) Ethyl acetate
! 25) o-Xylene
! 26) p-Xylene
! 27) Octane
! 28) 4-Methyl-2-pentanone
! 29) Acetylene
! 30) Nonane
! 31) 2-Methylpropene
! 32) Undecane
! 33) Methyl acetate
! 34) 1-butanol
! 35) Acetaldehyde
! 36) 2-Methylpentane
! 37) 2_Butoxyethanol
! 38) 1,3,5-Trimethylbenzene
! 39) 1-Propanol
! 40) 1,3-Butadiene
! 41) Dipentene
! 42) 2-Butene
! 43) 1-Methoxy-2-propanol
! 44) 1,2,3-Trimethylbenzene
! 45) Methylethylbenzene
! 46) 4-Methyldecane
! 47) 2-Pentene

DATA factor(1,:) / 0.,0.161834178477008,0.0181940700808625,0.462673027818842,0.,0.150546370097124,0.,0.,0.0291341616103246 /
DATA factor(2,:) / 0.0502257336343115,0.0512287088805097,0.117250673854447, &
0.0387800698830802,0.265133733432655,0.0706511794714153,0.085517824036047, &
0.0103301244010726,0.00245593759029182 /
DATA factor(3,:) / 0.0359292701279157,0.103107528279808,0.0630053908355795, &
0.0120447520494557,0.226545699869288,0.,0.0211213067048415,0.0113411578530924,0.290282192044688 /
DATA factor(4,:) / 0.031790820165538,0.0484549040003467,0.0589622641509434, &
0.0194110334632442,0.132447399120156,0.0137964853727979,0.00752153199877966,0.00817618356850851,0.274342675527304 /
DATA factor(5,:) / 0.0220090293453725,0.0317253933168639,0.045822102425876, &
0.0330096761188012,0.00119467596171415,0.0412211170395743,0.0962075520405529,0.0800914325904435,0.0168544736588655 /
DATA factor(6,:) / 0.,0.,0.,0.0174539712404247,0.,0.102006162674103,0.,0.,0.00756043532697679 /
DATA factor(7,:) / 0.00884123401053424,0.182334330169462,0.0350404312668464,&
0.0489349549791695,0.000189742652742835,0.,0.0924643871300838,0.171304233153106,0.0514302224790523/
DATA factor(8,:) / 0.0163656884875846,0.0474147271702856,0.0431266846361186, &
0.00839940868162881,0.0544912788654795,0.000179317714394455,0.111919457417099,0.0143303002329773,0.00158913608783589/
DATA factor(9,:) /0.0393152746425884,0.0334156806657132,0.117924528301887, &
0.0160344711732294,0.0918916639727895,0.00159556170359148,0.0540236089272723,0.00611015868829399,0.00211884811711451/
DATA factor(10,:) /0.00432656132430399,0.00034672561002037,0.00673854447439353, &
0.0164040451552211,0.,0.065637602558754,0.00713430804252423,0.00246164666578751,0.00014446691707599/
DATA factor(11,:) /0.0263355906696764,0.00429072942400208,0.0320080862533693, &
0.0364198360435425,0.0470561778802232,0.00900614071683171,0.0496937410527798,0.00382434392720559,0.0111239526148512/
DATA factor(12,:) /0.0940556809631302,0.00303384908767824,0.0252695417789757, &
0.0176975540921919,0.000428677845085665,0.0434058655190333,0.0289009880077914,0.0235614752296804,0.00775305788307811/
DATA factor(13,:) /0.00319789315274643,0.0127421661682486,0.00943396226415094, &
0.00187306813600323,0.0663466809090781,0.00345827020617878,0.0373964469268499, &
0.00549474702184711,0.000770490224405278/
DATA factor(14,:) /0.579194883370956,0.0319854375243791,0.201145552560647, &
0.00405691439322672,0.00107520836554273,0.0000878290846013657,0.0494942620450119,0.118466745791024,0.182895117018203/
!DATA factor(15,:) /0.,0.,0.,0.00762666308291896,0.,0.047698512028925,0.,0.,0.00611576615621689/
DATA factor(15,:) /0.0174943566591422,0.213409612967538,0.123989218328841, &
0.0137750302378713,0.00422353089994237,0.,0.034110910328319,0.107169545914106,0.0496484638351151/
DATA factor(16,:) /0.0112866817155756,0.0345858795995319,0.011455525606469, &
0.0528238811987636,0.0000983850791999888,0.,0.0410340053037948,0.0497604290298475,0.00293749398054512/
!DATA factor(18,:) /0.,0.,0.,0.0286755812390808,0.,0.0360282224125186,0.,0.,0.00683810074159684/
DATA factor(17,:) /0.,0.,0.,0.00561920440800968,0.,0.0437205864055215, &
0.00200652413695994,0.000351663809398215,0.00139651353173457/
DATA factor(18,:) /0.,0.,0.,0.00167148232764413,0.,0.0393840253533291,0.,0.,0.00216700375613984/
DATA factor(19,:) /0.,0.000433407012525463,0.,0.00683711866684585, &
0.000154605124457125,0.029667932869303,0.00738072328741405,0.0554310079563937,0./
DATA factor(20,:) /0.0276523702031603,0.000996836128808564,0.0181940700808625, &
0.0139346190028222,0.000126495101828557,0.0166472710771505,0.0241252258806412,0.0155171655896962,0.0124723105075604/
DATA factor(21,:) /0.,0.,0.,0.00436769251444698,0.0000281100226285682, &
0.0193772917901763,0.0311069911525193,0.0203085849927469,0./
DATA factor(22,:) /0.,0.000130022103757639,0.,0.00622396183308695,0.,0.0298911651259981,0.,0.,0.00163729172686122/
DATA factor(23,:) /0.00357411587659895,0.0146058163221081,0.000336927223719677, &
0.00251982260448864,0.0458825844354805,0.00541978642894261,0.0107483983009082,0.00470350345070113,0./
DATA factor(24,:) /0.,0.,0.,0.0111712135465663,0.,0.026345065835218,0.,0.,0.00231147067321583/
DATA factor(25,:) /0.0208803611738149,0.00173362805010185,0.00876010781671159, &
0.00608117188549926,0.00020379766405712,0.0109164233069114,0.0261082818990402,0.0207921227306695,0.00438216315130502/
DATA factor(26,:) /0.000752445447705041,0.0023403978676375,0.0148247978436658, &
0.0074334766832415,0.0000913575735428467,0.011728842339474,0.0223416488700101,0.0181986021363576,0.00601945487816623/
DATA factor(27,:) /0.,0.000780132622545833,0.,0.00155389060610133,0.0401832773475383, &
0.00466957966463928,0.00459975123794325,0.00167040309464152,0./
!DATA factor(30,:) 0.,0.,0.,0.00173027818841554,0.,0.0212436598379553,0.,0.,0.0126167774246364/
DATA factor(28,:) /0.,0.,0.,0.00551841150383013,0.,0.0209801725841512,0.,0.,0./
DATA factor(29,:) /0.00300978179082017,0.000303384908767824,0.045822102425876, &
0.00500604757425077,0.0000843300678857047,0.,0.0358358170425477,0.042946942722757,0./
DATA factor(30,:) /0.,0.000650110518788194,0.,0.00419130493213278,0.000407595328114239, &
0.0179061546231034,0.00176010889207012,0.0129676029715592,0./
DATA factor(31,:) /0.,0.00264378277640532,0.,0.00531682569547104,0.001286033535257, &
0.,0.0291474032526812,0.0569255791463361,0.000481556390253299/
DATA factor(32,:) /0.,0.,0.,0.00356134928101062,0.,0.0154579188898404,0.,0.0244845927293507,0./
DATA factor(33,:) /0.,0.,0.,0.0415098777046096,0.,0.,0.,0.,0./
DATA factor(34,:) /0.,0.,0.,0.00199905926622766,0.,0.0159226811291893,0.,0.,0.000674178946354618/
DATA factor(35,:) /0.,0.,0.,0.00615676656363392,0.,0.,0.0243599070662505,0.0601784693832696,0./
DATA factor(36,:) /0.000940556809631302,0.000216703506262731,0.00269541778975741, &
0.00769385835237199,0.00921305991651323,0.00430362514546692,0., &
0.000307705833223438,0.00524896465376096/
DATA factor(37,:) /0.,0.,0.,0.00084834027684451,0.,0.0133280635882573,0.,0.,0./
DATA factor(38,:) /0.,0.,0.,0.00158748824082785,0.,0.00670062724604586,0.0118044636361503,0.0112092839245681,0./
DATA factor(39,:) /0.,0.,0.,0.000529162746942615,0.,0.0127169195412394,0.,0.,0.00394876240007705/
DATA factor(40,:) /0.00018811136192626,0.,0.,0.00320017470770058, &
0.0000351375282857103,0.,0.0186688883152238,0.027385819156886,0.000722334585379948/
DATA factor(41,:) /0.,0.,0.,0.000109192312861175,0.,0.0120911373134547,0.,0.,0./
DATA factor(42,:) /0.00018811136192626,0.0104884497031162,0.,0.00115911839806478, &
0.00404784325851382,0.,0.0164863532890568,0.00413204976042903,0.00192622556101319/
DATA factor(43,:) /0.,0.,0.,0.000806343233436366,0.,0.0110518264790052,0.,0.,0./
DATA factor(44,:) /0.,0.,0.,0.00147829592796667,0.,0.00674820133353827,0.00704043556828049,0.0091432590443536,0./
DATA factor(45,:) /0.,0.,0.,0.00193186399677463,0.,0.00945992432060543,0.,0.,0./
DATA factor(46,:) /0.,0.,0.,0.00206625453568069,0.,0.00900248117163999,0.,0.,0./
DATA factor(47,:) /0.00244544770504138,0.00476747713778009,0., &
0.0000923934954979169,0.00713291824199918,0.,0.00993874821055596,0.000923117499670315,0./
END MODULE
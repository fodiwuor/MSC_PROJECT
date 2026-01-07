**Preparing a working directory
cd "C:/Users/THINKPAD/Desktop/DesktopDamp/MSC_STATISTICSNOTES/Personal_lerning notes/MSC_Project/documents/referencesTo_ReadForMethodFavouring"
**Reading the data
set more off
*use "C:\Users\THINKPAD\AppData\Local\Temp\eb1cd984-104b-4f5d-a114-8e20e6e8869c_all_admissions_temp.dta.zip.69c\all_admissions_temp.dta",clear
use "C:\Users\THINKPAD\Downloads\all_admissions_temp.dta\all_admissions_temp.dta", clear
// Mark PCV10 serotypes
gen pcv10=.
replace pcv10=1 if spn_serotype=="4" | spn_serotype=="6B" | spn_serotype=="9V" /*
*/  | spn_serotype=="14"| spn_serotype=="18C" | spn_serotype=="19F" | spn_serotype=="23F" /*
*/ | spn_serotype=="1" | spn_serotype=="5" | spn_serotype=="7F"
replace pcv10=0 if spn_serotype!="" & pcv10!=1
capture label define pcv10 0 "NVT" 1 "VT"
capture label values pcv10 pcv10
label var pcv10 "Indicator for PCV-10 vaccine serotype"
// Define VT/NVT where serotype 6A is excluded from NVT (this code fragment was added by JOjal on 13th Oct 18)
gen pcv10_6A=.
replace pcv10_6A=1 if spn_serotype=="4" | spn_serotype=="6B" | spn_serotype=="9V" /*
*/  | spn_serotype=="14"| spn_serotype=="18C" | spn_serotype=="19F" | spn_serotype=="23F" /*
*/ | spn_serotype=="1" | spn_serotype=="5" | spn_serotype=="7F" 
replace pcv10_6A=0 if spn_serotype!="" & spn_serotype!="6A" & pcv10_6A!=1
capture label values pcv10_6A pcv10
label var pcv10_6A "Indicator for PCV-10 vaccine serotype, with serotype 6A excluded from NVT"

***Generating all bacteria from blood_culture
gen Bacpositive=(ct_code1!="" & ct_code1!="NGR" & specimen_type=="blood")|(ct_code2!="" & ct_code2!="NGR" & specimen_type=="blood")|(ct_code3!="" & ct_code3!="NGR" & specimen_type=="blood")
gsort serial -Bacpositive // Bacteria from blood
by serial: replace Bacpositive=Bacpositive[1] // fill up the rest of the records with the first value
label var Bacpositive "Positive bacteria Isolates from blood"

// indicator for SPN by culture in blood only
gen blood_spn1=(ct_code1=="SPN" & specimen_type=="blood") | (ct_code2=="SPN" & specimen_type=="blood") | (ct_code3=="SPN" & specimen_type=="blood")
gsort serial -blood_spn1 // blood spn
by serial: replace blood_spn1=blood_spn1[1] // fill up the rest of the records with the first value
label var blood_spn1 "SPN isolated in blood specimen"



// Define IPD
gen ipd=.
replace ipd=1 if pcv10!=. 
replace ipd=0 if pcv10==.
lab var ipd "Indicator for IPD "
*save cleaned_datafiles/all_admissions_tempkl,replace

***Agreed Definitions***
// Define pneumococcal meningitis
// CSF SPN Culture Pos OR (Blood SPN Culture Pos AND CSF WBC >50) OR (Blood SPN Culture Pos AND CSF:plasma glucose ratio <0.1
gen spnmen=0 if specimen_type=="csf"
replace spnmen=1 if csf_spn==1 | (blood_spn==1 & csfwbc>50 & csfwbc!=.)|(blood_spn==1 & glucratio<0.1)
gsort serial -spnmen
by serial: replace spnmen=spnmen[1] // fillup in all records of an admission 
label var spnmen "Pneumococcal meningitis as agreed by Anthony and Laura"

// Define HIB meningitis
gen hibmng=0 if specimen_type=="csf" | specimen_type=="blood"
replace hibmng=1 if ((isolate1== "Haemophilus influenzae" & specimen_type=="csf")   | (isolate2=="Haemophilus influenzae" & specimen_type=="csf")   | (isolate3=="Haemophilus influenzae" & specimen_type=="csf") | csfhib_antigen==1 )| ///
                    (((isolate1=="Haemophilus influenzae" & specimen_type=="blood") | (isolate2=="Haemophilus influenzae" & specimen_type=="blood") | (isolate3=="Haemophilus influenzae" & specimen_type=="blood")) & gramhib==1) 
gsort serial -hibmng
by serial: replace hibmng=hibmng[1] // fillup in all records of an admission 
label var hibmng "HIB meningitis"

// Define probable bacterial meningitis	
gen probmng=0 if specimen_type=="csf"
replace probmng= 1 if (csfwbc>=50 & csfwbc!=.) | glucratio<0.1 
gsort serial -probmng
by serial: replace probmng=probmng[1] // fillup in all records of an admission 
label define probmng 0 "No probable bacterial meningitis" 1 "Probable bacterial meningitis"
label values probmng probmng
label var probmng "Probable bacterial meningitis by Cowgill"

// Define Suspected meningitis
gen suspmng=0 if specimen_type=="csf"
replace suspmng=1 if (csfwbc>=10 & csfwbc<.) | glucratio<0.67
gsort serial -suspmng
by serial: replace suspmng=suspmng[1] // fillup in all records of an admission
label var suspmng "Suspected bacterial meningitis by Cowgill"

// Define Confirmed bacterial meningitis
gen confmng=0 if specimen_type=="csf" | specimen_type=="blood"
tab1 isolate1 isolate2 isolate3 if specimen_type=="csf" , m // all of these isolates can be considered pathogens
replace confmng=1 if ((pk_isolate1<. & specimen_type=="csf") | (pk_isolate3<. & specimen_type=="csf") | (pk_isolate3<. & specimen_type=="csf")) | spnmen==1 | hibmng==1
gsort serial -confmng
by serial: replace confmng = confmng[1] // fillup in all records of an admission
label var confmng "Confirmed bacterial meningitis by Cowgill"

// Define Suspected, probable or confirmed bacterial meningitis
gen bactmng=0 if specimen_type=="csf" | specimen_type=="blood"
replace bactmng=1 if suspmng==1 | probmng==1 | confmng==1
gsort serial -bactmng
by serial: replace bactmng = bactmng[1] // fillup in all records of an admission
label var bactmng "Suspected, probable or confirmed bacterial meningitis"

// Define pneumonia
tab cough,m
replace cough=lower(cough)
tab  breathing_difficulty,m
tab cough breathing_difficulty, m
tab cough_drn if  cough=="y", m
generate byte baseline_ip=0 if cough!="" | breathing_difficulty!="" // indicator for baseline pneumonia
replace baseline_ip = 1 if (cough=="y" | breathing_difficulty=="y")
	
inspect resp_rate
codebook resp_rate
generate byte fast_br_ip=0 if resp_rate!=.
replace fast_br_ip=1 if ((resp_rate>50 & agem>1 & agem<12) | (agey>0 & agey<5 & resp_rate>40)) & resp_rate!=. // calculate number with fast breathing
tab fast_br_ip baseline_ip, m

tab indrawing
replace indrawing=lower(indrawing) // indrawing

tab conscious_level,m
replace conscious_level="prostrate" if strpos(conscious_level, "prost")>0
replace conscious_level="lethargic" if strpos(conscious_level, "leth")>0
replace conscious_level="unconscious" if strpos(conscious_level, "unco")>0
replace conscious_level="normal" if inlist(conscious_level, "norma","n","5") | strpos(conscious_level,"norma")>0
	
gen hypoxia=1 if saturation_oxy<90 // hypoxia
replace hypoxia=0 if saturation_oxy>=90 & saturation_oxy!=.
tab hypoxia,m

gen pneumonia=0 if baseline_ip!=. // gen pneumonia using WHO criteria(Include cynosis, Vomitting everthing)
replace pneumonia=1 if baseline_ip==1 // mild pneumonia
replace pneumonia=2 if baseline_ip==1 & indrawing=="y"  // severe pneumonia
replace pneumonia=3 if (baseline_ip==1) & /// // very severe pneumonia, missing head_nodding, vomits_everything
(hypoxia==1 | lower(cyanosis) == "y" | lower(convulsion) == "y"| lower(head_nodding) == "y" | ///
lower(conscious_level) == "unconscious" | lower(conscious_level) == "lethargic" | ///
lower(unable_to_drink) == "y" | lower(vomits_everything) == "y" | lower(conscious_level) == "prostrate")	
label define pneumo_stat_ip 0 "None" 1 "Mild" 2 "Severe" 3 "Very severe" 	
label values pneumonia pneumo_stat_ip
label var pneumonia "Pneumonia severity"

// Define Pneumococcal Pneumonia
gen spnpna = 0
replace spnpna = 1 if (pneumonia==2 | pneumonia==3) & spn==1
tab spnpna spn,m
ren spnpna pneumococcal_pneumonia
label var pneumococcal_pneumonia "Pneumococcal pneumonia as agreed by Both Anthony and Laura"

// Define Clinical Meningitis
replace neck_stiffness=lower(neck_stiffness)
replace fontanelle_bulge=lower(fontanelle_bulge)
gen pros_concious=1 if inlist(conscious_level,"coma", "irritable","prostrate","unconscious")
gen meningitis=0 if fontanelle_bulge!="" | neck_stiffness!="" | (temp_axilla!=. & pros_concious!=.) | (temp_axilla!=. & convulsion!="")
replace meningitis =1 if fontanelle_bulge=="y"| neck_stiffness=="y"| (temp_axilla > 37.5 & temp_axilla!=. & pros_concious == 1) | (temp_axilla > 37.5 & temp_axilla!=. & convulsion=="y")
label var meningitis "Clinical meningitis"

// Define Pneumococcal bacteremia
// Define as pneumococcal disease with no meningitis or pneumonia e.g. GDB "non-pn, non-mng"
// Meningitis: exclude suspected, probable or confirmed bacterial based on defs above
// Pneumonia: exclude sev or very sev based on WHO case defs

gen pneumo_bacteremia=0 if ipd!=. & pneumonia!=. & bactmng!=.
replace pneumo_bacteremia=1 if ipd==1 & pneumonia<2 & bactmng!=1
label variable pneumo_bacteremia "Pneumococcal bacteremia"

// Define Severe Malnutrition: by WAZ
gen sevmn = (waz06 <=-3.0)
replace sevmn = . if waz06 == .
label define sevmn 0 "Non case" 1 "Case"
label values sevmn sevmn
label var sevmn "Severe malnutrition (WAZ)"

// Define moderate malnutrition: by WAZ
gen modmn = waz06 >-3.0 & waz06<=-2.0 
replace modmn = . if waz06 == .
label define modmn 0 "Non case" 1 "Case"
label values modmn modmn
label var modmn "Moderate malnutrition (WAZ)"	

// Define moderate or severe malnutrition: by WAZ;i fred Will use this for my analysis
gen msmn = 0 if sevmn!=. | modmn!=.
replace msmn=1 if sevmn==1 | modmn==1
label values msmn modmn
label var msmn "Moderate or severe malnutrition (WAZ)"	
tab yoa msmn if resident==1 & agem<60,m 

// Define Severe Malnutrition: by MUAC
generate sevmnmuac=0 if muac!=. & agem<60
replace sevmnmuac =1 if muac<=7.5 & agem<6  & sevmnmuac!=.
replace sevmnmuac =1 if muac<=11.5 & agem>=6 & agem <60 & sevmnmuac!=. 
label values sevmnmuac modmn
label var sevmnmuac "Severe malnutrition (MUAC)"  
tab yoa sevmnmuac if resident==1 & agem<60,m // 1999 has alot of missing muac

// Define Severe Malnutrition: by WHZ
gen sevmnwhz = whz06 <=-3.0
replace sevmnwhz = . if whz06 == .
label values sevmnwhz modmn
label var sevmnwhz "Severe malnutrition (WHZ)"
tab yoa sevmnwhz if resident==1 & agem<60,m  

// Syndromes: hierarchical definition
gen synd=0 & ipd==1
replace synd=1 if spnmen==1 & synd!=.
replace synd=2 if pneumonia==3 & synd!=.
replace synd=3 if pneumonia==2 & synd!=.
replace synd=4 if pneumo_bacteremia==1 & synd!=.

tab synd
lab def synd 0 "other" 4 "bacteraemia" 3 "severe pneumonia"   2 "very severe pneumonia"   1 "meningitis"
lab val synd synd
label var synd "Syndromes: hierarchical definition"
tab synd ipd,m



// Define Eligibility for blood culture: Exclude cases with accidents/trauma/elective surgery
tab  diagnosis_1_disch 
replace  diagnosis_1_disch =lower(diagnosis_1_disch)
gen strep=1 if strpos(diagnosis_1_disch, "accidents")>0|strpos(diagnosis_1_disch,"elective surgery")>0 | ///
strpos(diagnosis_1_disch, "trauma")>0|strpos(diagnosis_1_disch, "burns")>0 |strpos(diagnosis_1_disch, "epilepsy")>0 ///
| strpos(diagnosis_1_disch, "snake bite")>0|  ///
strpos(diagnosis_1_disch, "poisoning")>0 | strpos(diagnosis_1_disch, "developmental")>0 |strpos(diagnosis_1_disch, "fracture")>0
gen eligible="no" if strep==1
replace eligible="yes" if eligible==""
tab eligible specimen_type,m
label var eligible "Eligible for blood culture"

// Drop duplicate records by serial
sort serial fk_specimen
duplicates drop serial, force
drop fk_specimen pk_specimen_id specimen_type
compress 
replace resident=2 if resident==. & ipd==1

**cleaning diagnosis discharge
rename diagnosis_1_disch diagnosis1_disch
rename diagnosis_2_disch diagnosis2_disch
tab diagnosis1_disch, sort	
	replace diagnosis1_disch=lower(diagnosis1_disch)
	replace diagnosis1_disch=trim(diagnosis1_disch)
	tab diagnosis1_disch, sort	
	replace diagnosis1_disch="anaemia" if strpos(diagnosis1_disch,"anaemia")
	replace diagnosis1_disch="dysentery" if strpos(diagnosis1_disch,"dysentry")
	replace diagnosis1_disch="malnutrition" if strpos(diagnosis1_disch,"malnutrition")
	replace diagnosis1_disch="poisoning - kerosene" if strpos(diagnosis1_disch,"parrafin")|strpos(diagnosis1_disch,"kerosine")
	replace diagnosis1_disch="cellulitis/pyomyositis/abscess" if strpos(diagnosis1_disch,"cecellulitis/pyomyositis/abscess")
	replace diagnosis1_disch="viral hapatitis" if strpos(diagnosis1_disch,"patitis")
	replace diagnosis1_disch="submandibular abscess" if strpos(diagnosis1_disch,"submandibula")
	replace diagnosis1_disch="prematurity" if diagnosis1_disch=="prematurity"|diagnosis1_disch=="preterm"
	replace diagnosis1_disch="poisoning wild fruits" if strpos(diagnosis1_disch, "fruits")
	replace diagnosis1_disch="infected post circumcised wound" if diagnosis1_disch=="infected post circumsed wound"
	replace diagnosis1_disch="lrti" if strpos(diagnosis1_disch,"lrti")
	replace diagnosis1_disch="unclassified disease" if strpos(diagnosis1_disch,"nclassified")
	replace diagnosis1_disch="nephritic syndrome" if diagnosis1_disch=="nephrotic syndrome"
	replace diagnosis1_disch="sickle  cell  disease" if diagnosis1_disch=="sickle cell disease"|diagnosis1_disch=="sickle cell trait"|diagnosis1_disch=="other sickle cell disorders"
	replace diagnosis1_disch="cellulitis/pyomyositis/abscess" if strpos(diagnosis1_disch,"cellulitis")
	order diagnosis1_disch
	
	tab diagnosis2_disch, sort
		replace diagnosis2_disch=lower(diagnosis2_disch)
	replace diagnosis2_disch=trim(diagnosis2_disch)
	tab diagnosis2_disch, sort	
	replace diagnosis2_disch="anaemia" if strpos(diagnosis2_disch,"anaemia")
	replace diagnosis2_disch="dysentery" if strpos(diagnosis2_disch,"dysentry")
	replace diagnosis2_disch="malnutrition" if strpos(diagnosis2_disch,"malnutrition")
	replace diagnosis2_disch="poisoning - kerosene" if strpos(diagnosis2_disch,"parrafin")|strpos(diagnosis2_disch,"kerosine")
	replace diagnosis2_disch="cellulitis/pyomyositis/abscess" if strpos(diagnosis2_disch,"cecellulitis/pyomyositis/abscess")
	replace diagnosis2_disch="viral hapatitis" if strpos(diagnosis2_disch,"patitis")
	replace diagnosis2_disch="submandibular abscess" if strpos(diagnosis2_disch,"submandibula")
	replace diagnosis2_disch="prematurity" if diagnosis2_disch=="prematurity"|diagnosis2_disch=="preterm"
	replace diagnosis2_disch="poisoning wild fruits" if strpos(diagnosis2_disch, "fruits")
	replace diagnosis2_disch="infected post circumcised wound" if diagnosis2_disch=="infected post circumsed wound"
	replace diagnosis2_disch="lrti" if strpos(diagnosis2_disch,"lrti")
	replace diagnosis2_disch="unclassified disease" if strpos(diagnosis2_disch,"nclassified")
	replace diagnosis2_disch="nephritic syndrome" if diagnosis2_disch=="nephrotic syndrome"
	replace diagnosis2_disch="sickle  cell  disease" if diagnosis2_disch=="sickle cell disease"|diagnosis2_disch=="sickle cell trait"|diagnosis2_disch=="other sickle cell disorders"
	replace diagnosis2_disch="cellulitis/pyomyositis/abscess" if strpos(diagnosis2_disch,"cellulitis")
**Keeping variables
**Malaria
	codebook malaria
	**merge with malarai_Anthony data to clear Malaria problem once and for all
	*br mps mps_100wbc mps_500rbc gametocytes pigment wbc rbc datah postive_fromSharefolderOnly
	preserve
	use "C:\Users\THINKPAD\Downloads\malaria_Anthony.dta", clear
	drop fk_specimen pk_specimen_id datah postive_fromSharefolderOnly
	keep serial mps gametocytes pigment mps_500rbc mps_100wbc
	generate mal2mps=1 if mps=="positive"
	replace mal2mps=0 if mps=="negative"
	replace mal2mps=9 if mps=="not requested"|missing(mps)
	save data/MalariaAnto,replace
	restore
	
	**check malaria count
	tab malaria if yoa>2016,m
	
	**merge admissions and updated Malaria Anthony
	merge 1:1 serial using data/MalariaAnto
	drop if _merge==2
	
	replace malaria=mal2mps if (!missing(mal2mps) & missing(malaria) & _merge==3 & yoa>2016 & !missing(yoa))
	
	
	drop _merge mal2mps
	*gen Malaria_mps=malaria
	drop mps
	**check malaria count
	tab malaria if yoa>2016,m
	ren malaria mps
	label variable mps "final malaria results"
	**Now I have may variable(malaria) and others to send to Lilian for formality mps gametocytes pigment mps_500rbc mps_100wbc
	
	**HIV results.The variable is hiv
	preserve
	keep serial fk_person diagnosis1_disch diagnosis2_disch
	order serial fk_person diagnosis1_disch diagnosis2_disch
	describe, short
	*export excel using "C:\Users\THINKPAD\Desktop\DesktopDamp\MSC_STATISTICSNOTES\Personal_lerning notes\MSC_Project\documents\referencesTo_ReadForMethodFavouring\data\diagnosis_Discharge.xls", sheetreplace firstrow(variables)
	export delimited using "C:\Users\THINKPAD\Desktop\DesktopDamp\MSC_STATISTICSNOTES\Personal_lerning notes\MSC_Project\documents\referencesTo_ReadForMethodFavouring\data\diagnosis_Discharge.csv", replace
	restore
	
	
	**lets compare chatGpt and mine
	* Keep originals, create cleaned versions
gen strL dx1 = lower(itrim(strtrim(diagnosis1_disch)))
gen strL dx2 = lower(itrim(strtrim(diagnosis2_disch)))
**chatGpt
gen byte measleschatGpt = .
replace measleschatGpt = 0 if !missing(dx1) | !missing(dx2)
replace measleschatGpt = 1 if dx1=="measles" | dx2=="measles"

label define yn 0 "non-case" 1 "case", replace
label values measleschatGpt yn

tab measleschatGpt, missing

	
	
	
	
	**Measles
	tab diagnosis1_disch,missing
	tab diagnosis1_disch if strpos(diagnosis1_disch,"mea")|strpos(diagnosis1_disch,"sles")
	gen Measles=1 if strpos(diagnosis1_disch,"mea")|strpos(diagnosis1_disch,"sles")
	
	
	tab diagnosis2_disch,missing
	tab diagnosis2_disch if strpos(diagnosis2_disch,"mea")|strpos(diagnosis2_disch,"sles")
	tab diagnosis2_disch,missing sort
	tab diagnosis2_disch if ///
    (strpos(diagnosis2_disch,"mea") | strpos(diagnosis2_disch,"sles")) ///
    & strpos(diagnosis1_disch,"measles") == 0
	
	replace Measles=1 if (strpos(diagnosis2_disch,"mea")|strpos(diagnosis2_disch,"sles")) & Measles!=1
	
	
	replace Measles=0 if (!missing(diagnosis2_disch)|!missing(diagnosis1_disch)) & missing(Measles)
	
	label define me 1 "case" 0 "non-case"
	label values Measles me
	order diagnosis1_disch diagnosis2_disch Measles
	gsort -Measles
	tab Measles,m
	
	**neonatal sepsi
	tab diagnosis1_disch if strpos(diagnosis1_disch,"sepsis")|strpos(diagnosis1_disch,"sep")
	*
	* Clean text (recommended)
replace diagnosis1_disch = lower(trim(diagnosis1_disch))

* Create neonatal sepsis indicator
*gen neonatal_sepsi = 0

* Replace =1 if diagnosis matches ANY sepsis string AND age <= 28 days
* Start clean
*replace diagnosis1_disch = lower(trim(diagnosis1_disch))

* Create neonatal sepsis indicator
gen neonatal_sepsi = 0

* Group 1
replace neonatal_sepsi = 1 if aged <= 28 & ///
    inlist(diagnosis1_disch, ///
        "neonatal sepsis", ///
        "nneonatal sepsis", ///
        "other septicaemia/septicaemia", ///
        "sepsis /iss", ///
        "septicemia")

* Group 2
replace neonatal_sepsi = 1 if aged <= 28 & ///
    inlist(diagnosis1_disch, ///
        "septicaemia/sepsis", ///
        "streptococcal septicaemia")

* Group 3
replace neonatal_sepsi = 1 if aged <= 28 & ///
    inlist(diagnosis1_disch, ///
        "septic wound", ///
        "septic arthritis", ///
        "septic arthritis right knee", ///
        "septic arthrits", ///
        "septic bedsores")

* Group 4
replace neonatal_sepsi = 1 if aged <= 28 & ///
    inlist(diagnosis1_disch, ///
        "septic glossitis", ///
        "septic i/site")

* Group 5
replace neonatal_sepsi = 1 if aged <= 28 & ///
    inlist(diagnosis1_disch, ///
        "septic incision", ///
        "sau sepsis")

		
		
		
	tab diagnosis2_disch,m	
	tab diagnosis2_disch if strpos(diagnosis2_disch,"sepsis")|strpos(diagnosis2_disch,"sep")
	*
	* Clean text (recommended)
replace diagnosis2_disch = lower(trim(diagnosis2_disch))

	* Add neonatal sepsis based on diagnosis2_disch + age <= 28 days

* Group 1: explicit neonatal + general septicemia/sepsis
replace neonatal_sepsi = 1 if aged <= 28 & ///
    inlist(diagnosis2_disch, ///
        "neonatal sepsis", ///
        "other septicaemia/septicaemia", ///
        "presumed sepsis", ///
        "sepsis", ///
        "septicaemia/sepsis", ///
        "septicemia")

* Group 2: localized septic conditions that can still be serious neonatal infections
replace neonatal_sepsi = 1 if aged <= 28 & ///
    inlist(diagnosis2_disch, ///
        "septic wound", ///
        "septic arthritis", ///
        "septic arthritis right knee", ///
        "septic arthrits", ///
        "septic bedsores", ///
        "septic glossitis", ///
        "septic i/site", ///
        "septic incision")

* Group 3: streptococcal septicemia + mixed label
replace neonatal_sepsi = 1 if aged <= 28 & ///
    inlist(diagnosis2_disch, ///
        "streptococcal septicaemia", ///
        "severe jaundice and neonatal sepsis")

		
	replace neonatal_sepsi=. if missing(diagnosis2_disch) & missing(diagnosis1_disch)	
	
	label values neonatal_sepsi me
	
	
	
	**Congenital abnomalities/defects/disorder
	* Generate birth_defects indicator (0 = no birth defect, 1 = congenital abnormality)
gen byte birth_defects = 0

replace birth_defects = 1 if inlist(lower(diagnosis1_disch), ///
    "congenital abnormality - other", ///
    "congenital anom", ///
    "congenital bowel obstruction", ///
    "congenital cns disease", ///
    "congenital syphillis", ///
    "heart disease - congenital", ///
    "club foot.")

	
	* Standardize text
replace diagnosis2_disch = lower(trim(diagnosis2_disch))

* Add birth defects based on diagnosis2 (structural congenital abnormalities)
* Mark additional birth defects from diagnosis2_disch
replace birth_defects = 1 if inlist(lower(diagnosis2_disch), ///
    "congenital abnormality - other", ///
    "congenital anom", ///
    "congenital cns disease", ///
    "congenital syphillis", ///
    "heart disease - congenital", ///
    "chromosomal  abnormality", ///
    "avsd")


	
	replace birth_defects=. if missing(diagnosis2_disch) & missing(diagnosis1_disch)	
	
	label values birth_defects me
	
	
	
	
	**Birth Asphyxia
	*--------------------------------------------------*
* 1. Create Birth_Asphyxia indicator              *
*--------------------------------------------------*

* Start with 0 = no birth asphyxia
generate byte birth_asphyxia = 0

* Flag any record where either diagnosis contains "asphyxia"
replace birth_asphyxia = 1 if ///
    (!missing(diagnosis1_disch) & strpos(lower(diagnosis1_disch), "asphyxia")) | ///
    (!missing(diagnosis2_disch) & strpos(lower(diagnosis2_disch), "asphyxia"))

*--------------------------------------------------*
* 2. Set to missing if both diagnosis fields empty *
*--------------------------------------------------*

replace birth_asphyxia = . if missing(diagnosis1_disch) & missing(diagnosis2_disch)

*--------------------------------------------------*
* 3. Label the variable                            *
*--------------------------------------------------*

*label define asphyx 0 "no birth asphyxia" 1 "birth asphyxia", replace
label values birth_asphyxia me

*--------------------------------------------------*
* 4. Quick checks                                  *
*--------------------------------------------------*

tab birth_asphyxia, missing

tab diagnosis1_disch if birth_asphyxia == 1
tab diagnosis2_disch if birth_asphyxia == 1

**Keep the required variables
keep if inrange(doa, d(01apr2011), d(31dec2019))
br agey if agey<=15
keep if agey<=15
preserve
keep serial mps gametocytes pigment mps_500rbc mps_100wbc hivresult Measles neonatal_sepsi dob doa dod sex weight height birth_defects birth_asphyxia fk_person oedema resident source 	
replace fk_person=fk_person+10
save data/AmissionMcsProjectFred,replace
restore
replace fk_person=fk_person+10
replace fk_person=(fk_person-10)
keep serial mps gametocytes pigment mps_500rbc mps_100wbc hivresult Measles neonatal_sepsi dob doa dod sex weight height birth_defects birth_asphyxia fk_person oedema resident msmn diagnosis1_disch diagnosis2_disch agem agey aged yoa source
sort fk_person doa
order serial fk_person doa
**check dups
duplicates tag fk_person doa,gen(dups)
br if dups>0
gsort -serial fk_person doa
duplicates drop fk_person doa,force
save data/AmissionMcsProjectFredMainAnalysis,replace	

*confounders:Malaria Inc,HIV Inc,measles Inc,neonatal sepsis Inc,congenital abnormalities Inc,birth Asphyxia
*Inc,Malnutrition In,strike,time,low birth weight Inc, and seasonality.
**UPTO NOW STRIKE PERIODS HAS NOT BEEN REMOVED FROM THIS DATA.PLEASE NOTE THAT!!!
codebook mps 
* 1=Positive for mps
codebook hivresult
* 1=Positive for HIV
codebook Measles
* 1=case for measles
codebook neonatal_sepsi
* 1=case for neonatal sepsis

codebook birth_defects
* 1=case for birth defect/congenital abnormalities/birth disorde

codebook birth_asphyxia
*1=case for birth asphia
 codebook  msmn
 *1=case for msmn
 **Will create low birthweight,strike variable,time(year), and seasonality no worries
 
 *Lets create data with strike periods removed then we will adjust mid-year population estimate by letting mid-year estimate be the proportion of days with no strike times the mid-year population. Peron time in month estimated as mid-year population divided by 12
 preserve
 count
 drop if doa>=date("01/09/2009", "DMY") & doa<=date("31/03/2010", "DMY") & source=="adult"  // 1979
drop if doa>=date("05/12/2011", "DMY") & doa<=date("12/12/2011", "DMY") // 76
drop if doa>=date("01/03/2012", "DMY") & doa<=date("15/03/2012", "DMY") // 170
drop if doa>=date("13/09/2012", "DMY") & doa<=date("04/11/2012", "DMY") // 802
drop if doa>=date("02/12/2012", "DMY") & doa<=date("22/12/2012", "DMY") // 82
drop if doa>=date("23/12/2012", "DMY") & doa<=date("13/01/2013", "DMY") & source=="adult" // 0
drop if doa>=date("16/01/2013", "DMY") & doa<=date("11/02/2013", "DMY") // 77
drop if doa>=date("10/12/2013", "DMY") & doa<=date("21/12/2013", "DMY") // 18
drop if doa>=date("05/12/2016", "DMY") & doa<=date("31/12/2016", "DMY") // 128
drop if doa>=date("01/01/2017", "DMY") & doa<=date("15/03/2017", "DMY")
drop if doa>=date("05/06/2017", "DMY") & doa<=date("02/11/2017", "DMY")
drop if doa>=date("24/03/2020", "DMY") & doa<=date("26/03/2020", "DMY") & source=="adult"
drop if doa>=date("08/12/2020", "DMY") & doa<=date("02/03/2021", "DMY") & source=="adult"
 save data/AmissionMcsProjectFredMainAnalysis_StrikePeriodsRemoved,replace
 count
 
 **tab msmn if agem < 12 & resident == 1 & (doa >= td(01apr2011) & doa <= td(31dec2019)),m
 restore
 

 **Lets do cleaning of mid_year population estimate
 // population data

insheet using "C:\Users\THINKPAD\Downloads\PCVIS Analysis-20251117T080637Z-1-001\PCVIS Analysis\Analysis_Files\raw_datafiles\Pop_estimate_full.csv", clear
keep if year<2023
gen midpop_12_23m=midpop_u2yrs - midpop_u2m - midpop_2_11m
gen midpop_2_23m=midpop_u2yrs - midpop_u2m
gen midpop_under1=midpop_u2m+midpop_2_11m

reshape long midpop, i(year) j(agecat) string
replace agecat="1" if agecat=="_u2m"
replace agecat="2" if agecat=="_2_11m"
replace agecat="3" if agecat=="_12_23m"
replace agecat="4" if agecat=="_2_4"
replace agecat="5" if agecat=="_5_14"
replace agecat="6" if agecat=="_15_49"
replace agecat="7" if agecat=="_50_64"
replace agecat="8" if agecat=="_over64"
replace agecat="9" if agecat=="_2_23m"
replace agecat="10" if agecat=="_u2yrs"
replace agecat="11" if agecat=="_u5"
replace agecat="12" if agecat=="_over15"
replace agecat="13" if agecat=="_over50"
replace agecat="14" if agecat=="_under1"
destring agecat, force replace
drop if agecat==.
ren year yoa
ren midpop midyrpop
**saving for methods **14 is the under1 year(infant)
br if agecat==14
*save cleaned_datafiles/temp14FREDPOP, replace
save "C:/Users/THINKPAD/Desktop/DesktopDamp/MSC_STATISTICSNOTES/Personal_lerning notes/MSC_Project/documents/referencesTo_ReadForMethodFavouring/data/temp14FREDPOP",replace
replace midyrpop= midyrpop-(122/365)*midyrpop if yoa==2009 & inlist(agecat,6,7,8,12,13) 
replace midyrpop= midyrpop-(90/365)*midyrpop if yoa==2010 & inlist(agecat,6,7,8,12,13)
replace midyrpop= midyrpop-(8/365)*midyrpop if yoa==2011 
replace midyrpop= midyrpop-(58/366)*midyrpop if yoa==2012
replace midyrpop= midyrpop-(8/366)*midyrpop if yoa==2012 & inlist(agecat,6,7,8,12,13)
replace midyrpop= midyrpop-(13/365)*midyrpop if yoa==2013 & inlist(agecat,6,7,8,12,13)
replace midyrpop= midyrpop-(39/365)*midyrpop if yoa==2013
replace midyrpop= midyrpop-(27/366)*midyrpop if yoa==2016
replace midyrpop= midyrpop-(223/365)*midyrpop if yoa==2017
**ward admission stopped 24th-26th march 2020
replace midyrpop= midyrpop-(27/366)*midyrpop if yoa==2020 & inlist(agecat,6,7,8,12,13)
***strike period 1st jan to 2nd march 2021
replace midyrpop= midyrpop-(61/365)*midyrpop if yoa==2021 & inlist(agecat,6,7,8,12,13)
*save cleaned_datafiles/temp14, replace
 save "C:/Users/THINKPAD/Desktop/DesktopDamp/MSC_STATISTICSNOTES/Personal_lerning notes/MSC_Project/documents/referencesTo_ReadForMethodFavouring/data/temp14fremat",replace
 preserve
 keep if agecat==14|agecat==5
 keep yoa midyrpop agecat
 
 reshape wide midyrpop,i(yoa) j(agecat)
 ren midyrpop14 midyrpopunder1
 ren midyrpop5  midyrpop5to14
  save "C:/Users/THINKPAD/Desktop/DesktopDamp/MSC_STATISTICSNOTES/Personal_lerning notes/MSC_Project/documents/referencesTo_ReadForMethodFavouring/data/temp14frematunder1",replace
 **Lets try and see Ben's mid_pop
 restore
 
 
 
 
 
 
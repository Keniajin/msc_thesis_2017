*********************************************************************************************************************************************************
****************** Script ffor analysis**************************************************************************
***********************   analysing admission data with surveilance data   *******************************************************************************************
****************************    Author: KMwai Feb 2016  ************************************************************************************************
***************** Data from KWTRP Kilifi (David Amadi Data Manager) . *****************************************************************************************
	clear all
	set more off
	cd "E:\Wits\Project"
	cap log close
	local time2=lower(word(c(current_time),1))
	local time=subinstr("`time2'", ":", "~",. )
	local tod=lower(word(c(current_date),1)+word(c(current_date),2)+word(c(current_date),3))
	local name = "logs\"+"masters_analysis"+"`tod'"+"~"+"`time'"
	log using `name',text replace	
//load haematology data from 1998-2008	from haematology folder
	***************set seed for the same random select any time******
	set seed 1221223

	****insheet data from R EVI
	use data/all_data.dta, clear
*****create categories of time to admission
	gen time_to_readmit_cat=0 if time_to_readmission==0
	replace time_to_readmit_cat=1 if time_to_readmission>=1 &  time_to_readmission<=30
	replace time_to_readmit_cat=2 if time_to_readmission>30 &  time_to_readmission<=90
	replace time_to_readmit_cat=3 if time_to_readmission>90 &  time_to_readmission<=180
	replace time_to_readmit_cat=4 if time_to_readmission>180 &  time_to_readmission<=365
	replace time_to_readmit_cat=5 if time_to_readmission>365 & time_to_readmission<=730
	replace time_to_readmit_cat=6 if time_to_readmission>730 & time_to_readmission<=1095
	replace time_to_readmit_cat=7 if time_to_readmission>1095 & time_to_readmission!=.
	
	cap label define ttr 0 "No readmission" 1 "1-30days" 2 "1-3 months" 3 "3 -6 months" 4 "6mnths - 1yr"  5 "1yr-2yrs" ///
	6 "2yr-3yrs" 7 ">3yrs"
	label values time_to_readmit_cat ttr

*****generate number of admissions per quarter
	gen quarters = qofd(ndad )
	format quarters %tq
	bys quarters fk_person : gen admission_qtr=_n if fk_person!=.
	
	bys quarters fk_person : gen admission_qtr_mal=_n if fk_person!=. & malnutKid==1
	replace admission_qtr_mal=0 if fk_person!=. & malnutKid==0
****generatte the number of malnutrition events per person
	gen malnut2=1 if malnutKid==0
	replace malnut2=0 if malnutKid==1

	bys  fk_person (ndad) : gen cases_count=sum(malnutKid) if fk_person!=. & malnutKid==1
	bys  fk_person (ndad) : gen cases_count_time=ndad[_n] - ndad[1] if fk_person!=. & malnutKid==1
	
	****cumulitive count for malnutrition cases
	gen cumulitive_count=cases_count
	bys  fk_person (ndad) : replace cumulitive_count=0 if ndad==ndad[1] & cumulitive_count==. & fk_person!=.
	bys  fk_person (ndad) : egen case_count_max=max(cases_count)
	
	****why 16
	forval i =  1/16 {
	bys  fk_person (ndad) : replace cumulitive_count=cases_count[_n-`i'] if cumulitive_count==. & fk_person!=. 
	}
	
	
	forval i =  1/16 {
	bys  fk_person (ndad) : replace cumulitive_count=cumulitive_count[_n-`i'] if cumulitive_count==. & fk_person!=.
	}
	
	***replace the vlaues after the maximum count to missing
	replace cumulitive_count=. if cases_count==. & cumulitive_count==case_count_max
	****count for control cases
	bys  fk_person (ndad) : gen control_count=sum(malnut2) if fk_person!=. & malnut2==1
	bys  fk_person (ndad) : gen control_count_time=ndad[_n] - ndad[1] if fk_person!=. & malnut2==1
	
	******
	gen cumulitive_time=control_count_time
	replace cumulitive_time=cases_count_time if cumulitive_time==. 
	replace cumulitive_time=cumulitive_time+1
	//replace admission_qtr_mal=0 if fk_person!=. & malnutKid==0
	// keep ndad fk_person serialno cumulitive_time control_count_time cases_count_time malnutKid admission_qtr	
*****generate number of admissions per quarter
	gen mnths = month(ndad )
	gen timem = mofd(ndad)
	format timem %tm
	//bys quarters fk_person : gen admission_qtr=_n if fk_person!=.
	
**************************generating collapse datas****		
******************save collapse*****************
	preserve
	collapse  (count) malnutT=malnutKid   ///
	(sum) mal_case_P=malnutKid mal_con_P=malnut2  severeD_P=severe_disease2 ///
	(mean) meanAge=nagem  meanEVI=EVI_VALUE if severe_disease2==1 	, by(quarters )
	save "data/colapse_seve.dta", replace
	restore
******************end save collapse*****************
	preserve
	collapse  (count) totalMalnut=malnutKid   (sum) mal_case_T=malnutKid severeD_T=severe_disease2 ///
	(mean) meanAge=nagem  meanEVI=EVI_VALUE , by(quarters )
	gen  mal_con_T=totalMalnut-mal_case_T
	merge 1:1 quarters using "data/colapse_seve.dta"
	drop _merge
	save "data/colapse_morbidity_quarters", replace
	restore
	

****save for monthly trend***
	******************save collapse*****************
	preserve
	collapse  (count) malnutT=malnutKid   ///
	(sum) mal_case_P=malnutKid mal_con_P=malnut2 severeD_P=severe_disease2  	///
	(mean) meanAge=nagem  meanEVI=EVI_VALUE if severe_disease2==1, by(timem )
	save "data/colapse_seve.dta", replace
	restore
******************end save collapse*****************
	preserve
	collapse  (count) totalMalnut=malnutKid   (sum) mal_con_T2=malnut2 mal_case_T=malnutKid severeD_T=severe_disease2 ///
	(mean) meanAge=nagem  meanEVI=EVI_VALUE , by(timem )
	gen  mal_con_T=totalMalnut-mal_case_T
	merge 1:1 timem using "data/colapse_seve.dta"
	drop _merge
	save "data/colapse_morbidity_mnths", replace
	restore
	
	
****save mortality information
	preserve
	collapse  (count) totalDiedC=noutcome   (sum) mal_con_Die=noutcome  ///
	(mean) meanAge_C=nagem  meanEVI_C=EVI_VALUE if malnut2==1, by(timem ) 
	save "data/colapse_morbidity_mnths_died_c", replace
	restore
	
****save mortality information
	preserve
	collapse  (count) totalDied=noutcome   (sum) mal_case_Die=noutcome  ///
	(mean) meanAge=nagem  meanEVI=EVI_VALUE if malnutKid==1, by(timem )
	merge 1:1 timem using "data/colapse_morbidity_mnths_died_c.dta"
	save "data/colapse_morbidity_mnths_died", replace
	restore
	
	
****save for monthly trend one***
	******************save collapse*****************
	preserve
	collapse  (count) malnutT=malnutKid   ///
	(sum) mal_case_P=malnutKid mal_con_P=malnut2 if severe_disease2==1 	, by(mnths )
	save "data/colapse_seve.dta", replace
	restore
******************end save collapse*****************
	preserve
	collapse  (count) totalMalnut=malnutKid   (sum) mal_case_T=malnutKid , by(mnths )
	gen  mal_con_T=totalMalnut-mal_case_T
	merge 1:1 mnths using "data/colapse_seve.dta"
	drop _merge
	save "data/colapse_morbidity_mnths_one", replace
	restore

******************save collapse file for SATSCAN*****************
	preserve
	collapse  (count) cases_cnt=malnutKid   (sum) case_T=malnutKid , by(yr sublocation)
	gen  con_T=cases_cnt-case_T
	outsheet using "data/malnut_admission_count.csv", replace c
	restore
	***take the data to R to save for the HOTSPOTS
	
	
******************save collapse file for logistic SATSCAN*****************
	preserve
	collapse  (count) cases_cnt=malnutKid   (sum) case_T=malnutKid if noutcome==1, by(yr sublocation) 
	outsheet using "data/malnut_died_count.csv", replace c
	restore
	
	
	******************save collapse file for logistic SATSCAN*****************
	preserve
	collapse  (count) con_cnt=malnut2   (sum) con_T=malnut2 if noutcome==1, by(yr sublocation) 
	outsheet using "data/non_malnut_died_count.csv", replace c
	restore
	***take the data to R to save for the HOTSPOTS
	
	
******************end save collapse*****************
******#########################*********************
**** time between the last and the first admision
	gen first_to_last=0	
	bysort fk_person (ndad): replace first_to_last= ndad[_N] - ndad[1] if  fk_person!=.

**** time between the last and the first admision by malnutrition
	gen first_to_last_case=.	
	bysort malnutKid fk_person (ndad) : replace first_to_last_case= ndad[_n] - ndad[1] if  fk_person!=.	& malnutKid==1	
	
	gen first_to_last_control=.	
	bysort malnutKid fk_person (ndad )  : replace first_to_last_control= ndad[_n] - ndad[1] if  fk_person!=.	& malnutKid==0
****count the total admissions
	bys fk_person (ndad) :gen count_adm=_n if fk_person!=.
	bys fk_person (ndad) :gen total_admission=_N if fk_person!=.
	
****
****generate the sublocations for numeric****
	encode sublocation, generate(sublocs)
	gen gender=1 if nsex==1
	replace gender=0 if nsex==2
  
	gen nbirthweight=birth_weight /100
	gen nagem_int=int(nagem)
	
	merge m:m sublocation using "data/sublocAdj_ids.dta", gen(_merge23)
	drop _merge23
	
	
*****drop the missing cumulitive count
	preserve 
	//837 obsevations deleted
	drop if cumulitive_count==.
	
	///14,755 admissions 14458
	duplicates tag fk_person , gen(tag)
	tab tag
	drop if tag==0
	tab cumulitive_count
	hist cumulitive_count
	tab count_adm
	merge m:m sublocation using "data/sublocAdj_ids.dta", gen(_merge23)
	save "data\morbidity.dta", replace
	restore

****keeep the variables of interest****
	preserve
	keep malnutKid serialno ndad fk_person  cumulitive_count cumulitive_time  ///
	first_to_last first_to_last_case first_to_last_control count_adm cases_count
	save "data/count_data", replace
	restore
	
	compress
	gen whz_all=whz06
	replace whz_all =_zbfa if whz_all==.
	save "data\mortality.dta" , replace
	outsheet using "data\mortality.csv" , replace c

******mortality outcome

////categorical variables for the outcome==1 (mortality)
	tabout  nsex  severe_disease      ///
	severe_anaemia    hypoglycaemia         ///
	malaria    ndiarrhoea  menegitis lrti  ///
	gastroenteritis ntransfused  bc_pos csf_pos ///
	malnutKid  using analysis/outputs/table1_morta.txt if noutcome==1, c(freq col) f(0c 1p) stats(chi2) rep
	
	tabout  nsex  severe_disease      ///
	severe_anaemia    hypoglycaemia         ///
	malaria    ndiarrhoea  menegitis lrti  ///
	gastroenteritis ntransfused  bc_pos csf_pos ///
	malnutKid  using analysis/outputs/table1_morta2.txt if noutcome==0, c(freq col) f(0c 1p) stats(chi2) rep

	
	 ///numeric variables for the outcome==1
	tabstat  nagem admdays  no_admissions total_admission nweight  nheight nmuac   haz06 waz06 whz06 ///
	EVI_VALUE  if noutcome==1 , by(malnutKid) ///
	stat(mean sd p50 min max) nototal long col(stat)	
	
		tabstat  nagem admdays  no_admissions total_admission nweight  nheight nmuac   haz06 waz06 whz06 ///
	EVI_VALUE  if noutcome==0 , by(malnutKid) ///
	stat(mean sd p50 min max) nototal long col(stat)	

  
	outsheet using "data/all_stata.csv" , c replace

  //Zero-inflated Negative Binomial Regression - Negative binomial regression does better with over dispersed data, i.e. variance much larger than the mean. 

 
****Model adding birth weigth
	/* ****

	 **model with birthweight the numbers go down
	 *mepoisson cumulitive_count EVI_VALUE   nagem i.nsex  nbirthweight ///
	 *|| sublocation: || yr: , irr cov(unstructured)
	 *estat ic
	 
	 xtmepoisson cumulitive_count EVI_VALUE   nagem i.nsex i.severe_disease nbirthweight ///
	 || sublocs: ,   exposure(cumulitive_time) cov(unstructured) nolog
	estat ic
	est store M5
	
	 xtmepoisson cumulitive_count EVI_VALUE   nagem i.nsex i.severe_disease nbirthweight ///
	 || fk_person: ,   exposure(cumulitive_time) cov(unstructured) nolog
	estat ic
	est store M5b

	****dont like this because it uses the pk_serial as the grouping variable
	 xtpoisson cumulitive_count  EVI_VALUE   nagem i.nsex i.severe_disease ///
	 nbirthweight , exposure(cumulitive_time)  normal nolog
	 estat ic
	est store M6*/
	
	
***negative binomial variance should be greater than the mean
***number of morbidity variance is less than mean
	 *xtnbreg no_morbidity i.malnutKid EVI_VALUE nagem nsex, exposure(cumulitive_time)  nolog
	 
	 gen ageyrs=nagem/12
	 *menbreg no_morbidity i.malnutKid EVI_VALUE ageyrs nnsex, exposure(cumulitive_time) irr || sublocation: || yr:


*#################################################################################****
****modelling the mortality variable
****multilevel model for the mortality
***### change the negative values to be positive values
	gen whz06new=0 if whz06>0
	replace whz06new=(whz06*-1) if whz06<=0
	gen rainX=rain_mm/50
	gen timeTR=time_to_readmission /30
	***when fitting the model for attributable fraction, we use power in the whz to get the best fit
	xtset fk_person count_adm
	xtdescribe
	****final model
	*gen whz_all=whz06
	*replace whz_all =_zbfa if whz_all==.
	
	***
	gen age_yrs_int = int(ageyrs)
	
	xtmelogit noutcome whz_all i.nsex  i.severe_disease   ///
	total_admission admdays timeTR age_yrs_int EVI_VALUE rainX || sublocs: ,    nolog 
	estat ic
	est store M7
	e
	
		
	*xtmelogit noutcome  i.nsex whz06new EVI_VALUE severe_disease || fk_person: ,  covariance(unstru) nolog
	*estat ic
	*est store M7b
	xtmelogit noutcome whz_all i.nsex malnutKid EVI_VALUE i.severe_disease ageyrs ///
	no_admissions || sublocs: ,  covariance(unstru) nolog or
	estat ic
	est store M8A
	
	xtmelogit noutcome  i.nsex malnutKid EVI_VALUE i.severe_disease ageyrs ///
	admdays no_admissions || sublocs: ,  covariance(unstru) nolog or
	estat ic
	est store M8
	
	
	*xtmelogit noutcome  i.nsex malnutKid EVI_VALUE i.severe_disease ageyrs ///
	*admdays no_admissions || fk_person: ,  covariance(unstru) nolog
	*estat ic
	*est store M8b
	
	
	
	xtmelogit noutcome  i.nsex malnutKid EVI_VALUE i.severe_disease ageyrs ///
	admdays no_admissions || sublocs:fk_person ,  covariance(unstru) nolog
	estat ic
	est store M8c
	
	
	*xtlogit noutcome  i.nsex malnutKid EVI_VALUE severe_disease ageyrs ///
	*admdays no_admissions  ,    nolog
	*estat ic
	*est store M9
	
	
	
	****bayesian model
	 bayesmh cumulitive_count i.gender i.severe_disease EVI_VALUE nagem, ///
 noconstant likelihood(poisson, exposure(cumulitive_time)) ///
 prior({cumulitive_count:i.gender}, normal(0,1)) ///
 prior({cumulitive_count:i.severe_disease}, normal(0,1)) ///
 prior({cumulitive_count:EVI_VALUE}, normal(0,1)) ///
 prior({cumulitive_count:nagem}, normal(0,1)) ///
 prior({cumulitive_count:_cons}, normal(0,1)) ///
 burnin(3000) mcmcsize(5000) dots 
 
 
 ****identity
   bayesmh cumulitive_count i.gender i.severe_disease EVI_VALUE nagem i.Adj_ID#i.fk_person, saving(poispostsl.dta, replace) ///
 noconstant likelihood(poisson, exposure(cumulitive_time)) ///
	prior({cumulitive_count:i.Adj_ID#i.fk_person}, normal({cumulitive_count:i.fk_person},{var_fk_person})) ///
	prior({var_0}, igamma(0.001, 0.001)) ///
	prior({var_Adj_ID}, igamma(0.001, 0.001)) ///
	prior({var_fk_person}, igamma(0.001, 0.001)) ///
	prior({cumulitive_count:i.gender}, normal(0,1)) ///
	prior({cumulitive_count:i.severe_disease}, normal(0,1)) ///
	prior({cumulitive_count:EVI_VALUE}, normal(0,1)) ///
	prior({cumulitive_count:nagem}, normal(0,1)) ///
	prior({cumulitive_count:_cons}, normal(0,1)) ///
	block({var_0}, gibbs) ///
	block({var_Adj_ID}, gibbs) ///
	block({var_fk_person}, gibbs) ///
	block({cumulitive_count:i.severe_disease}, gibbs) ///
	block({cumulitive_count:i.Adj_ID#i.fk_person}, gibbs) ///
	block({cumulitive_count:i.gender}, gibbs) ///
	block({cumulitive_count:_cons}, gibbs) ///
	block({cumulitive_count:nagem}, gibbs) ///
	block({cumulitive_count:EVI_VALUE}, gibbs) ///
 mcmcsize(5000) dots notable  
 
 **Unstructured 
  bayesmh cumulitive_count i.gender i.severe_disease EVI_VALUE nagem i.Adj_ID#i.fk_person, ///
 noconstant likelihood(poisson, exposure(cumulitive_time)) ///
 prior({cumulitive_count:i.Adj_ID }, ///
	mvnormal(2, {cumulitive_count:_cons}, {cumulitive_count:i.fk_person}, {Sigma,m})) ///
 prior({cumulitive_count:i.fk_person _cons}, normal(0, 1e2)) ///
 prior({var_0}, igamma(0.001,0.001)) ///
 prior({Sigma,m}, iwishart(2,3,I(2))) ///
 prior({cumulitive_count:i.gender}, normal(0,1)) ///
 prior({cumulitive_count:i.severe_disease}, normal(0,1)) ///
 prior({cumulitive_count:EVI_VALUE}, normal(0,1)) ///
 prior({cumulitive_count:nagem}, normal(0,1)) ///
 prior({cumulitive_count:_cons}, normal(0,1)) ///
 block({var_0}, gibbs) block({Sigma,m}, gibbs) ///
 block({cumulitive_count:i.gender}) block({cumulitive_count:i.severe_disease}) ///
 block({cumulitive_count:i.Adj_ID}, reffects) ///
 block({cumulitive_count:i.Adj_ID#i.fk_person}, reffects) ///
 noshow({cumulitive_count:i.Adj_ID#i.fk_person}) ///
 burnin(3000) mcmcsize(5000) dots notable nomodelsummary

 
 ****musenge bayesian
 bayesmh cumulitive_count i.gender i.severe_disease evi_value, reffects(spatial) initrandom ///
noconstant likelihood(poisson , exposure(age)) ///
prior({cumulitive_count:i.gender}, normal(0,1)) ///
prior({cumulitive_count:i.severe_disease}, normal(0,1)) ///
prior({cumulitive_count:evi_value}, normal(0,1)) ///
prior({cumulitive_count:i.spatial}, normal(3,{var_spatial})) ///
prior({cumulitive_count:_cons}, uniform(0,7)) ///
prior({var_spatial}, igamma(0.1, 0.1)) ///
prior({cumulitive_time}, normal(150,100000)) ///
block({var_spatial}, gibbs)



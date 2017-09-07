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
	local name = "logs\"+"morbidity_count"+"`tod'"+"~"+"`time'"
	log using `name',text replace	
//load haematology data from 1998-2008	from haematology folder
	***************set seed for the same random select any time******
	set seed 1221223

	****insheet data from stata_analysis excluding children with less than 1 admission
	use "data/morbidity.dta", clear
	*drop if count_adm>11
	outsheet using "data/morbidity.csv" , replace c
	
	******************save for AR check*****************
	preserve
	collapse  (count) totalMalnut=malnutKid   (sum) mal_con_T2=malnut2 mal_case_T=malnutKid severeD_T=severe_disease2 ///
	(mean) meanAge=nagem  meanEVI=EVI_VALUE , by(count_adm )
	gen  mal_con_T=totalMalnut-mal_case_T
	save "data/colapse_morbidity_count_adm", replace
	restore
	**********

	
///table
*****morbiditty objective	
	gen count_adm_cat=1 if count_adm==	1
	replace count_adm_cat=2 if count_adm==2
	replace count_adm_cat=3 if count_adm==3
	replace count_adm_cat=4 if count_adm==4
	replace count_adm_cat=5 if count_adm==5
	replace count_adm_cat=6 if count_adm>5 & count_adm!=.

	*** for under 5 analysis
	//drop if nagem>60.99
**summary continous
	xtset fk_person count_adm
	xtdescribe
	
	xtsum 	nagem admdays  no_admissions total_admission nweight  nheight nmuac   haz06 waz06 whz06 ///
	     EVI_VALUE rain_mm
	
	bys noutcome: xtsum 	nagem admdays  no_admissions total_admission nweight  nheight nmuac   haz06 waz06 whz06 ///
	     EVI_VALUE rain_mm  
	
*** summary categorical
foreach var of varlist nsex malnutKid severe_disease      ///
	severe_anaemia    hypoglycaemia         ///
	malaria    ndiarrhoea  menegitis lrti  ///
	gastroenteritis ntransfused  bc_pos csf_pos ///
	{
	xttab `var' 
}
	
foreach var of varlist nsex malnutKid severe_disease      ///
	severe_anaemia    hypoglycaemia         ///
	malaria    ndiarrhoea  menegitis lrti  ///
	gastroenteritis ntransfused  bc_pos csf_pos ///
	{
	xttab `var' if noutcome==0
}	
	
foreach var of varlist nsex malnutKid severe_disease      ///
	severe_anaemia    hypoglycaemia         ///
	malaria    ndiarrhoea  menegitis lrti  ///
	gastroenteritis ntransfused  bc_pos csf_pos ///
	{
	xttab `var' if noutcome==1
}		
	
	
foreach var of varlist nsex malnutKid severe_disease      ///
	severe_anaemia    hypoglycaemia         ///
	malaria    ndiarrhoea  menegitis lrti  ///
	gastroenteritis ntransfused  bc_pos csf_pos ///
	{
	bys noutcome: xttab `var' 
}	

****ancova models

anova cumulitive_count fk_person  count_adm, repeated(count_adm)

set mats

*anova cumulitive_count   nsex count_adm, repeated(count_adm) 
*anova cumulitive_count fk_person  noutcome count_adm, repeated(count_adm)
*anova cumulitive_count   nsex|fk_person count_adm, repeated(count_adm) 

	***testing for the differences using ttesti
	ttesti 3114 35.70618   32.35065   199 45.16001   38.68826
	ttesti 3113  5.149783   6.526309    189   6.513228   9.156795  
	ttesti 3114  3.500197   2.581411    199   2.959799    1.75188
	ttesti 3114  3.231859    2.28648     199   2.673367   1.370263
	ttesti 3108 10.96803   5.134796      195 10.53395   6.381912
	ttesti 3088 85.12218   18.97605     176   87.05341   22.44731
	ttesti 3112 13.87722   1.981982     192      12.35625   2.848948
	ttesti 2785 -1.651671   1.501167    127 -3.04622   1.606605
	ttesti 2776  -1.83153   1.480707    115 		-4.013391   1.693198
	ttesti 2742 -1.252816   1.506778      106    -3.388019   1.925638
	ttesti 3114 .3903826   .0821285       199    .3984442   .0828263
	ttesti 3068  31.77821   44.01779        190    37.26575   53.78804
	
	
	****chi
	
	********************88
	tabi 1751 199 \1363 82   ,chi2 exact
	tabi 2570 79 \987 116 ,chi2 exact
	tabi 1999  93 \2297  86 \308 10 \15 6 ,chi2 exact
	tabi 2984  160 \ 371 20 ,chi2 exact
	tabi 2800 131 \50 4 ,chi2 exact
	tabi 2634  155 \1022  26 ,chi2 exact
	tabi 2843 133 \1073 56 ,chi2 exact
	tabi 3062 178 \ 52  9 ,chi2 exact
	tabi 2698 152 \1351  35 ,chi2 exact
	tabi 2913 164 \775  23  ,chi2 exact
	tabi 3021 155 \367  31 ,chi2 exact
	tabi 3021 140 \220  41 ,chi2 exact
	tabi 1026 27 \11 5 ,chi2 exact

****table 1 output 
	tabout  nsex malnutKid severe_disease      ///
	severe_anaemia    hypoglycaemia         ///
	malaria    ndiarrhoea  menegitis lrti  ///
	gastroenteritis ntransfused  bc_pos csf_pos ///
	 count_adm_cat   using analysis/outputs/table1_morbi.txt , c(freq col) f(0c 1p) stats(chi2) rep
		
///numeric variables morbidit
	tabstat  nagem admdays  no_admissions total_admission nweight  nheight nmuac   haz06 waz06 whz06 ///
	     EVI_VALUE rain_mm , by(count_adm_cat) ///
	stat(mean sd p50 min max) nototal long col(stat)	
	
	
	gen rainX=rain_mm/50
	

	
	****
*#################################################################################****
 ****modelling the morbidity variable
 ****poisson regression model for the count of cases over time
	 mepoisson cumulitive_count  || sublocs: , irr
	 
	 ****declare the data to be a panel data
	**declaring fk_person to be the grouping variable
	*xtset fk_person
	****use this to specify 
	gen mnthyr= ym(  yr, mnths)
	format mnthyr  %tm
	xtset fk_person  count_adm
  
	*xtset pan_id
	*xtnbreg cumulitive_count total_admission,nolog exposure(nagem_int) irr 

	*menbreg cumulitive_count , nolog exposure(nagem_int) || sublocs: || fk_person:
*****model with just the random effects	
	 menbreg cumulitive_count total_admission   ,irr nolog  exposure(nagem_int)  || sublocs: 
	 estat ic


	 *menbreg cumulitive_count   ,  exposure(nagem_int)   nolog || sublocs: ||fk_person:
	 est store M1
	 estat ic
	


* mepoisson cumulitive_count EVI_VALUE   nagem i.nsex   ///
	 *|| sublocation:    , irr cov(unstructured)
*****model without severe disease
	* menbreg cumulitive_count EVI_VALUE   nagem i.nsex   total_admission,   exposure(nagem)   nolog || fk_person:
	* est store M2
	* estat ic
	 
	 menbreg cumulitive_count   i.nsex  total_admission , irr  exposure(nagem)   nolog || sublocs:
	 est store M2b
	 estat ic
	 
	* menbreg cumulitive_count EVI_VALUE    i.nsex   total_admission,   exposure(nagem)  nolog || sublocs:fk_person
	* est store M2c
	* estat ic
	* menbreg cumulitive_count EVI_VALUE    i.nsex i.severe_disease total_admission   , irr  exposure(nagem)  nolog || sublocs:
	 
****model with severe disease
	 menbreg cumulitive_count EVI_VALUE    i.nsex i.severe_disease total_admission   , irr  exposure(nagem)  nolog || sublocs:
	 est store M3
	 estat ic
	 
	 menbreg cumulitive_count EVI_VALUE rainX    i.nsex i.severe_disease total_admission admdays  ,  irr exposure(nagem)  nolog || sublocs:
	 est store M3b
	 estat ic
	 
	 menbreg cumulitive_count EVI_VALUE  rainX  i.nsex i.severe_disease total_admission admdays  ,  irr exposure(nagem)  nolog || sublocs:

	 
*****final model*****
	 menbreg cumulitive_count EVI_VALUE  rainX  i.nsex   total_admission admdays  nweight  ,  exposure(nagem)  nolog || sublocs:
	 est store M3c
	 estat ic
	 
	 menbreg cumulitive_count EVI_VALUE  rainX  i.nsex  i.severe_disease total_admission admdays  nweight  ,  exposure(nagem_int )  nolog || sublocs:
	 est store M3d
	 estat ic
	 
	 
		*estat icc

	 
	*menbreg cumulitive_count EVI_VALUE   nagem i.nsex i.severe_disease total_admission  ,   exposure(nagem_int)  nolog || fk_person:
	 *est store M3b
	 *estat ic
	 
	* menbreg cumulitive_count EVI_VALUE   nagem i.nsex i.severe_disease  total_admission,   exposure(nagem_int) nolog || sublocs:fk_person
	* est store M3c
	* estat icS
	 	 
	*xtnbreg cumulitive_count  EVI_VALUE   nagem i.nsex i.severe_disease total_admission, exposure(nagem_int)  normal nolog
	*est store M4
	*estat ic

	
*** bayesian model for the 	
	
	
	
*****bayesSMH
****musenge bayesian
	bayesmh cumulitive_count EVI_VALUE  rainX  i.nsex i.severe_disease total_admission admdays  nweight count_adm, reffects(sublocs)  ///
	noconstant likelihood(poisson , exposure(nagem)) ///
	prior({cumulitive_count:EVI_VALUE}, normal(0.38,0.01)) ///
	prior({cumulitive_count:rainX}, normal( 0.63, 0.88)) ///
	prior({cumulitive_count:i.nsex},  beta(0.5,0.5)) ///
	prior({cumulitive_count:i.severe_disease}, normal(0,1)) ///
	prior({cumulitive_count:total_admission}, normal( 3.12,4.08)) ///
	prior({cumulitive_count:admdays}, normal( 4.63,5.96)) ///
	prior({cumulitive_count:nweight}, normal( 10.65,4.92)) ///
	prior({cumulitive_count:count_adm}, normal(0, 1	)) ///
	prior({cumulitive_count:i.sublocs}, normal(3,{var_sublocs})) ///
	prior({cumulitive_count:_cons}, uniform(0,7)) ///
	prior({var_sublocs}, igamma(0.1, 0.1)) ///
	block({var_sublocs}, gibbs) ///
	burnin(3000) mcmcsize(10000) dots  thinning(10) ///
	saving("data/bayesmh_resul.dta" , replace) 
	
	est store bModel
	bayesstats ic 

	*prior({cumulitive_time}, normal(150,100000)) ///
	bayesstats summary _all 
	bayesstats ic 
	bayestest model 
	bayesgraph diagnostics _all 

	
	****musenge bayesian
	bayesmh cumulitive_count EVI_VALUE  rainX  i.nsex i.severe_disease total_admission admdays  nweight count_adm, reffects(sublocs)  ///
	noconstant likelihood(poisson , exposure(nagem)) ///
	prior({cumulitive_count:EVI_VALUE}, normal(0.38,0.01)) ///
	prior({cumulitive_count:rainX}, normal( 0.63, 0.88)) ///
	prior({cumulitive_count:i.nsex},  beta(0.5,0.5)) ///
	prior({cumulitive_count:i.severe_disease}, normal(0,1)) ///
	prior({cumulitive_count:total_admission}, normal( 3.12,4.08)) ///
	prior({cumulitive_count:admdays}, normal( 4.63,5.96)) ///
	prior({cumulitive_count:nweight}, normal( 10.65,4.92)) ///
	prior({cumulitive_count:count_adm}, normal(count_adm[_n-1], 1	)) ///
	prior({cumulitive_count:i.sublocs}, normal(3,{var_sublocs})) ///
	prior({cumulitive_count:_cons}, uniform(0,7)) ///
	prior({var_sublocs}, igamma(0.1, 0.1)) ///
	block({var_sublocs}, gibbs) ///
	burnin(3000) mcmcsize(10000) dots  thinning(10) ///
	saving("data/bayesmh_resul.dta" , replace) 

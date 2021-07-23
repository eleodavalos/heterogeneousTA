*******Databases*********
global finaldatabase  "C:\Users\zurit\Dropbox\MTEAgricola\Submissions_MTEAgricola\JAAE\Final_Version_JAAE\4github_heterogeneousTA"
global results "D:\Dropbox\MTEAgricola\Estimaciones\Results"

cd "$finaldatabase"

set more off

use "Base_estimaciones_Final_JAAE",clear

gen log_rend=log_p1-log_areasem1

***************Globals of control variables and instruments************

*Controles 
global X  tm cons_agro tot_tier P_S11P138
global P P_S15P168 sex_j edu_j edu edad_j P_S15P169 prod P_S15P165
global M area altura dis_capital dis_mayoris dis_bogota dis_otrosmer region_andina region_caribe region_pasifica 

**Instrumentos
global Z des log_cos_mean

**********Robustness estimations**********************

****Different instruments

**Cost instrument Bartik Median
ivreghdfe log_rend $X $P (trat=des log_cos_meadian) ,abs(COD_VEREDA) cl(COD_VEREDA)
estimates store iv1_fe
estadd scalar HJ1=e(jp)

**Cost instrument Bartik Weight
ivreghdfe log_rend $X $P (trat=des log_cos_weight) ,abs(COD_VEREDA) cl(COD_VEREDA)
estimates store iv2_fe
estadd scalar HJ2=e(jp)

**Cost instrument
ivreghdfe log_rend $X $P (trat=des log_cos_mean) ,abs(COD_VEREDA) cl(COD_VEREDA)
estimates store iv3_fe
estadd scalar HJ3=e(jp)

global labels trat AsisTéncinca tm TenenciaMaquinaria cons_agro ConstruccionesAgro tot_tier ÁreaUPA P_S11P138 TrabajoPermanente fenomenos DesastresNaturalez tenencia TenenciaPrivada raza Predominanciaétnica ///
P_S15P168 Hombresxhogar sex_j Hombresjefesdehogar edu_j Educaciónjefedehogar edu Educaciónhogar edad_j Edadjefedehogar P_S15P169 Edadhogar prod Fuerzadetrabajo P_S15P165 Tamañohogar ///
area Área altura Altura ind_rur Índiceruralidad dis_capital DistanciaCapital dis_mayoris DistanciaMayorista dis_bogota DistanciaBogota dis_otrosmer Distanciaotromercado region_andina RegiónAndina region_caribe RegiónCaribe region_pasifica RegiónPasifica _cons Constante ///
c_trat AsisTéncinca c_log_areasem1 ÁreaSembrada c_tm TenenciaMaquinaria c_cons_agro ConstruccionesAgro c_tot_tier ÁreaUPA c_P_S11P138 TrabajoPermanente c_fenomenos DesastresNaturalez c_tenencia TenenciaPrivada c_raza Predominanciaétnica ///
c_P_S15P168 Hombresxhogar c_sex_j Hombresjefesdehogar c_edu_j Educaciónjefedehogar c_edu Educaciónhogar c_edad_j Edadjefedehogar c_P_S15P169 Edadhogar c_prod Fuerzadetrabajo c_P_S15P165 Tamañohogar ///
c_area Área c_altura Altura c_ind_rur Índiceruralidad c_dis_capital DistanciaCapital c_dis_mayoris DistanciaMayorista c_dis_bogota DistanciaBogota c_dis_otrosmer Distanciaotromercado c_region_andina RegiónAndina c_region_caribe RegiónCaribe c_region_pasifica RegiónPasifica _cons Constante

esttab iv1_fe iv2_fe iv3_fe using "$results\Robustness_IV.csv",replace cells(b(star fmt(%9.3f)) se(par)) ///
           stats(N HJ1 HJ2 HJ3 r2_a , fmt(%9.0g %9.3f %9.3f %9.3f))      ///
            legend label mtitle("ETC+Bartik median" "ETC+Bartik weight" "ETC+Cost") varlabels($labels) star(* 0.1 ** 0.05 *** 0.01)


**Excluding main coffe producers

replace area_sem_cafe=0 if area_sem_cafe ==.
gen per_cafe=area_sem_cafe / AREA_SEMBRADA
cap drop id_cafe
gen id_cafe=(per_cafe>=.5)

**IV regressions
ivreghdfe log_rend $X $P (trat=des log_cos) if id_cafe==0,absorb(COD_VEREDA) cl(COD_VEREDA)
estimates store iv1_main

esttab iv1_main using "$results\Robustness_Coffe.csv",replace cells(b(star fmt(%9.3f)) se(par)) ///
           stats(N , fmt(%9.0g)) star(* 0.1 ** 0.05 *** 0.01)

**Average of variables at Vereda level		   
global V $X $P $Z
foreach x of global V{
bysort COD_VEREDA: egen v_`x'=mean(`x')
}

global V  v_tm v_cons_agro v_tot_tier v_P_S11P138 ///
  v_P_S15P168 v_sex_j v_edu_j v_edu v_edad_j v_P_S15P169 v_prod v_P_S15P165 ///
  v_des v_log_cos

cd "$results"
		   
**MTE estimation
movestay log_rend $X $P $V if id_cafe==0, select(trat= $Z $X $P $V) cl(COD_VEREDA)
cap drop sample 
gen sample=e(sample)

do "$finaldatabase\estimacion_parametros"

putexcel set  "parameters_main_cofe.xlsx", replace
putexcel B2=matrix(PAR), rownames

**Bootstrap

matrix Boot=J(100,5,.)

local Rep=100
forvalues r=1/`Rep'{
preserve 
count
keep if sample==1
count
bsample _N 
count

movestay log_rend $X $P $V if id_cafe==0, select(trat= $Z $X $P $V) cl(COD_VEREDA)
cap drop sample 
gen sample=e(sample)


cap drop Pz-ate_boot
cap drop MPRTE1 MPRTE2 ATT ATUT 
do "$finaldatabase\estimacion_parametros"

matrix Boot[`r',1]=PAR'
restore
}

cap drop Pz-ate_boot
cap drop MPRTE1 MPRTE2 ATT ATUT 
mat colnames Boot = MPRTE1 MPRTE2 ATT ATUT ATE 

svmat2 Boot,name(col)

matrix sd_main_cofe=J(5,1,.)
local i=1
foreach x in MPRTE1 MPRTE2 ATT ATUT ATE  {
sum `x'
matrix sd_main_cofe[`i',1]=r(sd)
local i=`i'+1
}

putexcel set  "sd_main_cofe.xlsx", replace
putexcel B2=matrix(sd_main_cofe), rownames

drop MPRTE1 MPRTE2 ATT ATUT ATE


**************Estimacion sin unidades con cultivos de cafe
qui ivreghdfe log_rend $X $P (trat=$Z) if area_cafe==0,absorb(COD_VEREDA) cl(COD_VEREDA)

movestay log_rend $X $P $V if e(sample), select(trat= $Z $X $P $V) cl(COD_VEREDA)
drop sample 
gen sample=e(sample)

cap drop MPRTE1 MPRTE2 ATT ATUT 
do "$finaldatabase\estimacion_parametros"

putexcel set  "parameters_all_cofe.xlsx", replace
putexcel B2=matrix(PAR), rownames

**Bootstrap

matrix Boot=J(100,5,.)

local Rep=100
forvalues r=1/`Rep'{
preserve 
count
keep if sample==1
count
bsample _N 
count

movestay log_rend $X $P $V if area_cafe==0, select(trat= $Z $X $P $V) cl(COD_VEREDA)
cap drop sample 
gen sample=e(sample)


cap drop Pz-ate_boot
cap drop MPRTE1 MPRTE2 ATT ATUT 
do "$finaldatabase\estimacion_parametros"

matrix Boot[`r',1]=PAR'
restore
}

cap drop Pz-ate_boot
cap drop MPRTE1 MPRTE2 ATT ATUT 
mat colnames Boot = MPRTE1 MPRTE2 ATT ATUT ATE 

svmat2 Boot,name(col)

matrix sd_all_cofe=J(5,1,.)
local i=1
foreach x in MPRTE1 MPRTE2 ATT ATUT ATE  {
sum `x'
matrix sd_all_cofe[`i',1]=r(sd)
local i=`i'+1
}

putexcel set  "sd_all_cofe.xlsx", replace
putexcel B2=matrix(sd_all_cofe), rownames


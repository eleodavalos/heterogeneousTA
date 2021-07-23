*******Databases*********
global finaldatabase "~\Dropbox\MTEAgricola\Estimaciones"
global results "~\Dropbox\MTEAgricola\Submissions_MTEAgricola\JAAE\R&R_submission2\Resultados"

cd "$finaldatabase"

set more off

use "Base_estimaciones_Final_JAAE",clear

gen log_rend=log_p1-log_areasem1

***************Globals of control variables and instruments************

*Controles 
global X tm cons_agro tot_tier P_S11P138 
global P P_S15P168 sex_j edu_j edu edad_j P_S15P169 prod P_S15P165

**Instrumentos
global Z des log_cos


***********MTE Parametric regression*************

**Average of variables at Vereda level
global V $X $P $Z
foreach x of global V{
bysort COD_VEREDA: egen v_`x'=mean(`x')
}

global V  v_tm v_cons_agro v_tot_tier v_P_S11P138 ///
  v_P_S15P168 v_sex_j v_edu_j v_edu v_edad_j v_P_S15P169 v_prod v_P_S15P165 ///
  v_des v_log_cos

***************Tamaño de la tierra

qui movestay log_rend $X $P $V , select(trat= $Z $X $P $V) cl(COD_VEREDA)

gen sample=e(sample)

gen area_upa=exp(tot_tier)

_pctile area_upa if sample==1, p(25 50 75)

gen per=1 if sample==1 & area_upa<=`r(r1)'
replace per=2 if sample==1 & area_upa>`r(r1)' & area_upa<=`r(r2)'
replace per=3 if sample==1 & area_upa>`r(r2)' & area_upa<=`r(r3)'
replace per=4 if sample==1 & area_upa>`r(r3)'

cd "$results"

forvalues i=2/4{
cap movestay log_rend $X $P $V if per== `i', select(trat= $Z $X $P $V) cl(COD_VEREDA)
cap drop sample 
gen sample=e(sample)

cap drop MPRTE1 MPRTE2 ATT ATUT 
do "$finaldatabase\estimacion_parametros"
putexcel set  "parameters_size`i'.xlsx", replace
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

cap movestay log_rend $X $P $V, select(trat= $Z $X $P $V) cl(COD_VEREDA)
drop sample 
gen sample=e(sample)


cap drop Pz-ate_boot
cap drop MPRTE1 MPRTE2 ATT ATUT 
do "$finaldatabase\estimacion_parametros"

matrix Boot[`r',1]=PAR'
restore
}

cap drop Pz-ate_boot

mat colnames Boot = MPRTE1 MPRTE2 ATT ATUT ATE 

cap drop MPRTE1 MPRTE2 ATT ATUT 

svmat2 Boot,name(col)

matrix sd_main_cofe=J(5,1,.)
local f=1
foreach x in MPRTE1 MPRTE2 ATT ATUT ATE  {
sum `x' if `x'<=2,d
matrix sd_main_cofe[`f',1]=r(sd)
local f=`f'+1
}

putexcel set  "sd_parameters_size`i'.xlsx", replace
putexcel B2=matrix(sd_main_cofe), rownames

drop MPRTE1 MPRTE2 ATT ATUT ATE

}


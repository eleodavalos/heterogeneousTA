*******Databases*********
global finaldatabase "C:\Users\zurit\Dropbox\MTEAgricola\Submissions_MTEAgricola\JAAE\Final_Version_JAAE\4github\GithubTorresDavalosMorales2021"
global results "D:\Dropbox\MTEAgricola\Estimaciones\Results "

cd "$finaldatabase"
*ssc install reghdfe
*ssc install center
*ssc install avar

set more off

use "Base_estimaciones_Final_JAAE",clear

gen log_rend=log_p1-log_areasem1

sum log_rend log_p1 log_areasem1 if trat!=.

***************Globals of control variables and instruments************

*Controles 
global X  tm cons_agro tot_tier P_S11P138
global P P_S15P168 sex_j edu_j edu edad_j P_S15P169 prod P_S15P165
global M area altura dis_capital dis_mayoris dis_bogota dis_otrosmer region_andina region_caribe region_pasifica 

**Instrumentos
global Z des log_cos_mean

reghdfe trat $X $P des if log_cos!=. & log_rend!=.,absorb(COD_VEREDA) cl(COD_VEREDA)

matrix  T=J(35,7,0)
local j=1
foreach var in $X $P des log_cos log_rend {
ttest `var' if e(sample) , by(trat) unequal
matrix T[`j',1]=r(mu_1)
matrix T[`j',2]=r(sd_1)
matrix T[`j',3]=r(N_1)
matrix T[`j',4]=r(mu_2)
matrix T[`j',5]=r(sd_2)
matrix T[`j',6]=r(N_2)
matrix T[`j',7]=r(p)
local j=`j'+1
}

matrix list T

sum area
ttest area      if e(sample) , by(trat) unequal
ttest tot_tier  if e(sample) , by(trat) unequal

gen Etot_tier=exp(tot_tier)
ttest Etot_tier if e(sample) , by(trat) unequal



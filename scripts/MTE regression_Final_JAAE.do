*******Databases*********
global finaldatabase "C:\Users\zurit\Dropbox\MTEAgricola\Submissions_MTEAgricola\JAAE\Final_Version_JAAE\4github\GithubTorresDavalosMorales2021"
global results "D:\Dropbox\MTEAgricola\Estimaciones\Results"

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


movestay log_rend $X $P $V , select(trat= $Z $X $P $V ) cl(COD_VEREDA)
estimate store margte
cap drop sample
gen sample=e(sample)
predict ps, psel 
matrix coef=e(b)
matrix vari=e(V)
matrix b_0=coef[1,"log_rend0:"]
matrix b_1=coef[1,"log_rend1:"]
matrix dif=b_1-b_0 
matrix rho=[-(exp(_b[lns1:_cons]) * tanh(_b[r1:_cons])),-(exp(_b[lns0:_cons]) * tanh(_b[r0:_cons])),-(exp(_b[lns1:_cons]) * tanh(_b[r1:_cons]))+(exp(_b[lns0:_cons]) * tanh(_b[r0:_cons]))]

estadd scalar rho1=rho[1,1]
estadd scalar rho2=rho[1,2]
estadd scalar rho3=rho[1,3]

matrix se_Rho=J(1,3,.)
nlcom -(exp(_b[lns1:_cons]) * tanh(_b[r1:_cons]))
matrix se_Rho[1,1]=r(V)
nlcom -(exp(_b[lns0:_cons]) * tanh(_b[r0:_cons]))
matrix se_Rho[1,2]=r(V)
nlcom -(exp(_b[lns1:_cons]) * tanh(_b[r1:_cons]))+(exp(_b[lns0:_cons]) * tanh(_b[r0:_cons]))
matrix se_Rho[1,3]=r(V)

estadd scalar se_rho1=se_Rho[1,1]^0.5
estadd scalar se_rho2=se_Rho[1,2]^0.5
estadd scalar se_rho3=se_Rho[1,3]^0.5

global X_bar $X $P $V

local i=0
foreach x of global X_bar{
local i=`i'+1
qui sum `x' if sample==1
local mte `mte' +(_b[log_rend1:`x']-_b[log_rend0:`x'])*`r(mean)'
}	

set matsize 1000
matrix var_mte=J(999,1,.)
matrix MTE=J(999,1,.)

forvalues i=1/999{
matrix mte_`i'= _b[log_rend1:_cons]-_b[log_rend0:_cons] `mte' -(-(exp(_b[lns1:_cons]) * tanh(_b[r1:_cons]))+(exp(_b[lns0:_cons]) * tanh(_b[r0:_cons])))*invnormal(`i'/1000)
matrix MTE[`i',1]=mte_`i'

nlcom _b[log_rend1:_cons]-_b[log_rend0:_cons] `mte' -(-(exp(_b[lns1:_cons]) * tanh(_b[r1:_cons]))+(exp(_b[lns0:_cons]) * tanh(_b[r0:_cons])))*invnormal(`i'/1000)
matrix var_mte[`i',1]=r(V)

}

svmat2 double MTE,name(mte)

svmat2 double var_mte,name(var_mte)	
gen evalgrid=_n/1000 if mte!=.
gen ci_ll=mte-var_mte^(0.5)*invnormal(0.975)
gen ci_ul=mte+var_mte^(0.5)*invnormal(0.975)

twoway (rarea ci_ul ci_ll evalgrid, fcolor(gs12) lcolor(gs11)) (line mte evalgrid, ///
 lcolor(black)),legend(order(1 "95% CI" 2 "MTE") position(6) cols(2)) xtitle("U_D") scheme(plottig) ///
 plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white))
graph export "$finaldatabase/MTE_estimation.eps",replace 
graph export "$finaldatabase/MTE_estimation.pdf",replace 
graph export "$finaldatabase/MTE_estimation.png",replace 


egen    gridv = fill(1/2) 
replace gridv = (gridv)/1000
replace gridv = . if gridv>0.9999

******************************ATE,ATT, ATUT, MPRTE computation***************

rename mte mte_p
rename ps Pz

capture drop h_mprte1-per_h_mprte2_wox MPRTE1 MPRTE2
global mte  "mte_p"


gen h_mprte1=.
gen h_mprte2=.

qui kdensity Pz, n(999) k(gauss) nograph
dis r(bwidth)

gen fp_v = .
local h = r(bwidth) /*optimal bandwith*/
dis `h'
qui levelsof gridv, local(xs)
local i 1
qui {
foreach l of local xs {
	dis `i'
	local j = `l'*1000
	gen double w = normalden((Pz-`l')/`h')/`h' 
	sum w
	replace fp_v = r(mean) in `i'
	drop w
	loc i=`i'+1
}
}

*
replace fp_v=fp_v/999
gen v_fp_v = fp_v*gridv

gen per_h_mprte1_wox = fp_v


egen temp1 = total(v_fp_v) in 1/999
gen per_h_mprte2_wox = (v_fp_v/temp1)
drop temp1

gen per_h_att_wox=.
egen temp1= total(v_fp_v) in 1/999
forvalue i=1/999{
qui total fp_v in `i'/999
qui replace per_h_att_wox=_b[fp_v]/temp1 in `i'
}
drop temp1
egen temp =total(per_h_att_wox) in 1/999
replace per_h_att_wox=per_h_att_wox/temp
drop temp

gen per_h_atut_wox=.
egen temp1= total(v_fp_v) in 1/999
forvalue i=1/999{
qui total fp_v in 1/`i'
qui replace per_h_atut_wox=_b[fp_v]/temp1 in `i'
}
drop temp1
egen temp =total(per_h_atut_wox) in 1/999
replace per_h_atut_wox=per_h_atut_wox/temp
drop temp

egen MPRTE1 = total(${mte}*per_h_mprte1_wox) in 1/999               
sum MPRTE1
matrix prob=r(mean)
qui svmat2 double prob,name(mprte1_boot)
estadd scalar MPRTE1=r(mean)

egen MPRTE2 = total(${mte}*per_h_mprte2_wox) in 1/999           
sum MPRTE2
matrix prob=r(mean)
qui svmat2 double prob,name(mprte2_boot)
estadd scalar MPRTE2=r(mean)

egen ATT = total(${mte}*per_h_att_wox) in 1/999  
sum ATT
matrix prob=r(mean)
qui svmat2 double prob,name(att_boot)
estadd scalar ATT=r(mean)

egen ATUT = total(${mte}*per_h_atut_wox) in 1/999   
sum ATUT
matrix prob=r(mean)
qui svmat2 double prob,name(atut_boot)
estadd scalar ATUT=r(mean)

sum ${mte}
matrix prob=r(mean)
qui svmat2 double prob,name(ate_boot)
estadd scalar ATE=r(mean)

matrix PAR=J(5,1,.)
local i=1
foreach x in mprte1_boot mprte2_boot att_boot atut_boot ate_boot {
qui sum `x'
matrix PAR[`i',1]=r(mean)
local i=`i'+1
}

* summarize parameters

drop mte_*

global labels trat AsisTéncinca log_areasem1 ÁreaSembrada tm TenenciaMaquinaria cons_agro ConstruccionesAgro tot_tier ÁreaUPA P_S11P138 TrabajoPermanente fenomenos DesastresNaturalez tenencia TenenciaPrivada raza Predominanciaétnica ///
P_S15P168 Hombresxhogar sex_j Hombresjefesdehogar edu_j Educaciónjefedehogar edu Educaciónhogar edad_j Edadjefedehogar P_S15P169 Edadhogar prod Fuerzadetrabajo P_S15P165 Tamañohogar ///
area Área altura Altura ind_rur Índiceruralidad dis_capital DistanciaCapital dis_mayoris DistanciaMayorista dis_bogota DistanciaBogota dis_otrosmer Distanciaotromercado region_andina RegiónAndina region_caribe RegiónCaribe region_pasifica RegiónPasifica _cons Constante

cd "$results"
esttab margte using "margte_fe.csv",replace cells(b(star fmt(%9.3f)) se(par)) ///
           scalars(N rho1 se_rho1 rho2 se_rho2 rho3 se_rho3) sfmt(%9.0f %9.3f)  ///
            legend label mtitle("Normal" ) varlabels($labels) star(* 0.1 ** 0.05 *** 0.01)

putexcel set "parameters_mte.xlsx",replace
putexcel B2 = matrix(PAR), rownames

			
*********************************Standard errors parameters (Bootstrap procedure)************************************
		
*****Bootstrap ATE, ATT y ATUT********************
cap drop att* atut* mprte* mte_p Pz gridv
cap drop h_mprte1 h_mprte2 fp_v v_fp_v per_h_mprte1_wox per_h_mprte2_wox per_h_att_wox per_h_atut_wox MPRTE1 MPRTE2 ATT ATUT
local Rep=100
matrix Boot=J(`Rep',5,.)
forvalues r=1/`Rep'{
preserve 
count
keep if sample==1
count
bsample _N 
count

qui movestay log_rend $X $P $V, select(trat=$Z $X $P $V ) cl(COD_VEREDA)
drop sample 
gen sample=e(sample)
do "$finaldatabase\estimacion_parametros"

matrix Boot[`r',1]=PAR'

restore
}

cap drop Pz-ate_boot

mat colnames Boot = MPRTE1 MPRTE2 ATT ATUT ATE 

svmat2 Boot,name(col)

matrix sd_main_cofe=J(5,1,.)
local i=1
foreach x in MPRTE1 MPRTE2 ATT ATUT ATE  {
sum `x'
matrix sd_main_cofe[`i',1]=r(sd)
local i=`i'+1
}

cd "$results"
putexcel set  "sd_parameters.xlsx", replace
putexcel B2=matrix(sd_main_cofe), rownames


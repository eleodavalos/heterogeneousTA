cap drop ps
cap drop mte
predict ps if sample==1, psel 
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
cap nlcom -(exp(_b[lns1:_cons]) * tanh(_b[r1:_cons]))
matrix se_Rho[1,1]=r(V)
cap nlcom -(exp(_b[lns0:_cons]) * tanh(_b[r0:_cons]))
matrix se_Rho[1,2]=r(V)
cap nlcom -(exp(_b[lns1:_cons]) * tanh(_b[r1:_cons]))+(exp(_b[lns0:_cons]) * tanh(_b[r0:_cons]))
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
}

svmat2 double MTE,name(mte)

***MPRTE Estimation using Carneiro et al 2011 method
cap drop gridv
egen    gridv = fill(1/2) 
replace gridv = (gridv)/1000
replace gridv = . if gridv>0.9999

cap drop mte_p
rename mte mte_p
cap drop Pz  
rename ps Pz

capture drop h_mprte1-per_h_atut_wox
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
	*gen double w=max(0,`h'-abs(Pz-`l'))       /*K=Triangular*/
	gen double w = normalden((Pz-`l')/`h')/`h' /*K=Gauss*/
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

matrix PAR=J(5,1,.)

egen MPRTE1 = total(${mte}*per_h_mprte1_wox) in 1/999                  
sum MPRTE1
matrix PAR[1,1]=r(mean)

egen MPRTE2 = total(${mte}*per_h_mprte2_wox) in 1/999                  
sum MPRTE2
matrix PAR[2,1]=r(mean)

egen ATT = total(${mte}*per_h_att_wox) in 1/999                  
sum ATT
matrix PAR[3,1]=r(mean)

egen ATUT = total(${mte}*per_h_atut_wox) in 1/999                  
sum ATUT
matrix PAR[4,1]=r(mean)

sum mte
matrix PAR[5,1]=r(mean)

mat rownames PAR = MPRTE1 MPRTE2 ATT ATUT ATE 

drop mte_*

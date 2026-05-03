cd "C:\Users\olofd\PhD\Mertrics 2\PS\4"

clear all 

**# 1.1-7

	
	
frame create measures
frame change measures
set obs 11

gen theta = -1 + (_n-1)*0.2

replace theta = round(theta, 0.1)


foreach t in " " "V"{
	forvalues i = 1/2{
		gen `t'D_`i' = 0
		gen `t'W_`i' = 0
	}
	
}

gen tao_D = 0
gen tao_W = 0

gen FEcoef = 0

gen DW = 0

gen VDW = 0

frame change default

local c = 0
	
set seed 15
forvalues p = 1/11{
	
	quietly{
	local theta = -1 
	//+ (`p' - 1)*0.2
	
	set obs 1000

	gen j = floor((_n+9)/10)

	bys j: gen i = _n

	gen D = _n>50
	bys j: gen W = _n > 5
	
	gen D_W = D*W

	gen Y_00 = rnormal(0,1)

	gen Y_10 = rnormal(1,1)

	bys j: gen mu_t = rnormal(2,1) if _n ==1
	bys j: egen mu_j = max(mu_t)
	drop mu_t
	
	gen X_ij = rnormal(`theta', 0.2)
	
	
	gen Y_01 = rnormal(mu_j,1)
	
	gen Y_11 = Y_10 + Y_01 + X_ij
	
	gen Y_ij = Y_00 + D*(Y_10 - Y_00) + W*(Y_01 - Y_00) + D*W*(Y_11 - Y_01 - Y_10 + Y_00) 
	
	gen tao_D = (Y_10 - Y_00)
	gen tao_W = (Y_01 - Y_00)
	
	
	
	local i = 0
	foreach v in " " "D_W"{
		local i = `i' + 1
		
		reg Y_ij D W `v' , robust
		
		
		local D_`i'_l = r(table)[1,1]
		local VD_`i'_l = r(table)[2,1]
		
		local W_`i'_l = r(table)[1,2]
		local VW_`i'_l = r(table)[2,2]
		
		if "`v'" == "D_W"{
			local DW_l = r(table)[1,3]
			local VDW_l = r(table)[2,3]
		}
		
		
	}
	
	reg Y_ij W i.j
	local FEcoef_l = r(table)[1,1]
	
	
	summarize tao_D 
	
	local tao_D_l = r(mean)
	
	summarize tao_W
	local tao_W_l = r(mean)
	
	frame change measures
	
	foreach t in "" "V"{
	forvalues i = 1/2{
		
		replace `t'D_`i' = ``t'D_`i'_l' if `p' == _n
		replace `t'W_`i' = ``t'W_`i'_l' if `p' == _n
		
	}
	
	replace tao_D = `tao_D_l' if `p' == _n
	replace tao_W = `tao_W_l' if `p' == _n
	
	replace FEcoef = `FEcoef_l' if `p' == _n
	
	replace DW = `DW_l' if `p' == _n
	
	replace VDW = `VDW_l' if `p' == _n
	
}
	frame change default
	drop _all
	}
}
frame change measures


foreach t in "D" "W"{
	forvalues i = 1/2{
		gen upper = 1.96*V`t'_`i' + `t'_`i' 

		gen lower  = -1.96*V`t'_`i' + `t'_`i'
		
		graph twoway (rcap upper lower theta) (scatter `t'_`i' theta) (scatter tao_`t' theta), legend(label(1 "95% CI") label(2 "Estimate") label(3 "True estimand"))
		
		graph export "bilder\\`t'_`i'.png", replace
		
		drop upper lower
	}
	
}




**# 1. 8-9

clear all


frame change default 
set seed 15


set obs 1000

gen j = floor((_n+9)/10)

bys j: gen i = _n

gen D = _n>50
bys j: gen W = _n > 5
	
gen D_W = D*W

gen Y_00 = rnormal(0,1)

gen Y_10 = rnormal(1,1)

bys j: gen mu_t = rnormal(2,1) if _n ==1
bys j: egen mu_j = max(mu_t)
drop mu_t
	
gen X_ij = rnormal(0, 0.2)
	
	
gen Y_01 = rnormal(mu_j,1)
	
gen Y_11 = Y_10 + Y_01 + X_ij
	
gen Y_ij = Y_00 + D*(Y_10 - Y_00) + W*(Y_01 - Y_00) + D*W*(Y_11 - Y_01 - Y_10 + Y_00)  
	
gen tao_D = (Y_10 - Y_00)
gen tao_W = (Y_01 - Y_00)

gen FE = 0

reg Y_ij i.j 

forvalues j = 1/100{
	local fe = r(table)[1,`j']
	
	local cons = r(table)[1, 101]
	
	replace FE = `fe' + `cons' if j == `j'
	
}

gen first_mu = mu_j if _n == 1
egen temp_1 = max(first_mu)

gen implied_mu = temp_1 + FE
gen implied_mu_2 = 2 + FE

drop first_mu temp_1

bys j: gen abs_error = abs(mu_j - implied_mu) if j != 1 & _n == 1

summarize abs_error

egen error_sum = sum(abs_error)


preserve
collapse (first) FE implied_mu mu_j implied_mu_2, by(j)

twoway(scatter FE mu_j) (line mu_j mu_j), xtitle("True estimand") ytitle("Estimate") legend(off)

graph export "bilder\scatterFE.png", replace
restore

//Implied lägger till mu(1) men denna kan vara felestimerad från första början. Därmed lägga till mean mu, då blir det bra. 

//Empirical Bayes

egen Y_bar_j = mean(Y_ij), by(j)

egen mu_hat = mean(Y_bar_j)

gen eps_sq = (Y_ij - Y_bar_j)^2

egen sum_eps_sq = sum(eps_sq)

gen mean_eps_sq = sum_eps_sq/(9*100)

bys j: gen alpha_sq = (Y_bar_j - mu_hat)^2 if _n == 1

egen sum_alpha_sq = sum(alpha_sq)

gen mean_alpha_sq = sum_alpha_sq/99 - mean_eps_sq/10

gen lambda = (mean_alpha_sq)/(mean_alpha_sq + mean_eps_sq/10)

gen FE_bayes = mu_hat + lambda*(Y_bar_j - mu_hat)

preserve
collapse (first) FE implied_mu mu_j implied_mu_2 FE_bayes, by(j)

twoway pcarrow FE mu_j FE_bayes mu_j, ytitle("Standard mean/Emprical Bayes mean") xtitle("True estimand")

graph export "bilder\scatterFE2.png", replace
restore





**# Q2


foreach T in 200 5{
	quietly{
	clear all

	local T = `T'
	local N = 10 

	local total = `T'*`N'

	set obs `total'

	set seed 400

	gen i = mod(_n-1,10) + 1
	bys i: gen t = _n
	gen random_m = runiform(1,21)
	gen j = floor(random_m)
	drop random_m


	bys i: gen temp = rnormal(0,1) if _n == 1
	egen a_i = max(temp), by(i)
	drop temp

	bys j: gen temp = rnormal(1,1) if _n == 1
	egen  q_j = max(temp), by(j)
	drop temp

	bys t: gen temp = rnormal(1,1) if _n == 1
	egen  g_t = max(temp), by(t)
	drop temp

	gen eps_it = rnormal(0,0.2)

	sort t i

	gen Y_it = a_i + g_t + q_j + eps_it
	
	gen est_a = 0
	gen est_q = 0

	reg Y_it i.i i.j 
	
	local j = 0
	
	forvalues c = 1/30{
		
		if `c'<11{
			
			replace est_a = r(table)[1, `c'] if i == `c'
			
			
		}
		
		else{
			local j = `j' + 1
			
			replace est_q = r(table)[1,`c'] if j == `j'
		}
		
		
	}
	
	gen temp = a_i if i == 1
	egen temp2 = max(temp)
	gen relative_a = a_i - temp2 
	drop temp temp2
	
	gen temp = q_j if j == 1
	egen temp2 = max(temp)
	gen relative_q = q_j - temp2
	drop temp temp2
	
	
	
	twoway (scatter est_a relative_a) (line relative_a relative_a), name(ind, replace)
	twoway (scatter est_q relative_q) (line relative_q relative_q), name(machine, replace)
	
	graph combine ind machine
	
	graph export "bilder\N_`T'.png", replace
	}
}


cd "C:\Users\olofd\PhD\Mertrics 2\PS\4"

clear all 

local theta = 1

	
	
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

frame change default

	
set seed 15
forvalues p = 1/11{
	
	quietly{
	local theta = -1 + (_n-1)*0.2
	

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
	
	gen Y_ij = Y_00 + D*(Y_10 - Y_00) + W*(Y_01 - Y_00) + D*W*(Y_11 - Y_00) 
	
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
		
		
		
	}
	
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
		
		graph twoway (rcap upper lower theta) (scatter `t'_`i' theta) (scatter tao_`t' theta)
		
		graph export "bilder\\`t'_`i'.png", replace
		
		drop upper lower
	}
	
}

clear all

set seed 17
frame change default 
 
local theta = 0
	

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
	
gen Y_ij = Y_00 + D*(Y_10 - Y_00) + W*(Y_01 - Y_00) + D*W*(Y_11 - Y_00) 
	
gen tao_D = (Y_10 - Y_00)
gen tao_W = (Y_01 - Y_00)
 
reg Y_ij W i.j 
reg Y_ij D W

summarize tao_W

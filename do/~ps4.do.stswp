cd "C:\Users\olofd\PhD\Mertrics 2\PS\4"

clear all 

local theta = 1

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
	
	reg Y_ij i.D i.W
	
	matrix list r(table)
	
	reg Y_ij i.D##i.W
	
	local i = 0
	foreach v in " " "D_W"{
		local `i' = `i' + 1
		
		reg Y_ij D W `v', robust
		
		
		
	}
	
	
	
	


frame create measures
frame change measures
set obs 11

gen theta = -1 + (_n-1)*0.2



foreach t in " " "V"{
	forvalues i = 1/2{
		gen `t'D_`i' = 0
		gen `t'W_`i' = 0
	}
	
}

gen tau_d = 0
gen tau_w = 0


	
	
forvalues theta = -1(0.2)1{
	local theta = round(`theta',0.1)
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
	
	gen X_ij = rnormal(`theta', 0.2)
	
	gen Y_01 = rnormal(mu_j,1)
	
	gen Y_11 = Y_10 + Y_01 + X_ij
	
	gen Y_ij = Y_00 + D*(Y_10 - Y_00) + W*(Y_01 - Y_00) + D*W*(Y_11 - Y_00) 
	
	local i = 0
	foreach v in " " "D_W"{
		`i' = `i' + 1
		
		reg Y_ij D W D_W, robust
		
		local D_`i' = r(table)[1,1]
		local VD_`i' = r(table)[2,1]
		
		local W_`i' = r(table)[1,2]
		local VW_`i' = r(table)[2,2]
		
	}
	
	frame change measures
	
	
	
}



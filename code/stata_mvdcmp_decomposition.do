// prepare data

use "2017-08-28-ideath.dta", clear
gen id = _n
replace survtime_h = 0.1 if survtime_h == 0
stset survtime_h, failure(death) id(id)
replace _st = 1 if _st == 0

// 4 different setups for 4 different stata sessions

    // hour-day
	
	// cut time into 3 intervals: [0, 1), [1, 24), [24, Inf)    
    stsplit time, at (1, 24)
	// calculate exposure time within interval as difference
	// between exit and entry times
    gen exposure = _t-_t0
	// disregard interval [24, Inf)
    drop if time == 24 
    // replace time = 0 if time == 0
    // replace time = 1 if time == 1
	// sum up death counts and exposures across intervals by strata
    collapse (sum) deaths = _d person_hours = exposure, by(time birthweight_c congenital_anomalies apgar5 age_of_mother_c sex education_of_mother race_and_hispanic_orig_of_mother)
    
    // day-week
    
    stsplit time, at (1, 24, 168)
    gen exposure = _t-_t0
    drop if time == 0
    drop if time == 168
    replace time = 0 if time == 1
    replace time = 1 if time == 24
    collapse (sum) deaths = _d person_hours = exposure, by(time birthweight_g_c congenital_anomalies apgar5 age_of_mother_c sex education_of_mother race_and_hispanic_orig_of_mother)
    
    // week-month
    
    stsplit time, at (24, 168, 672)
    gen exposure = _t-_t0
    drop if time == 0
    drop if time == 672
    replace time = 0 if time == 24
    replace time = 1 if time == 168
    collapse (sum) deaths = _d person_hours = exposure, by(time birthweight_g_c congenital_anomalies apgar5 age_of_mother_c sex education_of_mother race_and_hispanic_orig_of_mother)
    
    // month-year
    
    stsplit time, at (168, 672)
    gen exposure = _t-_t0
    drop if time == 0
    replace time = 0 if time == 168
    replace time = 1 if time == 672
    collapse (sum) deaths = _d person_hours = exposure, by(time birthweight_g_c congenital_anomalies apgar5 age_of_mother_c sex education_of_mother race_and_hispanic_orig_of_mother)

// decompose

poisson deaths i.time i.birthweight_g_c i.congenital_anomalies i.apgar5 i.age_of_mother_c i.sex i.education_of_mother i.race_and_hispanic_orig_of_mother, irr exposure(person_hours)

margins i.time, predict(ir)
gen rate = deaths/person_hours
tabstat rate [aw=person_hours], by(time)

gen log_person_hours = log(person_hours)

mvdcmp time, reverse: poisson deaths i.time i.birthweight_c i.congenital_anomalies i.apgar5 i.age_of_mother_c i.sex i.education_of_mother i.race_and_hispanic_orig_of_mother, offset(log_person_hours)


// test
    // cut time into 3 intervals: (0, 0.1], (0.1, 1], (1, Inf)
    // which given our coding actually translates into
	// [0, 1), [1, 24), [24, Inf)    
    stsplit time, at (0.1, 1)
	// calculate exposure time within interval as difference
	// between exit and entry times
    gen exposure = _t-_t0
	// disregard interval [24, Inf)
    drop if time == 1
	replace time = 0 if time == 0
    gen time2 = 1 if time == 0.1
	// sum up death counts and exposures across intervals by strata
    collapse (sum) deaths = _d person_hours = exposure, by(time apgar5)

gen log_person_hours = log(person_hours)
mvdcmp time, reverse: poisson deaths i.time i.apgar5, offset(log_person_hours)

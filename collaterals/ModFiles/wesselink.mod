TITLE Sensory Axon Node channels

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX wesselink
	NONSPECIFIC_CURRENT ina
	NONSPECIFIC_CURRENT ik

	RANGE gnapbar, gnabar, gkbar, gl, gkf, ena, ek, el, ekf, pna
	RANGE mp_inf, m_inf, h_inf, s_inf, n_inf
	RANGE tau_mp, tau_m, tau_h, tau_s, tau_n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

    : channel conductances
	pna	= .0704 
	gkbar     = 300 (mho/m2)

    : reversal potentials
	ena     = 50.0  (mV)
	ek     = -84.0  (mV)

    : variables read in from .hoc file
	celsius		(degC)
	dt              (ms)
	v               (mV)
	vtraub=-80

    : parameters determining rate constants

    F = 96485
    R = 8.3144

  	na_o = 154
  	na_i = 30
}

STATE {
	m h n
}

ASSIGNED {
	ina	(mA/cm2)
	ik	(mA/cm2)
	mp_inf
	m_inf
	h_inf
	n_inf
	tau_m
	tau_h
	tau_n
	q10_1
	q10_2
	q10_3
	T
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina =  0.1 * pna * m*m*m*h * v * (F * F) / (R * T) * (na_o - na_i*exp(v*F / (R*T))) / (1 - exp((v*F)/ (R*T))) : fast sodium, 0.1 converts to ma/cm2
	ik =  0.1 * gkbar * n*n*n*n * (v - ek): fast potassium, 0.1 converts to ma/cm2
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
       evaluate_fct(v)
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
	n' = (n_inf - n) / tau_n
}

UNITSOFF

INITIAL {
:
:	Q10 adjustment
:   Temperature dependence
:
	
	T = celsius + 273.15
	q10_1 = 2.2 ^ ((celsius-20)/ 10 )
	q10_2 = 2.9 ^ ((celsius-20)/ 10 )
	q10_3 = 3.0 ^ ((celsius-37)/ 10 )

	evaluate_fct(v)
	m = m_inf
	h = h_inf
	n = n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

    : fast Na+
	a = q10_3*vtrap6(v)
	b = q10_3*vtrap7(v)
	tau_m = 1 / (a + b)
	m_inf = a / (a + b)

	a = q10_2*vtrap8(v)
	b = q10_2*vtrap9(v)
	tau_h = 1 / (a + b)
	h_inf = a / (a + b)

    : fast K+
	a = q10_3*vtrap10(v)
	b = q10_3*vtrap11(v)
	tau_n = 1 / (a + b)
	n_inf = a / (a + b)

}

: vtrap functions to prevent discontinuity


FUNCTION vtrap6(x) {
	vtrap6 =  4600*(x + 18.4) / (1 - Exp( (-18.4 - x)/10.3 ))
}

FUNCTION vtrap7(x) {
	vtrap7 =  330*(-22.7 - x) / (1 - Exp( (x + 22.7)/9.16 ))
}

FUNCTION vtrap8(x) {
	vtrap8 =  210*(-111 - x) / (1 - Exp( (x + 111)/11 ))
}

FUNCTION vtrap9(x) {
	vtrap9 =  14100 / (1 + Exp( (-28.8 - x)/1.1 ))
}

FUNCTION vtrap10(x) {
	vtrap10 =  51.7*(x + 93.2) / (1 - Exp( (-93.2 - x)/1.1 ))
}

FUNCTION vtrap11(x) {
	vtrap11 =  92*(-76 - x) / (1 - Exp( (x + 76)/10.5 ))
}


FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}

UNITSON

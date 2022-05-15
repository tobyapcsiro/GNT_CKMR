//////////////////////////////////////////////////////
// glyphis NT rivers CKMR code ///////////////////////
// T. Patterson, R. Hillary CSIRO 2017 ///////////////
//////////////////////////////////////////////////////

#include <TMB.hpp>

template<class Type> Type inv_logit(Type x) {return(Type(1)/(Type(1)+exp(-x)));}

template<class Type>
Type objective_function<Type>::operator() ()
{
	PARAMETER_ARRAY(xi);		    // s x r - river probs parameters
	PARAMETER_VECTOR(lambda);		// r	 - river growth rate
	PARAMETER_VECTOR(lN0);			// r	 - log N0 river
	PARAMETER_VECTOR(iphi); 		// r   - probability by river.
	PARAMETER_VECTOR(lnu); 	    // litter effect (constant over rivers?)
  PARAMETER_VECTOR(izeta);    // sex-ratio parameter length r
  PARAMETER(itheta);          // multiple paternity parameter (relates to females)
  PARAMETER(lgamma);          // multiple female breeding partners per year (relates to males)

	DATA_INTEGER(nsex);
	DATA_INTEGER(nriv);
	DATA_INTEGER(tmax);  			// number of cohorts covered including tzero
  DATA_INTEGER(tzero); 			// cohort corresponding to time 0
	DATA_SCALAR(p_hsp); 			// critical PLOD (0.92).

	DATA_IVECTOR(c1); 				// cohort 1.
	DATA_IVECTOR(c2); 				// cohort 2
	DATA_IVECTOR(r1); 				// riv 1
	DATA_IVECTOR(r2); 				// riv 2
	DATA_IVECTOR(kcode); 			// ktype (0 UP, 1 HSP, 2 FSP)
	DATA_IVECTOR(k); 				  // number of kin in that comparison grouping
	DATA_VECTOR(plod); 				// PLOD for mtDNA stuff later on
  DATA_IVECTOR(h1);         // haplotype of fish 1
  DATA_IVECTOR(h2);         // haplotype of fish 2
  DATA_VECTOR(ph1);         // haplotype frequency of fish 1's haplo
  DATA_VECTOR(ph2);         // haplotype frequency of fish 2's haplo

	int nobs = r1.size();

	array<Type> N(nsex, nriv,tmax);
	N.setZero();

	array<Type> Ntil(nsex, nriv,tmax);
	Ntil.setZero();

	array<Type> omega(nsex,nriv,nriv);

	Type theta = inv_logit(itheta);
  Type gamma = exp(lgamma);

	vector<Type> nu(nriv);
	vector<Type> zeta(nriv);

	for(int i=0; i<nriv;i++) {

		nu(i) = exp(lnu(i));
		zeta(i) = inv_logit(izeta(i));

	}

	vector<Type> phi(nriv);
	array<Type> sexrat(nsex,nriv);

	// indices.
	// s = 1,2 = fem, male.
	// r = 1,2 = river1, river2 specified as whatevs (e.g. ad, al)

	//  POP DYNAM /

  // set up the sex-specific parameters
		for(int r=0; r < nriv; r++)	phi(r) = inv_logit(iphi(r));
	

  for(int s=0; s < nsex; s++) {

      omega(s,0,0) = inv_logit(xi(s,0));
      omega(s,0,1) = Type(1)-omega(s,0,0);
      omega(s,1,1) = inv_logit(xi(s,1));
      omega(s,1,0) = Type(1)-omega(s,1,1);

  }

  for(int r=0; r < nriv; r++) {
    sexrat(0,r) = zeta(r);
    sexrat(1,r) = Type(1)-zeta(r);
  }

  vector<Type> N0 = exp(lN0);
	for(int t=0; t < tmax; t++) { 			// time
		for(int r=0; r < nriv; r++){		// river
			for(int s=0; s < nsex; s++) {	// sex
			   N(s,r,t) = N0(r) * sexrat(s,r) * exp(lambda(r) * t);
				// calc repro adults in river at time
				for(int rp=0; rp < nriv; rp++){
					Ntil(s,r,t) += N(s,rp,t) * omega(s,r,rp);
				}
			}
		}
	}

	//  CK LIKE CALCS
	Type logl = Type(0);
  Type logl_ckmr = Type(0);
  Type logl_mtdna = Type(0);
  Type pkin;
  Type pfsp;
  Type pmhsp;
  Type pphsp;
  Type psum;
  vector<Type> praw(3);
  vector<Type> pmtdna(3); 
  int cmin,cmax,rmin,rmax;

  for(int i=0; i < nobs; i++) {           // loop over obs.

    // within-cohort cases

    if(c1(i)==c2(i)) {// within cohort comparisons

      cmin = c1(i);
      cmin -= tzero;

      // FSPs

      if(r1(i)==r2(i)) pfsp = (nu(r1(i)-1) * (Type(1)-theta))/ Ntil(0, r1(i)-1, cmin);
			else pfsp = Type(0);

      // MHSPs

      if(r1(i)==r2(i)) pmhsp = (p_hsp * nu(r1(i)-1) * theta)/ Ntil(0, r1(i)-1, cmin);
			else pmhsp = Type(0);

      // PHSPs

      if(r1(i)==r2(i)) pphsp = (p_hsp * gamma)/ Ntil(1, r1(i)-1, cmin);
			else pphsp = Type(0);

      pkin = pfsp + pmhsp + pphsp;

      if(pkin > Type(0)) {

        if(k(i) == 0) logl_ckmr += log(Type(1)-pkin);
        if(k(i) == 1) {
          
          logl_ckmr += log(pkin);
          
          // mtDNA bit needs calculating for an identified kin pair (PLOD > eta)

          praw(0) = pfsp;
          praw(1) = pmhsp + pphsp * ph2(i);
          praw(2) = pphsp * (Type(1)-ph2(i));
          psum = praw.sum();
          for(int m=0; m<3; m++) pmtdna(m) = praw(m)/psum;

          // is it an FSP?
          if(kcode(i) == 2) logl_mtdna += log(pmtdna(0));
          // is it an HSP which shares a haplo?
          if(kcode(i) == 1 & h2(i) == h1(i)) logl_mtdna += log(pmtdna(1));
          // is it an HSP which doesn't share a haplo?
          if(kcode(i) == 1 & h2(i) != h1(i)) logl_mtdna += log(pmtdna(2)); 

        }

      }
		}

    // across-cohort case

		if ( c1(i)!=c2(i) ) {

      cmin = c1(i) < c2(i) ? c1(i) : c2(i);
      cmax = c1(i) < c2(i) ? c2(i) : c1(i);
      rmin = c1(i) < c2(i) ? r1(i) : r2(i);
      rmax = c1(i) < c2(i) ? r2(i) : r1(i);
      cmin -= tzero;
      cmax -= tzero;
      rmin--;
      rmax--;

      // FSPs

      pfsp = Type(0);

      // MHSPs

      pmhsp = Type(0);
      for(int rp=0; rp < nriv; rp++) pmhsp += (N(0,rp,cmin) * omega(0,rmin,rp))/ (Ntil(0,rmin,cmin)) * pow(phi(rp),cmax-cmin) * omega(0,rmax,rp)/(Ntil(0,rmax,cmax));
      pmhsp *= p_hsp;

      // PHSPs

      pphsp = Type(0);
      for(int rp=0; rp < nriv; rp++) pphsp += (N(1,rp,cmin) * omega(1,rmin,rp))/ (Ntil(1,rmin,cmin)) * pow(phi(rp),cmax-cmin) * omega(1,rmax,rp)/(Ntil(1,rmax,cmax));
      pphsp *= p_hsp;

      pkin = pfsp + pmhsp + pphsp;

      if(pkin > Type(0)) {

        if(k(i) == 0) logl_ckmr += log(Type(1)-pkin);
        if(k(i) == 1) {
          
          logl_ckmr += log(pkin);

          // mtDNA bit needs calculating for an identified kin pair (PLOD > eta)

          praw(0) = pfsp;
          praw(1) = pmhsp + pphsp * ph2(i);
          praw(2) = pphsp * (Type(1)-ph2(i));
          psum = praw.sum();
          for(int m=0; m<3; m++) pmtdna(m) = praw(m)/psum;

          // is it an FSP?
          if(kcode(i) == 2) logl_mtdna += log(pmtdna(0));
          // is it an HSP which shares a haplo?
          if(kcode(i) == 1 & h2(i) == h1(i)) logl_mtdna += log(pmtdna(1));
          // is it an HSP which doesn't share a haplo?
          if(kcode(i) == 1 & h2(i) != h1(i)) logl_mtdna += log(pmtdna(2));  
        }
      }

    }
  }

  // overall log-likelihood

  logl = logl_ckmr + logl_mtdna;

  // report stuff on natural scale
  ADREPORT(N0);
  ADREPORT(nu);
  ADREPORT(zeta);
  ADREPORT(theta);
  ADREPORT(gamma);
  ADREPORT(omega); 

  return(-logl);
}

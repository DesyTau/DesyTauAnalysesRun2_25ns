// author : F. Costanza
// date : 20.01.16
//
// Function to calculte weights for HT binned samples
// Returned weights must be multiplied by instantaneous luminosity[pb-1].


enum samples{ DYM5to50=1, DYM50, W};

double DYM5to50_HTBinned_Weights( float ht){
  float xslo_incl = 71310.;
  float xslo_binned = 0.;
  float k = 1.;

  float n_incl = 9404398.;
  float n_bin = 0.;

  if(ht < 100.){
    xslo_binned = 71043.;
    n_bin = xslo_binned / xslo_incl * n_incl;
  }
  else if (ht < 200.){
    xslo_binned = 224.2;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 1013479.;
  }
  else if (ht < 400.){
    xslo_binned = 37.2;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 1011756.;
  }
  else if (ht < 600.){
    xslo_binned = 3.581;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 998751.;
  }
  else{
    xslo_binned = 1.124;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 1007309.;
  }
  
  return (k * xslo_binned / n_bin);
}
double DYM50_HTBinned_Weights( float ht){
  float xslo_incl = 4895.;
  float xslo_binned = 0.;
  float xslnnlo_incl = 6025.2;
  float k = xslnnlo_incl / xslo_incl;

  float n_incl = 9052671;
  float n_bin = 0.;

  if(ht < 100.){
    xslo_binned = 4705.;
    n_bin = xslo_binned / xslo_incl * n_incl;
  }
  else if (ht < 200.){
    xslo_binned = 139.4;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 2725655.;
  }
  else if (ht < 400.){
    xslo_binned = 42.75;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 973937.;
  }
  else if (ht < 600.){
    xslo_binned = 5.497;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 1067758.;
  }
  else{
    xslo_binned = 2.21;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 998912.;
  }
  
  return (k * xslo_binned / n_bin);
}

double W_HTBinned_Weights( float ht){
  float xslo_incl = 50690.;
  float xslo_binned = 0.;
  float xslnnlo_incl = 61526.7;
  float k = xslnnlo_incl / xslo_incl;

  float n_incl = 72207128;
  float n_bin = 0.;

  if(ht < 100.){
    xslo_binned = 48917.6;
    n_bin = xslo_binned / xslo_incl * n_incl;
  }
  else if (ht < 200.){
    xslo_binned = 1345.;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 10152718.;
  }
  else if (ht < 400.){
    xslo_binned = 359.7;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 5221599.;
  }
  else if (ht < 600.){
    xslo_binned = 48.91;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 1745914.;
  }
  else{
    xslo_binned = 18.77;
    n_bin = xslo_binned / xslo_incl * n_incl;
    n_bin += 1039152.;
  }
  
  return (k * xslo_binned / n_bin);
}

double HTBinnedWeights( int sample, float ht){
  if(sample == DYM5to50)
    return DYM5to50_HTBinned_Weights(ht);
  else if(sample == DYM50)
    return DYM50_HTBinned_Weights(ht);
  else if(sample == W)
    return W_HTBinned_Weights(ht);
  else
    return 1.;
}


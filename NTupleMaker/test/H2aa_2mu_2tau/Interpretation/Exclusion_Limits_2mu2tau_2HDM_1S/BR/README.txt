The file 

BR_type_tan(beta).dat 

contains the relevant branching ratios (BRs) of the pseudoscalar a in the 2HDM+S model. Here "type" denotes the type of the 2HDM, i.e. I, II, III or IV, and "tan(beta)" indicates the value of tan(beta) that has been used to obtain the BRs --- in the case of the type-I model the BRs do not depend on tan(beta) and therefore only the file BR_I.dat is provided. The calculation of the BRs is based on the formulas presented in Appendix A and B of 

https://arxiv.org/pdf/1802.02156.pdf

which should be cited if using the data files. 

Each line in each data file has 11 entries. These entries are: 

1) type of the 2HDM 
2) mass of the pseudoscalar a in GeV
3) tan(beta) value
4) BR(a -> bb)
5) BR(a -> cc)
6) BR(a -> tautau)
7) BR(a -> mumu)
8) BR(a -> gluon gluon) 
9) BR(a -> photon photon)
10) BR(a -> ss + dd + uu)
11) BR(a -> hadronic) 

We emphasise that BR(a -> ss + dd + uu) only includes "perturbative" contributions to the decay channels a -> ss, dd and uu but no non-perturbative effects or contributions that are associated to the mixing of the pseudoscalar a with QCD resonances. Since non-perturbative and mixing effects are generically important for m_a values below the bb-threshold, the given BR(a -> ss + dd + uu) numbers are in practice of little use. 

Non-perturbative and mixing effects in the decays of the pseudoscalar a to light hadrons are instead (indirectly) included in BR(a -> hadronic), since this quantity is defined as follows 

BR(a -> hadronic) = 1 - sum_{x = b, c, tau, mu, gluon, photon} BR(a -> xx) 

















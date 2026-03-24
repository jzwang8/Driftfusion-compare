
Ef_Eg_p = -.89;
Ef_Eg_n = -.026;

[Nd_p,Na_p] = getDopingConcFromEf(Ef_Eg_p);
[Nd_n,Na_n] = getDopingConcFromEf(Ef_Eg_n);

[mu_n_p,mu_p_p] = getMobilitiesFromDopingConc(Na_p)
[mu_n_n,mu_p_n] = getMobilitiesFromDopingConc(Nd_n)
function alpha_norm = NormalizeAbsorbance(lambda, alpha, laserpower)
    Isun = lightsource('AM15', lambda);
    Ilaser = laserpower * 1e-3;
    dK = Isun ./ Ilaser .* alpha;
    K = trapz(lambda,dK);
    alpha_norm = alpha ./ K;
end
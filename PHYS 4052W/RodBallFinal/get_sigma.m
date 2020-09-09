function s = get_sigma(obtained, accepted, obtained_err)
    p = normcdf(obtained, accepted, obtained_err);
    lower = p/2;
    upper = 1 - (p/2);
    bounds = norminv([lower upper]);
    s = bounds(2);
end


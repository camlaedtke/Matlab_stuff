function vals = ChiVals(y, fit_y, y_err)
   chi = (y - fit_y)./(y_err);
   chi2 = chi.^2;
   vals = [chi, chi2];
end

  function thechi=chilinear(x,y,sigma,m,b)
  chis=(y-(m .* x + b)) ./ sigma;
  chi2s=chis.^2;
  thechi=sum(chi2s);


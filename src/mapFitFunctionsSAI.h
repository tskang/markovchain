List _mcFitMap(CharacterVector stringchar, bool byrow, double confidencelevel, NumericMatrix hyperparam) {
  // get mapEstMatr and freqMatr 
  CharacterVector elements = stringchar;
  elements = unique(elements).sort();
  int sizeMatr = elements.size();
  
  // if no hyperparam argument provided, use default value of 1 for all 
  if(hyperparam.nrow() == 1 && hyperparam.ncol() == 1){
    NumericMatrix temp(sizeMatr, sizeMatr);
    for(int i = 0; i < sizeMatr; i++)
      for(int j = 0; j < sizeMatr; j++)
        temp(i, j) = 1;
    hyperparam = temp;
  }
  
  // validity check
  if(hyperparam.nrow() != sizeMatr || hyperparam.ncol() != sizeMatr) 
    stop("Dimensions of the hyperparameter matrix are inconsistent");
  
  NumericMatrix mapEstMatr(sizeMatr), expMatr(sizeMatr);
  NumericMatrix freqMatr(sizeMatr);
  mapEstMatr.attr("dimnames") = List::create(elements, elements); 
  expMatr.attr("dimnames") = List::create(elements, elements); 

  NumericMatrix lowerEndpointMatr = NumericMatrix(mapEstMatr.nrow(), mapEstMatr.ncol());
  NumericMatrix upperEndpointMatr = NumericMatrix(mapEstMatr.nrow(), mapEstMatr.ncol());
  NumericMatrix varianceMatr = NumericMatrix(mapEstMatr.nrow(), mapEstMatr.ncol());

  // populate frequeny matrix for old data; this is used for inference 
  int posFrom, posTo;
  for(int i = 0; i < stringchar.size() - 1; i ++) {
    for (int j = 0; j < sizeMatr; j ++) {
      if(stringchar[i] == elements[j]) posFrom = j;
      if(stringchar[i + 1] == elements[j]) posTo = j;
    }
    freqMatr(posFrom,posTo)++;
  }
 
  // sanitize and to row probs
  for (int i = 0; i < sizeMatr; i++) {
    double rowSum = 0, paramRowSum = 0;
    for (int j = 0; j < sizeMatr; j++)
      rowSum += freqMatr(i, j), paramRowSum += hyperparam(i, j);
    
    // toRowProbs
    for (int j = 0; j < sizeMatr; j++) {
      // confidence intervals and bounds
      double p = freqMatr(i, j) + hyperparam(i, j), q = rowSum + paramRowSum - freqMatr(i, j) - hyperparam(i, j);
      
      // expected value of the transition parameters
      expMatr(i, j) = p / (p + q);
      
      if(p + q == sizeMatr)
              mapEstMatr(i, j) = 1 / sizeMatr;
      else
              // maximum a posteriori estimate
              mapEstMatr(i, j) = (p - 1) / (p + q - sizeMatr);
              
      double beta = lbeta(p, q);
      double cdf = betain(double(mapEstMatr(i, j)), p, q, beta);

      if(cdf + confidencelevel / 2 > 1.){
        upperEndpointMatr(i, j) = 1.;
        lowerEndpointMatr(i, j) = xinbta(p, q, beta, 1 - confidencelevel);
      }
      else if(cdf - confidencelevel / 2 < 0.){
        lowerEndpointMatr(i, j) = 0.;
        upperEndpointMatr(i, j) = xinbta(p, q, beta, confidencelevel);
      }
      else{
        lowerEndpointMatr(i, j) = xinbta(p, q, beta, cdf - confidencelevel / 2);
        upperEndpointMatr(i, j) = xinbta(p, q, beta, cdf + confidencelevel / 2);
      }

      varianceMatr(i, j) = p * q / (p + q) / (p + q) / (1 + p + q);
    }
  }

  if(byrow==false) mapEstMatr = _transpose(mapEstMatr);

  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = mapEstMatr;
  outMc.slot("name") = "Bayesian Fit";  
  
  return List::create(_["estimate"] = outMc
    , _["confidenceInterval"] = List::create(_["confidenceLevel"]=confidencelevel, 
              _["lowerEndpointMatrix"]=lowerEndpointMatr, 
              _["upperEndpointMatrix"]=upperEndpointMatr),
              _["expectedValue"]=expMatr,
              _["varianceMatrix"]=varianceMatr
  );
}

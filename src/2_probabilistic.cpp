// [[Rcpp::depends(RcppArmadillo)]]

//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::vec colSums(arma::mat X){
   int nCols = X.n_cols;
   arma::vec out(nCols);
   for(int i = 0; i < nCols; i++){
      out(i) = sum(X.col(i));
   }
   return(out);
}

//List commclassesKernel(NumericMatrix P){
// [[Rcpp::export(.commclassesKernelRcpp)]]
extern "C" SEXP commclassesKernel(NumericMatrix P){
//  Rf_PrintValue(P);
  unsigned int m = P.ncol(), n;
  CharacterVector stateNames = rownames(P);
  std::vector<int> a;
  arma::vec b, c, d;
  arma::mat T = arma::zeros(m, m);
//  std::cout << T << std::endl;
  unsigned int i = 0;
  int oldSum, newSum, ai;
  while(i < m) {
    a.resize(0);
    a.push_back(i);
//    a.resize(1);
//    a[0] = i;
//    std::cout << "i = " << i << ", a=" << a ;
//    b = NumericVector(m);
    b = arma::zeros<arma::vec>(m);
    b[i] = 1;
//    std::cout << "b = " << b << ", b[i]=" << b[i] << std::endl;
    newSum = 0;
    oldSum = 1;
//    bool first = true;
    while(oldSum != newSum) {
      oldSum = 0;
      for(int j = 0; j < b.size(); j ++)
        if(b[j] > 0) oldSum += (j + 1);
//      std::cout << "b = " << b << std::endl;
//      std::cout << "find(b>0) = " << find(b > 0) << std::endl;
      n = a.size();
//      std::cout << "oldSum = " << oldSum << ", n=" << n << std::endl;
      NumericVector temp; 
      NumericMatrix matr(n, m);
      for(unsigned int j = 0; j < n; j ++) {
        temp = P.row(a[j]);
//        Rcout << "temp a[" << j << "] " << a[j] << " = ";
//        Rf_PrintValue(temp);
        for(int k = 0; k < temp.size(); k++) 
          matr(j, k) = temp[k];
      }
//      Rcout << "matr = " << std::endl;
//      Rf_PrintValue(matr);
      c = arma::zeros<arma::vec>(m);
      for(int j = 0; j < m; j++) 
        for(int k = 0; k < n; k++)
          c[j] += matr(k, j);
//      c = colSums(matr);
//      Rcout << "colSums = " << c << std::endl;
      newSum = 0;
      a.resize(0);
//      Rcout << "a.size = " << a.size() << std::endl;
      for(int j = 0; j < b.size(); j++) {
        if(c[j] > 0) {
          b[j] = 1; a.push_back(j);
        }
        if(b[j] > 0) newSum += (j + 1);
      }
//      Rcout << "a.size = " << a.size() << std::endl;
    }
    for(unsigned int j = 0; j < b.size(); j ++)
      T(i, j) = b[j];
    i++;
  }
//  std::cout << "T = " << T << std::endl;
  arma::mat F = arma::trans(T);
//  std::cout << "T = " << T.n_rows << " x " << T.n_cols << std::endl;
//  std::cout << "F = " << F.n_rows << " x " << F.n_cols << std::endl;
  LogicalMatrix C;
  arma::mat Ca(T.n_rows, T.n_cols);// = (T > 0);// & (F > 0);
//  std::cout << "C created: " << Ca.n_rows << " x " << Ca.n_cols << std::endl;
  for(i = 0; i < T.n_rows; i ++) {
//   std::cout << "i = " << i << std::endl;
   for(unsigned int j = 0; j < T.n_cols; j++) {
      Ca(i, j) = (T(i, j) > 0 && F(i, j) > 0);
   }
//   std::cout << "i end = " << i << std::endl;
  }
//  std::cout << "C filled" << std::endl;
  LogicalVector v(T.n_cols);
  arma::mat tC = Ca.t();
  arma::mat tT = T.t();
  unsigned int sums[tC.n_cols];
  for(unsigned int j = 0; j < T.n_cols; j++) {
    sums[j] = 0;
    for(i = 0; i < T.n_rows; i ++)
      if(tC(i, j) == tT(i, j)) sums[j] ++;
    v[j] = (sums[j] == m);
  }
//  std::cout << "v created" << std::endl; 
//  std::cout << "Ca: " << Ca << std::endl;
  C = as<LogicalMatrix>(wrap(Ca));
  C.attr("dimnames") = List::create(stateNames, stateNames);
  v.names() = stateNames;
//  Rf_PrintValue(C);
//  Rf_PrintValue(v);
  return List::create(_["C"] = C, _["v"] = v);
}

//returns the underlying communicating classes
// [[Rcpp::export(.communicatingClassesRcpp)]]
List communicatingClasses(LogicalMatrix adjMatr)
{
  int len = adjMatr.nrow();
  List classesList;
//  Rf_PrintValue(adjMatr);
//  Rf_PrintValue(rownames(adjMatr));
  CharacterVector rnames = rownames(adjMatr);
  for(int i = 0; i < len; i ++) {
    bool isNull = false;
    LogicalVector row2Check = adjMatr(i, _);
//    Rf_PrintValue(row2Check);    
    CharacterVector proposedCommClass;
    for(int j = 0; j < row2Check.size(); j++) {
      if(row2Check[j] == true) {
        String rname = rnames[j];
//        Rcout << (std::string)rname << std::endl;
        proposedCommClass.push_back(rname);
      }
    }
    if (i > 0) {
      for(int j = 0; j < classesList.size(); j ++) {
        bool check = false;        
        CharacterVector cv = classesList[j];
//        Rf_PrintValue(cv);
        std::set<std::string> s1, s2;
//        std::cout << "before insert: " << cv.size() << " " << proposedCommClass.size() << std::endl;
        for(int k = 0; k < cv.size(); k ++) {
          s1.insert(as<std::string>(cv[k]));
//          std::cout << "s1 inserted " << cv[k] << std::endl;
          if(proposedCommClass.size() > k) {
            s2.insert(as<std::string>(proposedCommClass[k]));
//            std::cout << "s2 inserted " << proposedCommClass[k] << std::endl;
          }
        }
//        for(std::set<std::string>::iterator it = s1.begin(); it != s1.end(); it++)
//          std::cout << *it << " ";
//        std::cout << std::endl;
//        for(std::set<std::string>::iterator it = s2.begin(); it != s2.end(); it++)
//          std::cout << *it << " ";
//        std::cout << std::endl;
        check = std::equal(s1.begin(), s1.end(), s2.begin());
//        std::cout << check << std::endl;
        if(check) {
          isNull = true;
          break;
        }
      }
    }
//    if(!Rf_isNull(proposedCommClass) && proposedCommClass.size() > 0) 
    if(!isNull) 
      classesList.push_back(proposedCommClass);
  }
//  Rf_PrintValue(classesList);
  return classesList;
}

arma::mat _pow(arma::mat A, int n) {
  arma::mat R = arma::eye(A.n_rows, A.n_rows);
  for(int i = 0; i < n; i ++) 
    R = A*R;
  return R;
}

//communicating states
// [[Rcpp::export(.commStatesFinderRcpp)]]
NumericMatrix commStatesFinder(NumericMatrix matr)
{
  //Reachability matrix
  int dimMatr = matr.nrow();
  arma::mat X(matr.begin(), dimMatr, dimMatr, false);
  arma::mat temp = arma::eye(dimMatr, dimMatr) + arma::sign(X);
  temp = _pow(temp, dimMatr - 1);
  arma::mat m;
  arma::mat R = arma::sign(temp);
  return wrap(R);
}
/*
#@ Tae: to be fully moved in Rcpp
#communicating states
.commStatesFinder<-function(matr)
{
  #Reachability matrix
  dimMatr<-dim(matr)[1]
  Z<-sign(matr)
  temp<-(eye(dimMatr)+Z)%^%(dimMatr-1)
  R<-sign(temp)
  return(R)
}
*/

bool _intersected(CharacterVector v1, CharacterVector v2) {
  CharacterVector::iterator first1 = v1.begin();
  CharacterVector::iterator last1 = v1.end();
  CharacterVector::iterator first2 = v2.begin();
  CharacterVector::iterator last2 = v2.end();
  while(first1!=last1 && first2!=last2) {
    if(*first1 == *first2) return true;
    else if(*first1 < *first2) ++first1;
    else ++first2;    
  }
  return false;
}

// [[Rcpp::export(.summaryKernelRcpp)]]
List summaryKernel(S4 object)
{
  NumericMatrix matr = object.slot("transitionMatrix");
  List temp = commclassesKernel(matr);
  List communicatingClassList = communicatingClasses(temp["C"]);
  List v = temp["v"];
  CharacterVector ns = v.names();
  CharacterVector transientStates; //<-names(which(temp$v==FALSE))
  for(int i = 0; i < v.size(); i++) {
    if(v[i] == false)
      transientStates.push_back(ns[i]);
  }
  List closedClasses, transientClasses;

  for(int i = 0; i < communicatingClassList.size(); i ++)
  {
    CharacterVector class2Test = communicatingClassList[i];
    if(_intersected(class2Test,transientStates)) 
        transientClasses.push_back(class2Test); 
      else 
        closedClasses.push_back(class2Test);
  }
  List summaryMc = List::create(_["closedClasses"] = closedClasses,
                                _["transientClasses"] = transientClasses);
  return(summaryMc);
}
/*
#@ Tae: fully move in Rcpp
.summaryKernel<-function(object)
{
  matr<-object@transitionMatrix
  temp<-.commclassesKernel(matr)
  communicatingClassList<-.communicatingClasses(temp$C)
  transientStates<-names(which(temp$v==FALSE))
  closedClasses<-list()
  transientClasses<-list()
  for(i in 1:length(communicatingClassList))
  {
    class2Test<-communicatingClassList[[i]]
    if(length(intersect(class2Test,transientStates))>0) 
      transientClasses[[length(transientClasses)+1]]<-class2Test 
    else 
      closedClasses[[length(closedClasses)+1]]<-class2Test
  }
  summaryMc<-list(closedClasses=closedClasses, 
                  transientClasses=transientClasses)
  return(summaryMc)
}
*/

//here the kernel function to compute the first passage
// [[Rcpp::export(.firstpassageKernelRcpp)]]
NumericMatrix firstpassageKernel(NumericMatrix P, int i, int n){
  arma::mat G(P.begin(), P.nrow(), P.ncol());
  arma::mat Pa = G;
  arma::mat H(n, P.ncol()); //here Thoralf suggestion
  H.insert_rows(0, G.row(0)); //initializing the first row  
  arma::mat E = 1 - arma::eye(P.ncol(), P.ncol());
  
  for (int m = 1; m < n; m++) {
    G = Pa * (G*E);
    //H<-rbind(H,G[i,]) //removed thanks to Thoralf 
    H.insert_rows(m, G.row(i)); //here Thoralf suggestion
  }
  return wrap(H);
}
/*
#@ Tae: move in Rcpp
#here the kernel function to compute the first passage
.firstpassageKernel<-function(P,i,n){
  G<-P
  H <- matrix(NA, ncol=dim(P)[2], nrow=n) #here Thoralf suggestion
  H[1,]<-P[i,] #initializing the first row
  E<-1-diag(size(P)[2])
  for (m in 2:n) {
    G<-P%*%(G*E)
    #H<-rbind(H,G[i,]) #removed thanks to Thoralf 
    H[m,] <- G[i,] #here Thoralf suggestion
  }
  return(H)
}
*/

// greatest common denominator: to be moved in Rcpp
// [[Rcpp::export(.gcdRcpp)]]
double gcd (int f, int s) {
  int g, n, N, u;
  f = abs(f);
  s = abs(s);
    
  n = std::min(f,s);
  N = std::max(f,s);
  
  if (n==0) {
		g=N;
	}
	else {
		u=1;
		while (u!=0) {
			u=N%n;
			if (u==0) {
				g=n;
			}
			N=n;
			n=u;
		}
	}
	return g;
}
/* 
# greatest common denominator: to be moved in Rcpp
.gcd = function(f,s) {
  
  f <- abs(f)
  s <- abs(s)
  
  n <- min(f,s)
  N <- max(f,s)
  
	if (n==0) {
		g=N
	}
	else {
		u=1
		while (u!=0) {
			u=N%%n
			if (u==0) {
				g=n
			}
			N=n
			n=u
		}
	}
	return(g)
}
*/

/*** R
library(matlab)
i <- 1
size(i)
n <- size(i)[2]
n
matr <- ones(3, 3)
c <- colSums(matr)
c
d <- find(c)
d
find(d > 1)
a=i
a
a <- d
a
size(d)
n <- size(d)[2]
n
b <- zeros(1, 5)
b[1,d] <- ones(1, n)
b
matr <- ones(4, 4)
matr[a,]
as.numeric(matr[a,])
apply(t(matr) == t(matr), 2, sum)
apply(t(matr) == t(matr), 2, sum) == 3

#matr > 0
#(matr > 0) & (matr < 0)
#dim(matr)
#dim(matr)[1]
#sign(matr)
#1-diag(3)
#matr[2, 1] = 21
#matr[2, 2] = 22
#matr[[2]]

matr <- matrix(
c(TRUE,  TRUE, FALSE, FALSE, FALSE
,  TRUE,  TRUE, FALSE, FALSE, FALSE
, FALSE, FALSE,  TRUE,  TRUE, FALSE
, FALSE, FALSE,  TRUE,  TRUE, FALSE
, FALSE, FALSE, FALSE, FALSE,  TRUE)
, nrow = 5, ncol = 5, byrow =TRUE
)
dimnames(matr) <- list(c("a","b","c", "d", "e"), c("a","b","c", "d", "e"))
              
library(markovchain)
#.commStatesFinderRcpp(matr)
#.commclassesKernel(zeros(2))
matr <- matrix(
  c(0.0, 0.3333333, 0.0, 0.6666667, 0.0
  , 0.5, 0.0000000, 0.0, 0.0000000, 0.5
  , 0.0, 0.0000000, 0.5, 0.5000000, 0.0
  , 0.0, 0.0000000, 0.5, 0.5000000, 0.0
  , 0.0, 0.0000000, 0.0, 0.0000000, 1.0
  ), nrow = 5, byrow = TRUE,
  , dimnames = list(c("a","b","c", "d", "e"), c("a","b","c", "d", "e"))
)
.commclassesKernelRcpp(matr)
#.communicatingClassesRcpp(matr)

*/

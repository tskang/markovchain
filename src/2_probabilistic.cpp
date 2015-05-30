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

// [[Rcpp::export(.commclassesKernelRcpp)]]
extern "C" SEXP commclassesKernel(NumericMatrix P){
  unsigned int m = P.ncol(), n;
  CharacterVector stateNames = rownames(P);
//  NumericMatrix T(m);
//  NumericVector b;
  arma::vec a, b, c, d;
  arma::mat T = arma::zeros(m);
  unsigned int oldSum, newSum, i = 0;
  while(i <= m) {
    a = i;
//    b = NumericVector(m);
    b = arma::vec(m);
    b[i] = 1;
//    Rf_PrintValue(b);
    newSum = 0;
    oldSum = 1;
    while(oldSum != newSum) {
      oldSum = 0;
//      for(arma::vec::iterator it = b.begin(); it < b.end(); it ++)
//        if(*it > 0) oldSum += *it;
      oldSum = sum(find(b > 0));
      n = a.size();
      NumericVector r, temp; // = P(i, _);
//      arma::vec r, temp;
      for(unsigned int j = 0; j < n; j ++) {
        temp = P.row(a[j]);
        for(NumericVector::iterator it = temp.begin(); it != temp.end(); it ++)
          r.push_back(*it);
//        r = P(a[j], _);
//        r.insert(r.end(), P.row(a[j]));
      }
//      std::cout << r[0] << std::endl;
//      Rf_PrintValue(r);
      NumericMatrix matr(n, m, r.begin());
//      arma::mat matr(n, m, r.begin());
//      Rf_PrintValue(matr);
      c = sum(matr);
//      c = colSums(matr);
      arma::vec d = c;
			n = d.size();
      for(arma::vec::iterator it = d.begin(); it != d.end(); it++)
        b[*it] = 1;
//      b = arma::ones(1, n);
//      for(int j = 0; j < n; j++) 
//        b[d[j]] = 1;
      newSum = 0;
      for(arma::vec::iterator it = b.begin(); it != b.end(); it ++)
        if(*it > 0) newSum += *it;
//      for(int j = 0; j < b.size(); j ++)
//        if(b[j] > 0) newSum += b[j];
      newSum = sum(find(b > 0));
	    a = d;
    }
    T.insert_rows(i, b);
    i++;
  }

  arma::mat F = arma::trans(T);
  NumericMatrix C;
  arma::mat Ca(T.n_rows, T.n_cols);// = (T > 0);// & (F > 0);
  for(i = 0; i < T.n_rows; i ++)
    for(unsigned int j = 0; j < T.n_cols; j++)
      Ca(i, j) = (T(i, j) > 0 && F(i, j) > 0);
  LogicalVector v(T.n_cols);
  arma::mat tC = Ca.t();
  arma::mat tT = T.t();
  unsigned int sums[tC.n_cols];
//  arma::umat equalMat = (arma::trans(Ca) == arma::trans(T));
//  arma::vec colsums = colSums(equalMat);
  for(unsigned int j = 0; j < T.n_cols; j++) {
    sums[j] = 0;
    for(i = 0; i < T.n_rows; i ++)
      if(tC(i, j) == tT(i, j)) sums[j] ++;
    v[j] = (sums[j] == m);
  }
      
  C.attr("dimnames") = List::create(stateNames, stateNames);
//  C.names() = List::create(stateNames, stateNames);  
  v.names() = stateNames;
  return List::create(_["C"] = C, _["v"] = v);
}
//  m <- ncol(P)
//	stateNames <- rownames(P)
//	T <- zeros(m) 
//	i <- 1
//	while (i<=m) { 
//		a=i 
//		b<-zeros(1,m)
//		b[1,i]<-1
//		old<-1
//		new<-0
//		while (old != new) {
//			old <- sum(find(b>0))
//			n <- size(a)[2]
//			matr <- matrix(as.numeric(P[a,]),ncol=m,nrow=n) #fix
//			c <- colSums(matr)
//			d <- find(c)
//			n <- size(d)[2]
//			b[1,d] <- ones(1,n)
//			new <- sum(find(b>0))
//			a<-d
//		}
//		T[i,] <- b
//		i <- i+1 
//	}
//	F <- t(T)  
//	C <- (T>0)&(F>0)
//	v <- (apply(t(C)==t(T),2,sum)==m)
//	colnames(C) <- stateNames
//	rownames(C) <- stateNames
//	names(v) <- stateNames
//	out <- list(C=C,v=v)
//	return(out)


//returns the underlying communicating classes
// [[Rcpp::export(.communicatingClassesRcpp)]]
List communicatingClasses(LogicalMatrix adjMatr)
{
  int len = adjMatr.nrow();
  List classesList;
  for(int i = 0; i < len; i ++) {
    LogicalVector row2Check = adjMatr(i, _);
    CharacterVector rnames = row2Check.names();
    CharacterVector proposedCommClass;// = names(which(row2Check == true));
    for(int j = 0; j < row2Check.size(); j++) 
      if(row2Check[j] == true) 
          proposedCommClass.push_back(rnames(j));
    if (i > 0) {
      for(int j = 0; j < classesList.size(); j ++) {
        bool check = false;        
        CharacterVector cv = classesList[j];
        std::set<std::string> s1, s2;
        for(int k = 0; k < cv.size(); k ++) {
          s1.insert(as<std::string>(cv[k]));
          s2.insert(as<std::string>(proposedCommClass[k]));
        }
        check = std::equal(s1.begin(), s1.end(), s2.begin());
        if(check) {
          proposedCommClass = R_NilValue; break;
        }
      }
    }
    if(!Rf_isNull(proposedCommClass) && proposedCommClass.size() > 0) 
      classesList.push_back(proposedCommClass);    
  }
  return classesList;
}
/*
.communicatingClasses<-function(adjMatr)
{
  len <- dim(adjMatr)[1]
  classesList <- list()
  for(i in 1:len)
  {
    row2Check <- adjMatr[i,]
    proposedCommClass <- names(which(row2Check==TRUE))
    if(i>1) 
    {
      for(j in 1:(length(classesList)))
      {
        check <- FALSE
        check <- setequal(classesList[[j]],proposedCommClass)
        if(check==TRUE) {proposedCommClass<-NULL; break}
      }
    }
    if(!is.null(proposedCommClass)) classesList[[length(classesList)+1]]<-proposedCommClass
  }
  return(classesList)
*/

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
#ones(2)
#ones(1, 2)
#dim(matr)
#dim(matr)[1]
#sign(matr)
#1-diag(3)
#matr[2, 1] = 21
#matr[2, 2] = 22
#matr[[2]]

#.commStatesFinderRcpp(matr)
#.commclassesKernelRcpp(zeros(2))
*/

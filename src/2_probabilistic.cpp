// [[Rcpp::depends(RcppArmadillo)]]

//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export(.commclassesKernelRcpp)]]
List commclassesKernel(NumericMatrix P){
  int m = P.ncol(), n;
  String stateNames = rownames(P);
  NumericMatrix T(m), b;
  int a, d, old, newSum, i = 0;
  double oldSum, c;
  while(i <= m) {
    a = i;
    b = NumericMatrix(1, m);
    b(0, i) = 1;
//    Rf_PrintValue(b);
    old = 1;
    newSum = 0;
    oldSum = 0.0;
    while(old != newSum) {
      for(int j = 0; j < b.ncol(); j ++)
        if(b(0, j) > 0) oldSum += b(0, j);
      old = oldSum;
//      n = size(a)[2];
      n = 1;
      NumericVector r = P(a, _);
//      std::cout << r[0] << std::endl;
      Rf_PrintValue(r);
      NumericMatrix matr(m, n, r.begin());
//      Rf_PrintValue(matr);
        c = sum(matr.column(0));
        //  		d <- find(c)
        d = c;
//			n <- size(d)[2]
//			b[1,d] <- ones(1,n)
//        for(int j = 0; j < n; j++) 
//          b[0, j] = 1;
//			new <- sum(find(b>0))
        for(int j = 0; j < b.ncol(); j ++)
          if(b(0, j) > 0) newSum += b(0, j);
			  a = d;
    }
    i++;
  }
  List out;
  return out;
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
List communicatingClasses(NumericMatrix adjMatr)
{
  List classesList;
//  len <- dim(adjMatr)[1]
//  classesList <- list()
//  for(i in 1:len)
//  {
//    row2Check <- adjMatr[i,]
//    proposedCommClass <- names(which(row2Check==TRUE))
//    if(i>1) 
//    {
//      for(j in 1:(length(classesList)))
//      {
//        check <- FALSE
//        check <- setequal(classesList[[j]],proposedCommClass)
//        if(check==TRUE) {proposedCommClass<-NULL; break}
//      }
//    }
//    if(!is.null(proposedCommClass)) classesList[[length(classesList)+1]]<-proposedCommClass
//  }
//  return(classesList)
  return classesList;
}

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
size(d)
n <- size(d)[2]
n
dim(matr)
dim(matr)[1]
sign(matr)
1-diag(3)
matr[2, 1] = 21
matr[2, 2] = 22
matr[[2]]

#.commStatesFinderRcpp(matr)
#.commclassesKernelRcpp(zeros(2))
*/

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
  List out;
  return out;
}


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
arma::mat expm(arma::mat x) {
    arma::mat z(x.n_rows, x.n_cols);
    (*expmat)(x.begin(), x.n_rows, z.begin(), Ward_2);
    return z;
}
*/

//communicating states
// [[Rcpp::export(.commStatesFinderRcpp)]]
NumericMatrix commStatesFinder(NumericMatrix matr)
{
  //Reachability matrix
  double dimMatr = matr.nrow();
  int nrow = matr.nrow();
  int ncol = matr.ncol();

  arma::mat Z(nrow, ncol);
  arma::mat X(matr.begin(), nrow, ncol, false);
  arma::mat temp = arma::eye(dimMatr, ncol) + arma::sign(X);
  // raise to the power (dimMatr - 1)
  // temp = pow(temp, dimMatr - 1);
  
  for(int i = 0; i < nrow; i ++) 
    for(int j = 0; j < ncol; j ++) {
      
//      double temp = (eye(dimMatr) + Z(i, j)) % (dimMatr - 1);  
      arma::mat m;
    }
//  temp<-(eye(dimMatr)+Z)%^%(dimMatr-1)
//  double temp = (eye(dimMatr) + sign(matr)) % (dimMatr - 1);

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

#@ Tae: its upt to you to decide to move or to keep in R just calling the .commStatesFinder Rcpp version

is.accessible<-function(object, from, to)
{
  out<-FALSE
  statesNames<-states(object)
  fromPos<-which(statesNames==from)
  toPos<-which(statesNames==to)
  R<-.commStatesFinder(object@transitionMatrix)
  if(R[fromPos,toPos]==TRUE) out<-TRUE
  return(out)
}

#@ Tae: as above

#a markov chain is irreducible if is composed by only one communicating class
is.irreducible<-function(object)
{
  out<-FALSE
  tocheck<-.communicatingClasses(.commclassesKernel(object@transitionMatrix)$C)
  if(length(tocheck)==1) out<-TRUE
  return(out)
}

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
    if(length(intersect(class2Test,transientStates))>0) transientClasses[[length(transientClasses)+1]]<-class2Test else closedClasses[[length(closedClasses)+1]]<-class2Test
  }
  summaryMc<-list(closedClasses=closedClasses, 
                  transientClasses=transientClasses)
  return(summaryMc)
}

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

firstPassage<-function(object,state,n)
{
  P<-object@transitionMatrix
  stateNames<-states(object)
  i<-which(stateNames==state)
  outMatr<-.firstpassageKernel(P=P,i=i,n=n)
  colnames(outMatr)<-stateNames
  rownames(outMatr)<-1:n
  return(outMatr)
}
#periodicity


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

#@TAE: probably could be moved in Rcpp

#function to  get the period of a DTMC
period<-function(object) {
	check<-is.irreducible(object)
	if(check==FALSE){
		warning("The matrix is not irreducible")
		return(0)
	} else {
	P<-object@transitionMatrix
	n=size(P,2)
	v=zeros(1,n)
	v[1,1]=1
	w=numeric()
	d=0
	T=c(1)
	m=size(T,2)
	while (m>0 & d!=1) {
		i <- T[1]
		T <- T[-1]
		w <- c(w,i)
		j <- 1
		while (j<=n) {
			if (P[i,j]>0) {
				r=c(w,T)
				k=sum(r==j)
				if (k>0) {
					b=v[1,i]+1-v[1,j]
					d=.gcd(d,b)
				}
				else {
					T=c(T,j)
					v[1,j]=v[1,i]+1
				}
			}
			j=j+1
		}
		m=size(T,2)
	}
	v=v%%d
	return(d)
	}
}
*/

/*** R
library(matlab)
i <- 1
size(i)
n <- size(i)[2]
n
zeros(1, 3)
b<-zeros(1,4)
b[1,1]<-1
b[1, 3:4] <- ones(1, 2)
b
matr <- ones(3, 4)
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

.commStatesFinderRcpp(matr)
#.commclassesKernelRcpp(zeros(2))
*/

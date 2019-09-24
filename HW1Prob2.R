#Vivek Sathyanarayana
#CWID:  10442999
#FE 621- Fall 2019
#HW1-Problem 2

install.packages('timeDate',repos = "http://cran.us.r-project.org")
library('timeDate')
install.packages('zoo',repos = "http://cran.us.r-project.org")
library('zoo')

#Load HW1DATA2.RData
load('HW1DATA2.RData')

#Trim data to use only monthly expiries
AMZNopCh20190918 = AMZNopCh20190918[c('Oct.18.2019','Nov.15.2019','Dec.20.2019')]

#Compress Option chains to 10 strike prices using moneyness condition
#Create temp Spot-Strike column in chains (deleted at end of process)
Spot <- as.numeric(AMZN[1,4])
b=1
for (item in AMZNopCh20190918) {
  c=1
  for (callput in item) {
    Imp.Vol <- matrix(0,nrow = nrow(callput),ncol = 1,byrow = TRUE)
    BlSch.Op.Price <- matrix(0,nrow = nrow(callput),ncol = 1,byrow = TRUE)
    BiTree.Op.Price <- matrix(0,nrow = nrow(callput),ncol = 1,byrow = TRUE)
    SpotStrike <- matrix(0,nrow = nrow(callput),ncol = 1,byrow = TRUE)
    callput <- cbind(callput,Imp.Vol,BlSch.Op.Price,BiTree.Op.Price,SpotStrike)
    item[[c]] <- callput
    c=c+1
    }
  AMZNopCh20190918[[b]] = item
  b=b+1
}

#Compute Spot Strike difference to select 10 best strike rows
b=1
for (item in AMZNopCh20190918) {
  c=1
  for (callput in item) {
    callput$SpotStrike <- abs(Spot-as.numeric(callput$Strike))
    callput <- callput[order(callput$SpotStrike),]
    item[[c]] <- callput[1:10,-dim(callput)[2]]
    c=c+1
  }
  AMZNopCh20190918[[b]] <- item
  b=b+1
}

#Effective federal funds rate from http://www.federalreserve.gov/releases/H15/Current/
rate20190918 = 2.25/100

#Define Black Scholes Option pricing function
BSEuro <- function(K,S,r,sig,t,Type) {
  d1  <-  (log(S/K) + ((r + (sig^2)/2)*t)) / (sig*(t^0.5))
  d2  <-  d1 - (sig*(t^0.5))
  if(Type=='calls') {
    price =  S*pnorm(d1)-K*exp(-r*(t))*pnorm(d2)
  }
  if(Type=='puts') {
    price = S*(pnorm(d1)-1)+K*exp(-r*(t))*(1-pnorm(d2))
  }
  return(price)
}

#Define vega
vega <- function(K,S,r,t,sig) {
  d1 = (log(S/K) + (r+(sig^2)/2)*(t))/(sig*sqrt(t))
  v <- S*dnorm(d1)*sqrt(t)
  return(v)
}

#Define implied vol. function
impVol <- function(K,S,r,t,Type,m) {
  err = as.numeric(10^-3)
  x0 = as.numeric(0.4)
  delta = BSEuro(K,S,r,x0,t,Type) - m
  while (delta>err) {
    if(vega(K,S,r,t,x0)==0){
      x0=NA
      break
    }
    x0 = x0 - ((BSEuro(K,S,r,x0,t,Type)-m)/vega(K,S,r,t,x0))
    delta = abs(BSEuro(K,S,r,x0,t,Type) - m)
  }
  return(x0)
}

#Function to count number of days 't'
Nweekdays<- function(D1,D2){
  t = as.numeric(sum(!weekdays(seq(Date1, Date2, "days")) %in% c("Saturday", "Sunday")))
  return(t)
}

#Compute Implied Vol for Option chain
Date1 <- as.Date(index(AMZN),format = '%Y-%m-%d')
i=1
for (maturity in AMZNopCh20190918) {
  j=1
  Date2 <- as.Date(names(AMZNopCh20190918)[i],format = '%b.%d.%Y') 
  t = abs(Nweekdays(Date1,Date2))/360
  
  for (cp in maturity) {
    for (row in seq(1,dim(cp)[1],1)) {
      if ((is.na(cp$Bid[row])==TRUE)||(is.na(cp$Ask[row])==TRUE)){
        cp$Imp.Vol[row]=NA
        next
      }
      cp$Imp.Vol[row] = impVol(as.numeric(cp$Strike[row]),Spot,rate20190918,t,names(maturity)[j]
                                    ,(as.numeric(cp$Bid[row])+as.numeric(cp$Ask[row]))/2)
    }
    maturity[[j]] <-cp
    j=j+1
  }
  AMZNopCh20190918[[i]] <- maturity
  i=i+1
}

#Binomial tree European option pricing (as implemented in fOptions package) 
BinomialTreeOption = 
  function(TypeFlag = c("calls", "puts"), S, X, Time, r, b, sigma, n) {
    # Check Flags:
    TypeFlag = TypeFlag[1]
    if (TypeFlag == "calls") z = +1
    if (TypeFlag == "puts") z = -1    
    
    # Parameters:
    dt = Time / n
    u  = exp(sigma*sqrt(dt))
    d  = 1 / u
    p  = (exp(b*dt) - d) / (u - d)
    Df = exp(-r*dt)
    
    # Algorithm:
    OptionValue = z*(S*u^(0:n)*d^(n:0) - X)
    offset = 1
    Tree = OptionValue = (abs(OptionValue)+OptionValue)/2   
    
    # European Type:
    for (j in (n-1):0) {
      Tree <-c(Tree, rep(0, times=n-j))
      for (i in 0:j) {         
        OptionValue[i+offset] = (p*OptionValue[i+1+offset] + (1-p)*OptionValue[i+offset]) * Df 
        Tree = c(Tree, OptionValue[i+offset]) } } 
    Tree = matrix(rev(Tree), byrow = FALSE, ncol = n+1)
    # Return Value:
    invisible(Tree)
    }

#Create reporting table
report <- data.frame(matrix(ncol = 6))
colnames(report)<- c('Option','Maturity','Op.Type','Imp.Vol.','BS Price','BiTree Price')
op.type<- names(AMZNopCh20190918[[1]])
maturity <- names(AMZNopCh20190918)

#Compute option prices using Black Scholes & Binomial Tree
#Nested loop to navigate list of lists
b=1 #counter
for (item in AMZNopCh20190918) {
  c=1
  Date2 <- as.Date(names(AMZNopCh20190918[b]),format = '%b.%d.%Y') #Maturity date to compute 't'
  t = abs(Nweekdays(Date1,Date2))/360
  for (callput in item) {
    for (row in 1:dim(callput)[1]) {
      #Compute Option Price for each strike price
      sig <- AMZNopCh20190918[[b]][[c]]$Imp.Vol[row]
      if(is.na(sig)==TRUE)
        next
        #Compute black Scholes price
        callput$BlSch.Op.Price[row] = BSEuro(as.numeric(callput$Strike[row]),as.numeric(Spot),
                                             rate20190918,sig,t,names(item)[c])
        #Compute binomial tree price
        callput$BiTree.Op.Price[row] = BinomialTreeOption(b=0,r=rate20190918,S = Spot, X = as.numeric(callput$Strike[row]),
                                                          sigma = sig,Time = t,n=round(t*360),TypeFlag = names(item)[c])
        #Apped row to reporting table
        report <- rbind(report,c(rownames(callput)[row],names(AMZNopCh20190918)[b],names(item)[c],
                                 sig,callput$BlSch.Op.Price[row],callput$BiTree.Op.Price[row]))
    }
    item[[c]] <- callput  #Update list of lists with new dataframe
    c=c+1
    } 
  AMZNopCh20190918[[b]] = item #Update list of list with new list
  b=b+1
  }

#Delete initial NA row in report and display report
report <- report[-1,]
report
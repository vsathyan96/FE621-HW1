#Vivek Sathyanarayana
#CWID:  10442999
#FE 621- Fall 2019
#HW1-Problem 1

#Install packages 
install.packages('quantmod',repos = "http://cran.us.r-project.org")
library('quantmod')
install.packages('timeDate',repos = "http://cran.us.r-project.org")
library('timeDate')

#Load Data (Downloaded historical data since re-downloading would result in new data that have later expiry)
load('HW1Data.RData')

#Problem 1 part 1
#Function to download data from Yahoo
getData <- function(source="yahoo",symbol,startdate,enddate,option=FALSE,expiry=NULL) {
  if(option==FALSE)
    data = getSymbols(Symbols = symbol,src = source,from = startdate,to = enddate,auto.assign=FALSE, verbose=TRUE)
  if(option==TRUE)
    data = getOptionChain(Symbols = symbol,src = source,Exp = expiry)
return(data)
}

#Problem 1 part 2
#Extract data--  Code is commented out as data was downloaded and saved locally since it reflects option data as
#                 on 18.Sep.2019  and 20.Sep.2019

# AMZN = getData(source = "yahoo",symbol = "AMZN",startdate = '2019-09-18',enddate = '2019-09-21',option = FALSE)
# BA = getData(source = "yahoo",symbol = "BA",startdate = '2019-09-18',enddate = '2019-09-21',option = FALSE)
# SPY = getData(source = "yahoo",symbol = "SPY",startdate = '2019-09-18',enddate = '2019-09-21',option = FALSE)
#VIX = getData(source = "yahoo",symbol = "^VIX",startdate = '2019-09-18',enddate = '2019-09-21',option = FALSE)
#  
# AMZNopCh20190918 = getData(source = "yahoo",symbol = "AMZN",startdate = '2019-09-18', expiry = NULL, option = TRUE)
# BAopCh20190918 = getData(source = "yahoo",symbol = "BA",startdate = '2019-09-18', expiry = NULL, option = TRUE)
# 
# AMZNopCh20190920 = getData(source = "yahoo",symbol = "AMZN",startdate = '2019-09-20', expiry = NULL, option = TRUE)
# BAopCh20190920 = getData(source = "yahoo",symbol = "BA",startdate = '2019-09-20', expiry = NULL, option = TRUE)
# 
# SPY = getData(source = "yahoo",symbol = "SPY",startdate = '2019-09-18',enddate = '2019-09-20',option = FALSE)
# VIX = getData(source = "yahoo",symbol = "^VIX",startdate = '2019-09-18',enddate = '2019-09-20',option = FALSE)
# VXX = getData(source = "yahoo",symbol = "VXX",startdate = '2019-09-18',enddate = '2019-09-20',option = FALSE)
# 
# SPYopCh20190918 = getData(source = "yahoo",symbol = "SPY",startdate = '2019-09-18', expiry = NULL, option = TRUE)
# VIXopCh20190918 = getData(source = "yahoo",symbol = "^VIX",startdate = '2019-09-18', expiry = NULL, option = TRUE)
# VXXopCh20190918 = getData(source = "yahoo",symbol = "VXX",startdate = '2019-09-18', expiry = NULL, option = TRUE)
# 
# SPYopCh20190920 = getData(source = "yahoo",symbol = "SPY",startdate = '2019-09-20', expiry = NULL, option = TRUE)
# VIXopCh20190920 = getData(source = "yahoo",symbol = "^VIX",startdate = '2019-09-20', expiry = NULL, option = TRUE)
# VXXopCh20190920 = getData(source = "yahoo",symbol = "VXX",startdate = '2019-09-20', expiry = NULL, option = TRUE)

#Trim data to use only monthly expiries
AMZNopCh20190918 = AMZNopCh20190918[c('Oct.18.2019','Nov.15.2019','Dec.20.2019')]
AMZNopCh20190920 = AMZNopCh20190920[c('Oct.18.2019','Nov.15.2019','Dec.20.2019')]
BAopCh20190918 = BAopCh20190918[c('Oct.18.2019','Nov.15.2019','Dec.20.2019')]
BAopCh20190920 = BAopCh20190920[c('Oct.18.2019','Nov.15.2019','Dec.20.2019')]
SPYopCh20190918 = SPYopCh20190918[c('Oct.18.2019','Nov.15.2019','Dec.20.2019')]
SPYopCh20190920 = SPYopCh20190920[c('Oct.18.2019','Nov.15.2019','Dec.20.2019')]

DATA1 <- list(AMZNopCh20190918,BAopCh20190918,SPYopCh20190918)
names(DATA1) <- c('AMZNopCh20190918','BAopCh20190918','SPYopCh20190918')

#Note: BA does not have option data for December 20, 2019  maturity
DATA2 <- list(AMZNopCh20190920,BAopCh20190920,SPYopCh20190920)
names(DATA2) <- c('AMZNopCh20190920','BAopCh20190920','SPYopCh20190920')

#Problem 1 part 3
cat('SPY is an ETF run by State Street, which replicates the holdings of the S&P500 index. The VIX Index is a calculation 
designed to produce a measure of constant, 30-day expected volatility of the U.S. stock market,
    derived from real-time, mid-quote prices of SP 500 Index call and put options.
VXX (iPath S&P 500 VIX Short Term Futures ETN) is the largest and most liquid in the volatility ETF/ETN universe.')

#Problem 1 part 4
#Effective federal funds rate from http://www.federalreserve.gov/releases/H15/Current/
rate20190918 = 2.25/100
rate20190920 = 1.90/100

#Problem 1 part 5
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

#Problem 1 part 6
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

#Create implied vol columns in DATA1
a = 1
for (OpCh in DATA1) {
  b=1
  for (item in OpCh) {
    c=1
    for (callput in item) {
      IV <- matrix(0,nrow = nrow(callput),ncol = 1,byrow = TRUE)
      callput <- cbind(callput,IV)
      item[[c]] <- callput
      c=c+1
      }
    OpCh[[b]] = item
    b=b+1
    }
  DATA1[[a]] = OpCh
  a=a+1
}

#Function to compute number of days (t)
Nweekdays<- function(D1,D2){
  t = as.numeric(sum(!weekdays(seq(Date1, Date2, "days")) %in% c("Saturday", "Sunday")))
  return(t)
  }

#Compute Implied Vol for DATA1
#Create reporting table
report1 <- data.frame(matrix(ncol = 8))
colnames(report1)<- c('Equity','Option','Maturity','Op.Type','Moneyness','Imp.Vol.','Avg.Imp.Vol','VIX')

#Nested loops to navigate list of lists of option chains
a = 1  #counter
eq<- c('AMZN','BA','SPY')
moneyness<- c('In','Out')
op.type<- names(DATA1$AMZNopCh20190918[[1]])
maturity <- names(DATA1[[1]])

#Nested loop to navigate list of lists
for (OpCh in DATA1) {
  if (a==1) {
    Spot = AMZN[1,4]  #extract spot price for respective  equity
    Date1 <- as.Date(index(AMZN[1]),format = '%Y-%m-%d') #extract date of Spot to compute 't'
    }
  if (a==2) {
    Spot = BA[1,4] #extract spot price for respective  equity
    Date1 <- as.Date(index(BA[1]),format = '%Y-%m-%d') #extract date of Spot to compute 't'
    }
  if (a==3) {
    Spot = SPY[1,4] #extract spot price for respective  equity
    Date1 <- as.Date(index(SPY[1]),format = '%Y-%m-%d') #extract date of Spot to compute 't'
    }
  b=1 #counter
  for (item in OpCh) {
    c=1 #counter
    Date2 <- as.Date(names(OpCh[b]),format = '%b.%d.%Y')  #Extract date of option maturity to compute 't'
    t = abs(Nweekdays(Date1,Date2))/360
    
    for (callput in item) {
      for (row in 1:dim(callput)[1]) { #For loop to iterate through each row of strike price in option chain
        #Condition to print NA when required values do not exist
        #(or) print NA if "moneyness" condition fails
        if ((is.na(callput$Bid[row])==TRUE)||(is.na(callput$Ask[row])==TRUE)){
          callput$IV[row]=NA
          next
        }
        if (((as.numeric(Spot)/as.numeric(callput$Strike[row]))>1.05)||(as.numeric(Spot)/as.numeric(callput$Strike[row]))<0.95){
          callput$IV[row]=NA
          next
        }
        #Compute implied vol for each strike price
        callput$IV[row] = impVol(as.numeric(callput$Strike[row]),as.numeric(Spot),rate20190918,t,names(item)[c],
                                 ((as.numeric(callput$Bid[row]+callput$Ask[row])/2)))
        if (callput$Strike[row]<Spot)
          mon = moneyness[1] #Update moneyness column
        if (callput$Strike[row]>Spot)
          mon = moneyness[2] #Update moneyness column
        
        #Report for problem 1 part 7
        report1<-rbind(report1,c(eq[a],rownames(callput)[row],names(OpCh)[b],op.type[c],mon,callput$IV[row],0,0))
      }
      item[[c]] <- callput #Update list with updated dataframe
      c=c+1
      } 
    OpCh[[b]] = item #Update list of list with new list
    b=b+1
    }
  DATA1[[a]] = OpCh #Update list of list with new list
  a=a+1
}

#Problem 1 part 7
#Delete empty first row
report1<-report1[-1,]

#Compute average implied volitility
avgVol <- data.frame(matrix(ncol = 6))
colnames(avgVol) <- c("Equity","Maturity","Moneyness","Op.Type","Avg.Imp.Vol.","VIX")
for (item in eq) {
  i = 1 #counter
  for (t in maturity) {
    for (ty in op.type) {
      for (mo in moneyness) {
        avg = 0
        avgvol = report1[which(report1$Equity==item),]
        avgvol = report1[which(report1$Maturity==t),]
        avgvol = report1[which(report1$Op.Type==ty),]
        avgvol = as.numeric(report1[which(report1$Moneyness==mo),]$Imp.Vol.)
        avg = mean(avgvol)
        avgVol <- rbind(avgVol,c(item,t,mo,ty,avg,VIX[i,]$VIX.Close/100))
      }
    }
    i = i + 1
  }
}
#Delete empty first row and display average volitility
#Update report with average volitility
avgVol <- avgVol[-1,]
avgVol

#Update report table with avg.vol., VIX and display report
for (row in 1:dim(report1)[1]) {
  report1$Avg.Imp.Vol[row] <- avgVol[which((report1$Equity[row]==avgVol$Equity) & (report1$Maturity[row]==avgVol$Maturity)
                                          & (report1$Op.Type[row]==avgVol$Op.Type) &
                                            (report1$Moneyness[row]==avgVol$Moneyness)),]$Avg.Imp.Vol
  report1$VIX[row] <- avgVol[which((report1$Equity[row]==avgVol$Equity) & (report1$Maturity[row]==avgVol$Maturity)
                                          & (report1$Op.Type[row]==avgVol$Op.Type) &
                                            (report1$Moneyness[row]==avgVol$Moneyness)),]$VIX
  }
report1

#Comment on report observations (problem 1 part 7)
cat('AMZN implied vol ranges from 0.21 to 0.31; BA implied vol ranges from 0.25 to 0.33;
SPY implied vol ranges from 0.09 to 0.20. VIX value is ~0.15 which is much less than implied vol of AMZN and BA
and is almost the average of SPYs implied vol. From the report, it is evident that implied volatility increases
as time to maturity increases. When the option strike price increases,the option transitions from
in the money to out of the money and the implied volatility decreases')

#Problem 1 part 8
#Create empty computed option price column for DATA2
a = 1
for (OpCh in DATA2) {
  b = 1
  for (item in OpCh) {
    c = 1
    for (callput in item) {
      Comp.Op.Price <- matrix(NA,nrow = nrow(callput),ncol = 1,byrow = TRUE)
      callput <- cbind(callput,Comp.Op.Price)
      item[[c]] <- callput
      c = c+1
    }
    OpCh[[b]] = item
    b = b+1
  }
  DATA2[[a]] = OpCh
  a = a+1
}

#Compute Option Price using Black Scholes
a = 1
op.type<- names(DATA2$AMZNopCh20190920$Oct.20.2019)
maturity <- names(DATA2$AMZNopCh20190920)
for (OpCh in DATA2) {
  if (a==1) {
    Spot = AMZN[3,4]
    Date1 <- as.Date(index(AMZN[3]),format = '%Y-%m-%d')
  }
  if (a==2) {
    Spot = BA[3,4]
    Date1 <- as.Date(index(BA[3]),format = '%Y-%m-%d')
  }
  if (a==3) {
    Spot = SPY[3,4]
    Date1 <- as.Date(index(SPY[3]),format = '%Y-%m-%d')
  }
  b=1
  for (item in OpCh) {
    c=1
    Date2 <- as.Date(names(OpCh[b]),format = '%b.%d.%Y') 
    t = abs(Nweekdays(Date1,Date2))/360
    for (callput in item) {
      for (row in 1:dim(callput)[1]) {
        #Compute Option Price for each strike price
        sig <- DATA1[[a]][[b]][[c]]$IV[row]
        if(is.na(sig)==TRUE)
          next
        callput$Comp.Op.Price[row] = BSEuro(as.numeric(callput$Strike[row]),as.numeric(Spot),rate20190920,sig,t,names(item)[c])
      }
      item[[c]] <- na.omit(callput) #delete rows with NA
      c=c+1
    } 
    OpCh[[b]] = item
    b=b+1
  }
  DATA2[[a]] = OpCh
  a=a+1
}

#Display DATA2 with computed Black Scholes option price
DATA2

#HW1-Problem 2

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
report2 <- data.frame(matrix(ncol = 6))
colnames(report2)<- c('Option','Maturity','Op.Type','Imp.Vol.','BS Price','BiTree Price')
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
      report2 <- rbind(report2,c(rownames(callput)[row],names(AMZNopCh20190918)[b],names(item)[c],
                               sig,callput$BlSch.Op.Price[row],callput$BiTree.Op.Price[row]))
    }
    item[[c]] <- callput  #Update list of lists with new dataframe
    c=c+1
  } 
  AMZNopCh20190918[[b]] = item #Update list of list with new list
  b=b+1
}

#Delete initial NA row in report and display report
report2 <- report2[-1,]
report2
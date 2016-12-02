# Input parameters
Ne = 5								    #effective population size
u = 0.0025								#mutation rate
r = 0.01								#recombination rate
nGenerations = 400						#number of generations of random mating
interval = 50							#output frequency


# Setting up the matrices
n = 2*Ne								#there are 2*Ne gamets in the population
x = matrix(ncol=4,nrow=0)				#x stores all combinations of 4 haplotypes possible with n gametes 
for (i in 0:n){
	for (j in 0:(n-i)){
		for (k in 0:(n-i-j)){
			l = n-i-j-k
			x = rbind(x,c(i,j,k,l))
		}
	}
}
	
# Possibilities of mutations	
a = (1-u)*(1-u)							#no mutation and locus 1 abd no mutation at locus 1
b = (1-u)*u								#mutation at one locus
c = u*u								    #mutation at both locci
	
Mu = matrix(nrow=4,ncol=4,				#possibilities to get a haplotype given mutation
	                c(a,b,b,c,			#00 can come without mutations, from single mutation in 01, from single mutation in 10, from 2 mutations in 11
	                  b,a,c,b,			#01 can come from single mutation in 00, without mutations in 01 ,from 2 mutations in 10, from single mutation in 11 
	                  b,c,a,b,			#10 can come from single mutation in 00, from 2 mutations in 01, without mutations in 10 ,from single mutation in 11 
	                  c,b,b,a),byrow=T) #11 from 2 mutations in 00, from single mutation in 01, from single mutation in 10, without mutations 11
	                  
sizeX = nrow(x)							#number of combinations
A = matrix(nrow=sizeX,ncol=sizeX)
nMinus1Inv = 1/(n-1)
rsqr = c()    
for (i in 1:sizeX){ 
	   
    x00 = x[i,1]						#count for 00 haplotype from x matrix
    x01 = x[i,2]						#count for 01 haplotype from x matrix
    x10 = x[i,3]						#count for 10 haplotype from x matrix
    x11 = x[i,4]						#count for 11 haplotype from x matrix
    p00 = x00/n							#probability for 00 haplotype from x matrix
    p01 = x01/n							#probability for 01 haplotype from x matrix
    p10 = x10/n							#probability for 10 haplotype from x matrix
    p11 = x11/n							#probability for 11 haplotype from x matrix
    
    p1 = p10 + p11						#total probability of 1 at locus 1
    p2 = p01 + p11						#total probability of 1 at locus 2
    cov = p11 - p1*p2					#covariance between the locci
    rsqr = c(rsqr,cov^2/(p1*(1-p1)*p2*(1-p2))) #squared correlation between the locci
 
# Adding recombinations 
    p.recomb = c(p00*(1-r) + r*nMinus1Inv*( p00*(x00-1 + x10) + p01*(x00   + x10) ),	#00 can come from 00 without recombination, from recombination of 00 with 00 or 10, recombination of 01 with 00 or 10		
                 p01*(1-r) + r*nMinus1Inv*( p00*(x01   + x11) + p01*(x01-1 + x11) ),	#01 can come from 01 without recombination, from recombination of 00 with 01 or 11, recombination of 01 with 01 or 11
                 p10*(1-r) + r*nMinus1Inv*( p10*(x10-1 + x00) + p11*(x10   + x00) ),	#10 can come from 10 without recombination, from recombination of 10 with 10 or 00, recombination of 11 with 10 or 00
                 p11*(1-r) + r*nMinus1Inv*( p10*(x11   + x01) + p11*(x11-1 + x01) )	      #11 can come form 11 without recombination, from recombination of 10 with 11 or 01, recombination of 11 with 11 or 01
                )

#putting together mutations and recombinations				
    p.mut = Mu%*%p.recomb
    #print(x[i,])
    #print(p.mut)

#getting multinomial probabilities for all possible combinations of haplotypes
#transform matrx 286*286 -----> A%*%P to get new P in new generation, p is a vectot of 286 elements    
    for (j in 1:sizeX){
    	A[j,i] = dmultinom(x[j,],prob=p.mut,size=n)
    }
}

#Getting transition matrix probabilities, starting with equal probabilities for all haplotypes
p = c()
for (i in 1:sizeX){
	p = c(p,dmultinom(x[i,],prob=c(1,1,1,1),size=n))
}


count = interval
quartz(1)
par(mfrow=c(2,5))
expr2Vec = c();
for (i in 1:nGenerations){
		data = cbind(rsqr,p)
		data=data[complete.cases(data),]	#exclude the combinations where there is no segregation at both loci
		r2=round(data[,1],4)
		probs = data[,2]/sum(data[,2])	    #standardize probabilities to sum to one for segregating locci
		cat("generation",i,"\n")
		sumPr = by(probs,r2,sum)	#summing probabilities for given rsq
		res =cbind(sort(unique(r2)),as.matrix(sumPr)[,1])
		expr2=sum(as.matrix(res[,1]*res[,2])[,1]) #expected rsq is weighted mean of rsqs
		print(expr2)
		expr2Vec = cbind(expr2Vec,c(i,expr2))	#storing expected rsq for plotting
	if (count==interval){
		plot(sort(unique(r2)),sumPr,yaxp=c(0,1,10),ylim=c(0,1)) #plotting rsq for generations determined by interval
		count = 0 
	}
	count = count + 1
	p = A%*%p			 #recalculating new probabilities every generation										
}
cat("generation",nGenerations,"\n")
print(cbind(x,rsqr,p))		#output with haplotype combination,rsq and probability									

quartz(2)
windows(2)
plot(expr2Vec[1,],expr2Vec[2,]) #plotting expected rsq every generation  

          							
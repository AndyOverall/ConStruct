# *******************************
# R Script for ConStruct.Option.3
# By ADJ Overall
# University of Brighton
# Last Update October 2015
# *******************************
#
#
# To run this script, please copy ConStruct.Option.3.r to working directory
# containing your input files. Then type the following lines into R:
# >source("ConStruct.Option.3.r")
#	
# This option simulates a data file of diploid microsatellite-type genotypes
# for a specified degree of Fst between two subpopulations and proportion (c)
# of consanguineous individuals with parents related as 2r (e.g., offspring of
# 1st cousins r = 0.0625). 
#
# You will be prompted to input the population size (N), number of loci and 
# number of alleles for each locus. 
#
#
# Initialise variables
fst = 0
r = 0
c = 0
#
#
#
writeLines("Simulation of Dataset")
writeLines("Population Size (N):  ")
N = readLines(n = 1)
N = as.numeric(N)
#
#
#
writeLines("Number of Loci:  ")
num.loc = readLines(n = 1)
num.loc = as.numeric(num.loc)
#
#
#
writeLines("Required value of Fst:  ")
fst=readLines(n=1)
fst=as.numeric(fst)
#
#
# This is the actual degree of consanguinity simulated: 
writeLines("Required degree of Consanguinity between Parents (eg, 0.25; 0.125 etc):  ")
r = readLines(n = 1)
r = as.numeric(r)
#
#
# This is the proportion of the simulated population that are the offspring of consanguineous
# parents
writeLines("Required Proportion(%) of Population that is Consanguineous (eg, 0.1; 0.5 etc):  ")
c = readLines(n = 1)
c = as.numeric(c)
#
#
# 
writeLines("Analysis of Simulated of Dataset")  #  Simulated dataset is then analysed using ConStruct method
writeLines("Value of consanguinity being investigated (eg, 0.25; 0.125 etc): ")
rr = readLines(n = 1)
rr = as.numeric(rr)
#
#
#
max.alleles = 1000
f.resolution = 100
c.resolution = 100
resolution = 100
num.alleles = array(0, dim = c(max.alleles))
#
#
#  Specify  the number of alleles per locus:
writeLines("Type Number of Alleles for each Locus on each line: ")
for(i in 1:num.loc){
#  print("For Locus",quote=FALSE)
  print(i)
  num.alleles[i] = readLines(n = 1)
}
num.alleles = as.numeric(num.alleles)
#
#
# Specify the number of iterations for the algorithm that generates 
# two populations with allele frequencies at the specified Fst
writeLines("Type Number of iterations: ")
iteration = readLines(n = 1)
iteration = as.numeric(iteration)
#
#
# Generate random allele frequencies for each locus for two pops
# Initialise variables
limit = array(0, dim = c(num.loc))
x.old.1 = matrix(0, num.loc, max.alleles)   
x.old.2 = matrix(0, num.loc, max.alleles) 
x.new.1 = matrix(0, num.loc, max.alleles)   
x.new.2 = matrix(0, num.loc, max.alleles) 
allele.freq.1 = matrix(0, num.loc, max.alleles)
allele.freq.2 = matrix(0, num.loc, max.alleles)
homo.1.F = matrix(0, num.loc, max.alleles)
homo.1.O = matrix(0, num.loc, max.alleles)
homo.cumulative.1.F = matrix(0, num.loc, max.alleles)
homo.cumulative.1.O = matrix(0, num.loc, max.alleles)
het.1.F = array(0, dim = c(num.loc, max.alleles, max.alleles))
het.1.O = array(0, dim = c(num.loc, max.alleles, max.alleles))
het.cumulative.1.F = array(0, dim = c(num.loc, max.alleles, max.alleles))
het.cumulative.1.O = array(0, dim = c(num.loc, max.alleles, max.alleles))
homo.2.F = matrix(0, num.loc, max.alleles)
homo.2.O = matrix(0, num.loc, max.alleles)
homo.cumulative.2.F = matrix(0, num.loc, max.alleles)
homo.cumulative.2.O = matrix(0, num.loc, max.alleles)
het.2.F = array(0, dim = c(num.loc, max.alleles, max.alleles))
het.2.O = array(0, dim = c(num.loc, max.alleles, max.alleles))
het.cumulative.2.F = array(0, dim = c(num.loc, max.alleles, max.alleles))
het.cumulative.2.O = array(0, dim = c(num.loc, max.alleles, max.alleles))
genotype.1.allele.1.F = matrix(0, N, num.loc)
genotype.1.allele.2.F = matrix(0, N, num.loc)
genotype.2.allele.1.F = matrix(0, N, num.loc)
genotype.2.allele.2.F = matrix(0, N, num.loc)
genotype.1.allele.1.O = matrix(0, N, num.loc)
genotype.1.allele.2.O = matrix(0, N, num.loc)
genotype.2.allele.1.O = matrix(0, N, num.loc)
genotype.2.allele.2.O = matrix(0, N, num.loc)
all.count.1 = matrix(0, num.loc, max.alleles)
all.count.2 = matrix(0, num.loc, max.alleles)
all.freq = matrix(0, num.loc, max.alleles)
x.average = array(0, dim = c(max.alleles)) 
x.average.new = array(0, dim = c(max.alleles)) 
numerator = array(0, dim = c(max.alleles))
denominator = array(0, dim = c(max.alleles)) 
sum.allele = array(0, dim = c(num.loc))
fst.old = array(0, dim = c(num.loc))
fst.new = array(0, dim = c(num.loc))
likelihood.individual.c = array(0, dim = c(N))
likelihood.individual.s = array(0, dim = c(N))
ln.likelihood.individual = array(0, dim = c(N))
ln.likelihood.cs = matrix(0, resolution, resolution)
e.likelihood.cs = matrix(0, resolution, resolution)
num = 0
#
#
# Initial allele frequencies are equal
for(i in 1:num.loc){
  for(j in 1:num.alleles[i]){
    x.old.1[i, j] = 1 / num.alleles[i]
    x.old.2[i, j] = 1 / num.alleles[i]
  }
}
#
# 
# Specifies maximum number of alleles	
for(i in 1:num.loc){
  limit[i] = num.alleles[i]
#
  while(sum(x.old.1[i, ]) >= 1){
    for(j in 1:limit[i]){
      x.old.1[i, j] = runif(1, 0, 1)  #  Generate random allele frequencies until they sum to 1
    }
  }
  x.old.1[i, limit[i]] = 1 - sum(x.old.1[i, 1:limit[i] - 1])
#                                                                                                                                                                                                                                                                                                                      
  while(sum(x.old.2[i, ]) >= 1){
    for(j in 1:limit[i]){
      x.old.2[i, j] = runif(1, 0, 1)
    }
  }	
  x.old.2[i, limit[i]] = 1 - sum(x.old.2[i, 1:limit[i] - 1])	
#
#
# Calculates initial Fst between two simulated populations
  for(j in 1:num.alleles[i]){
    x.average[j] = (x.old.1[i, j] + x.old.2[i, j]) / 2
    numerator[j] = (x.old.1[i, j] - x.average[j]) ^ 2
    denominator[j] = x.average[j] * (1 - x.average[j])
  }
  fst.old[i] = sum(numerator) / sum(denominator)
  for(k in 1:iteration){            #  Number of iterations to find desired Fst
    for(j in 1:(limit[i]-1)){       #  For each allele 
      sign=runif(1, 0, 1)           #  Random number between 0 - 1
      if(sign>0.5) x=1 else x=-1   
      if(x.old.1[i,j]+(x*0.01)>0 & x.old.1[i,j]+(x*0.01)<1) x.new.1[i,j]=x.old.1[i,j]+(x*0.01)
      sign=runif(1, 0, 1)
      if(sign>0.5) x=1 else x=-1
      if(x.old.2[i,j]+(x*0.01)>0 & x.old.2[i,j]+(x*0.01)<1) x.new.2[i,j]=x.old.2[i,j]+(x*0.01)		
    }
    if(sum(x.new.1[i,1:(limit[i]-1)])<=1){
      x.new.1[i,limit[i]]=1-sum(x.new.1[i,1:limit[i]-1])
    }
    if(sum(x.new.2[i,1:(limit[i]-1)])<=1){
      x.new.2[i,limit[i]]=1-sum(x.new.2[i,1:limit[i]-1])
    }
    if(sum(x.new.1[i,])==1 & sum(x.new.2[i,])==1){
      for(j in 1:limit[i]){
        x.average.new[j]=(x.new.1[i,j]+x.new.2[i,j])/2
	numerator[j]=(x.new.1[i,j]-x.average.new[j])^2
        denominator[j]=x.average.new[j]*(1-x.average.new[j])
      }
      fst.new[i]=sum(numerator)/sum(denominator)    #fst across all alleles
      if((fst-fst.new[i])^2<(fst-fst.old[i])^2){
        for(j in 1:limit[i]){
          x.old.1[i,j]=x.new.1[i,j]  
	  x.old.2[i,j]=x.new.2[i,j]
          fst.old[i]=fst.new[i]
        }
      }
    }	
  }
  for(j in 1:limit[i]){
    allele.freq.1[i,j]=x.old.1[i,j]
    allele.freq.2[i,j]=x.old.2[i,j]
  }
}
print("Simulated Fst values for each locus")
print(data.frame(fst.new))
#
#
N.O = round((1 - c) * N / 2)   #  Proportion of simulated population non-inbred
N.F = round(c * N / 2)         #  Proportion of simulated population inbred
#
#
#
# Generate population number 1
# Generate homozygotes in inbred population(homo.F) and non-inbred(homo.O)
for(i in 1:num.loc){
  for(j in 1:limit[i]){
    homo.1.F[i, j] = allele.freq.1[i, j] * (r + (1 - r) * allele.freq.1[i, j])
    homo.1.O[i, j] = (allele.freq.1[i, j]) ^ 2
  }
# Generate heterozygous genotypes
  for(j in 1:(limit[i] - 1)){
    for(k in (j + 1):limit[i]){
      het.1.F[i, j, k] = 2 * allele.freq.1[i, j] * allele.freq.1[i, k] * (1 - r)
      het.1.O[i, j, k] = 2 * allele.freq.1[i, j] * allele.freq.1[i, k]
    }
  }
# Generate cumulative sum
  for(j in 1:limit[i]){
    homo.cumulative.1.F[i, j] = sum(homo.1.F[i, 1:j])
    homo.cumulative.1.O[i, j] = sum(homo.1.O[i, 1:j])
  }
  het.cumulative.1.F.temp = homo.cumulative.1.F[i, limit[i]]
  het.cumulative.1.O.temp = homo.cumulative.1.O[i, limit[i]]
  for(j in 1:(limit[i] - 1)){
    for(k in (j + 1):limit[i]){
      het.cumulative.1.F[i, j, k] = het.1.F[i, j, k] + het.cumulative.1.F.temp
      het.cumulative.1.O[i, j, k] = het.1.O[i, j, k] + het.cumulative.1.O.temp
      het.cumulative.1.F.temp = het.cumulative.1.F[i, j, k]
      het.cumulative.1.O.temp = het.cumulative.1.O[i, j, k]
    }
  }
}
if(N.F > 0){      
  for(i in 1:N.F){            #inbred population		
    for(j in 1:num.loc){
      random.number = runif(1, 0, 1)
      check = FALSE
      for(k in 1:limit[j]){				
        if(random.number <= homo.cumulative.1.F[j, k] & check == FALSE){
          genotype.1.allele.1.F[i, j] = k        
          genotype.1.allele.2.F[i, j] = k
          all.count.1[j, k] = all.count.1[j, k] + 2
          check = TRUE
        }	
      }
      for(l in 1:(limit[j] - 1)){
        for(m in (l + 1):limit[j]){
          if(random.number <= het.cumulative.1.F[j, l, m] & check == FALSE){
            genotype.1.allele.1.F[i, j] = l
            genotype.1.allele.2.F[i, j] = m
            all.count.1[j, l] = all.count.1[j, l] + 1
            all.count.1[j, m] = all.count.1[j, m] + 1
            check = TRUE
          }
        }
      }
    }
  }
}
if(N.O > 0){
  for(i in (2 * N.F + 1):(2 * N.F + N.O)){		#non-inbred population
    for(j in 1:num.loc){
      random.number = runif(1, 0, 1)
      check = FALSE
      for(k in 1:limit[j]){				
        if(random.number <= homo.cumulative.1.O[j, k] & check == FALSE){
          genotype.1.allele.1.O[i, j] = k        
          genotype.1.allele.2.O[i, j] = k
          all.count.1[j, k] = all.count.1[j, k] + 2
          check = TRUE
        }	
      }    
      for(l in 1:(limit[j] - 1)){
        for(m in (l + 1):limit[j]){
          if(random.number <= het.cumulative.1.O[j, l, m] & check == FALSE){
            genotype.1.allele.1.O[i, j] = l
            genotype.1.allele.2.O[i, j] = m
            all.count.1[j, l] = all.count.1[j, l] + 1
            all.count.1[j, m] = all.count.1[j, m] + 1
            check = TRUE
          }
        }
      }
    }
  }
}
#
# Generate population number 2
# Generate homozygotes in inbred population(homo.F) and non-inbred(homo.O)
for(i in 1:num.loc){
  for(j in 1:limit[i]){
    homo.2.F[i, j] = allele.freq.2[i, j] * (r + (1 - r) * allele.freq.2[i, j])
    homo.2.O[i, j] = (allele.freq.2[i, j]) ^ 2
  }
# Generate heterozygous genotypes
  for(j in 1:(limit[i] - 1)){
    for(k in (j + 1):limit[i]){
      het.2.F[i, j, k] = 2 * allele.freq.2[i, j] * allele.freq.2[i, k] * (1 - r)
      het.2.O[i, j, k] = 2 * allele.freq.2[i, j] * allele.freq.2[i, k]
    }
  }
# Generate cumulative sum
  for(j in 1:limit[i]){
    homo.cumulative.2.F[i, j] = sum(homo.2.F[i, 1:j])
    homo.cumulative.2.O[i, j] = sum(homo.2.O[i, 1:j])
  }
  het.cumulative.2.F.temp = homo.cumulative.2.F[i, limit[i]]
  het.cumulative.2.O.temp = homo.cumulative.2.O[i, limit[i]]
  for(j in 1:(limit[i] - 1)){
    for(k in (j + 1):limit[i]){
      het.cumulative.2.F[i, j, k] = het.2.F[i, j, k] + het.cumulative.2.F.temp
      het.cumulative.2.O[i, j, k] = het.2.O[i, j, k] + het.cumulative.2.O.temp
      het.cumulative.2.F.temp = het.cumulative.2.F[i, j, k]
      het.cumulative.2.O.temp = het.cumulative.2.O[i, j, k]
    }
  }
}
if(N.F > 0){     
  for(i in (N.F + 1):(2 * N.F)){            #inbred population		
    for(j in 1:num.loc){
      random.number = runif(1, 0, 1)
      check = FALSE
      for(k in 1:limit[j]){				
        if(random.number <= homo.cumulative.2.F[j, k] & check == FALSE){
          genotype.2.allele.1.F[i, j] = k        
          genotype.2.allele.2.F[i, j] = k
          all.count.2[j, k] = all.count.2[j, k] + 2
          check = TRUE
        }	
      }  
      for(l in 1:(limit[j] - 1)){
        for(m in (l + 1):limit[j]){
          if(random.number <= het.cumulative.2.F[j, l, m] & check == FALSE){
            genotype.2.allele.1.F[i, j] = l
            genotype.2.allele.2.F[i, j] = m
            all.count.2[j, l] = all.count.2[j, l] + 1
            all.count.2[j, m] = all.count.2[j, m] + 1
            check = TRUE
          }
        }
      }
    }
  }
}
if(N.O > 0){
  for(i in (2 * N.F + N.O + 1):(2 * N.F + 2 * N.O)){		#non-inbred population
    for(j in 1:num.loc){
    random.number = runif(1, 0, 1)
      check = FALSE
      for(k in 1:limit[j]){				
        if(random.number <= homo.cumulative.2.O[j, k] & check == FALSE){
          genotype.2.allele.1.O[i, j] = k        
          genotype.2.allele.2.O[i, j] = k
          all.count.2[j, k] = all.count.2[j, k] + 2
          check = TRUE
        }	
      } 
      for(l in 1:(limit[j] - 1)){
        for(m in (l + 1):limit[j]){
          if(random.number <= het.cumulative.2.O[j, l, m] & check == FALSE){
            genotype.2.allele.1.O[i, j] = l
            genotype.2.allele.2.O[i, j] = m
            all.count.2[j, l] = all.count.2[j, l] + 1
            all.count.2[j, m] = all.count.2[j, m] + 1
            check = TRUE
          }
        }
      }
    }
  }
} 
#
#
# Calculate allele frequencies
#
for(i in 1:num.loc){
  for(j in 1:max.alleles){
    all.freq[i, j] = (all.count.1[i, j] + all.count.2[i, j]) / (2 * N)
  }
}
#
# Calculate maximum Fst
homozygous.frequency = 0
for(i in 1:num.loc){
  for(j in 1:max.alleles){
    if(all.freq[i,j] > 0){
      homozygous.frequency = homozygous.frequency + (all.freq[i, j]) ^ 2
    }
  }
}
homozygous.frequency = homozygous.frequency / num.loc
fst.max = homozygous.frequency / (1 - homozygous.frequency) 
#
if(fst.max > 1){
  fst.max = 1
}
writeLines("Maximum value of Fst = ")
print(fst.max)  
#
# ANALYSIS
#
# Calculate probability of multilocus genotypes for range of f values
#
substructure.homozygous = function(f,all.freq.1){
  all.freq.1 * (f + (1 - f) * all.freq.1)
}
substructure.heterozygous = function(f, all.freq.1, all.freq.2){
  2 * all.freq.1 * all.freq.2 * (1 - f)
}
consanguinity.homozygous = function(f, rr, all.freq.1){
  all.freq.1 * (rr + (1 - rr) * (f + (1 - f) * all.freq.1))
}
consanguinity.heterozygous = function(f, rr, all.freq.1, all.freq.2){
  2 * all.freq.1 * all.freq.2 * (1 - rr) * (1 - f)
}
res = 1
ln.likelihood = array(0, resolution)
likelihood.locus.res = matrix(0, resolution, N, num.loc)
#
# Iterate through values of f and c to estimate maximum likelihood values given simulated genotypes
f.resolution = 1
f = 0
repeat{
  c.resolution = 1
  c = 0
  repeat{
    likelihood.locus.s = matrix(0, N, num.loc)
    likelihood.locus.c = matrix(0, N, num.loc)
    if(N.F > 0){
      for(i in 1:N.F){                    #inbred population number 1, if homozygous 
        for (j in 1:num.loc){
          if(genotype.1.allele.1.F[i, j] > 0){
            if(genotype.1.allele.1.F[i, j] == genotype.1.allele.2.F[i, j]){
              likelihood.locus.s[i, j] = substructure.homozygous(f, all.freq[j, genotype.1.allele.1.F[i, j]])
              likelihood.locus.c[i, j] = consanguinity.homozygous(f, rr, all.freq[j, genotype.1.allele.1.F[i, j]])
            }
          }
        }
      }
      for(i in (N.F + 1):(2 * N.F)){                    #inbred population number 2, if homozygous
        for (j in 1:num.loc){
          if(genotype.2.allele.1.F[i, j] > 0){
            if(genotype.2.allele.1.F[i, j] == genotype.2.allele.2.F[i, j]){
              likelihood.locus.s[i, j] = substructure.homozygous(f, all.freq[j, genotype.2.allele.1.F[i, j]])
              likelihood.locus.c[i, j] = consanguinity.homozygous(f, rr, all.freq[j, genotype.2.allele.1.F[i, j]])
            }
          }
        }
      }
      for(i in 1:N.F){                    #inbred population number 1, if heterozygous
        for (j in 1:num.loc){
          if(genotype.1.allele.1.F[i, j] > 0){
            if(genotype.1.allele.1.F[i, j] != genotype.1.allele.2.F[i, j]){
              likelihood.locus.s[i, j] = substructure.heterozygous(f, all.freq[j, genotype.1.allele.1.F[i, j]], all.freq[j, genotype.1.allele.2.F[i, j]])
              likelihood.locus.c[i, j] = consanguinity.heterozygous(f, rr, all.freq[j, genotype.1.allele.1.F[i, j]], all.freq[j, genotype.1.allele.2.F[i, j]])
            }
          }
        }
      }
      for(i in (N.F + 1):(2 * N.F)){                    #inbred population number 2, if heterozygous
        for (j in 1:num.loc){
          if(genotype.2.allele.1.F[i, j] > 0){
            if(genotype.2.allele.1.F[i, j] != genotype.2.allele.2.F[i, j]){
              likelihood.locus.s[i, j] = substructure.heterozygous(f, all.freq[j, genotype.2.allele.1.F[i, j]], all.freq[j, genotype.2.allele.2.F[i, j]])
              likelihood.locus.c[i, j] = consanguinity.heterozygous(f, rr,all.freq[j, genotype.2.allele.1.F[i, j]],all.freq[j, genotype.2.allele.2.F[i, j]])
            }
          }
        }
      }
    }
    if(N.O > 0){
      for(i in (2 * N.F + 1):(2 * N.F + N.O)){                    #non-inbred population number 1, if homozygous
        for (j in 1:num.loc){
          if(genotype.1.allele.1.O[i, j] > 0){
            if(genotype.1.allele.1.O[i, j] == genotype.1.allele.2.O[i, j]){
              likelihood.locus.s[i, j] = substructure.homozygous(f, all.freq[j, genotype.1.allele.1.O[i, j]])
              likelihood.locus.c[i, j] = consanguinity.homozygous(f, rr, all.freq[j, genotype.1.allele.1.O[i, j]])
            }
          }
        }
      }
      for(i in (2 * N.F + N.O + 1):(2 * N.F + 2 * N.O)){                    #non-inbred population number 2, if homozygous
        for (j in 1:num.loc){
          if(genotype.2.allele.1.O[i, j] > 0){
            if(genotype.2.allele.1.O[i, j] == genotype.2.allele.2.O[i, j]){
              likelihood.locus.s[i, j] = substructure.homozygous(f, all.freq[j, genotype.2.allele.1.O[i, j]])
              likelihood.locus.c[i, j] = consanguinity.homozygous(f, rr, all.freq[j, genotype.2.allele.1.O[i, j]])
            }
          }
        }
      }
      for(i in (2 * N.F + 1):(2 * N.F + N.O)){                    #non-inbred population number 1, if heterozygous
        for (j in 1:num.loc){
          if(genotype.1.allele.1.O[i, j] > 0){
            if(genotype.1.allele.1.O[i, j] != genotype.1.allele.2.O[i, j]){
              likelihood.locus.s[i, j] = substructure.heterozygous(f, all.freq[j, genotype.1.allele.1.O[i, j]], all.freq[j, genotype.1.allele.2.O[i, j]])
              likelihood.locus.c[i, j] = consanguinity.heterozygous(f, rr, all.freq[j, genotype.1.allele.1.O[i, j]], all.freq[j, genotype.1.allele.2.O[i, j]])
            }
          }
        }
      }
      for(i in (2 * N.F + N.O + 1):(2 * N.F + 2 * N.O)){                    #non-inbred population number 2, if heterozygous 
        for (j in 1:num.loc){
          if(genotype.2.allele.1.O[i, j] > 0){
            if(genotype.2.allele.1.O[i, j] != genotype.2.allele.2.O[i, j]){
              likelihood.locus.s[i, j] = substructure.heterozygous(f, all.freq[j, genotype.2.allele.1.O[i, j]], all.freq[j, genotype.2.allele.2.O[i, j]])
              likelihood.locus.c[i, j] = consanguinity.heterozygous(f, rr, all.freq[j, genotype.2.allele.1.O[i, j]], all.freq[j, genotype.2.allele.2.O[i, j]])
            }
          }
        }
      }
    }
    if(N.F > 0){
      for(i in 1:N.F){                     #inbred population number 1, if missing genotype
        for (j in 1:num.loc){
          if(genotype.1.allele.1.F[i, j] == 0){
            likelihood.locus.s[i, j] = 1
            likelihood.locus.c[i, j] = 1
          }
        }
      }
      for(i in (N.F + 1):(2 * N.F)){                     #inbred population number 2, if missing genotype
        for (j in 1:num.loc){
          if(genotype.2.allele.1.F[i, j] == 0){
            likelihood.locus.s[i, j] = 1
            likelihood.locus.c[i, j] = 1
          }
        }
      }
    }
    if(N.O > 0){
      for(i in (2 * N.F + 1):(2 * N.F + N.O)){                     #non-inbred population number 1, if missing genotype
        for (j in 1:num.loc){
          if(genotype.1.allele.1.O[i, j] == 0){
            likelihood.locus.s[i, j] = 1
            likelihood.locus.c[i, j] = 1
          }
        }
      }
      for(i in (2 * N.F + N.O + 1):(2 * N.F + 2 * N.O)){                     #non-inbred population number 2, if missing genotype
        for (j in 1:num.loc){
          if(genotype.2.allele.1.O[i, j] == 0){
            likelihood.locus.s[i, j] = 1
            likelihood.locus.c[i, j] = 1
          }
        }
      }
    }
    for(i in 1:N){
      if(prod(likelihood.locus.s[i, ]>0)){
        likelihood.individual.s[i] = prod(likelihood.locus.s[i, ])
      } else {
        likelihood.individual.s[i] = 10 ^ (-100)
      }
      if(prod(likelihood.locus.c[i, ]>0)){ 
        likelihood.individual.c[i] = prod(likelihood.locus.c[i, ])
      } else{
        likelihood.individual.c[i] = 10 ^ (-100)
      }
      ln.likelihood.individual[i] = log(c * likelihood.individual.c[i] + (1 - c) * likelihood.individual.s[i])
    }
#
    ln.likelihood.cs[f.resolution, c.resolution] = sum(ln.likelihood.individual)
    c.resolution = c.resolution + 1
    c = c + 0.01
    if(c.resolution == 101){
      break
    }
  }  
  f.resolution = f.resolution+1
  f.values = fst.max / resolution
  f = f + f.values
  if(f.resolution == 101){
    break
  }
}
max.value = max(ln.likelihood.cs)
temp = matrix(0, resolution, resolution)
for(i in 1:resolution){
  for(j in 1:resolution){      
    temp[i, j] = ln.likelihood.cs[i,j] - max.value
  }
}
for(i in 1:resolution){
  for(j in 1:resolution){
    ln.likelihood.cs[i, j] = temp[i, j]
  }
}
for(i in 1:resolution){
  for(j in 1:resolution){
    if(ln.likelihood.cs[i, j] == 0){
      ML.F = ((i - 1) / resolution) * fst.max
      ML.C = (j / resolution) - 0.01
    }
  }
}
writeLines("Maximum Likelihood value: ")
writeLines("Fst=")
print(ML.F)
writeLines("C = ")
print(ML.C)
#
#
#
for(i in 1:resolution){
  for(j in 1:resolution){
    e.likelihood.cs[i,j]=exp(ln.likelihood.cs[i,j])
  }
}
total=sum(e.likelihood.cs)
for(i in 1: resolution){
  for(j in 1: resolution){
    e.likelihood.cs[i, j] = e.likelihood.cs[i, j] / total
  }
}
write.table(e.likelihood.cs,file = "ConStruct.Sim.Outfile.txt")
#
x = seq(0, fst.max - f.values, by = f.values)
y = seq(0, 0.99, by = 0.01)
contour(x, y, e.likelihood.cs,xlab = "FST",ylab = "% Consanguinity")
writeLines("G = ")
G = exp(log(max(e.likelihood.cs))-2)print(G)

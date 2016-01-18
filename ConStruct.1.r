# **********************************************************************
# R Script for ConStruct.1
# By ADJ Overall
# University of Brighton
# Last Update November 2015
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# **********************************************************************
#
#
# To run this script, copy into the working directory containing your input files and type 
# > source("ConStruct.1.r"). 
# There are three functions:
#
#  (1) max.likelihood = function(data, max.alleles, resolution)
  #  Arguments:
  #    data is the input file
  #    max.alleles places an uppermost limit on the number of alleles considered
  #    resolution is the resolution of the F parameter	
  #  Example of use:
  #  > max.likelihood(data=”infile.txt", max.alleles=1000, resolution=100)
# 
#  (2) construct = function(data, max.alleles, f.resolution, c.resolution, r)
  #  Arguments:
  #    data is the input file
  #    max.alleles places an uppermost limit on the number of alleles considered
  #    f.resolution is the resolution of the Fst parameter
  #    c.resolution is the resolution on the c parameter
  #    r is the value of the inbreeding coefficient being considered for the      
  #      analysis of the dataset 
  #  Example of use:
  #  > construct(data=”infile.txt", max.alleles=1000, f.resolution=100, c.resolution=100, r=0.0625)
#
#  (3) simulate = function(N, num.loc, fst, r.actual, c, r.consider, max.alleles, f.resolution, c.resolution, iteration)
  # Arguments:
  #   N is the total sample size
  #   num.loc is the number of loci
  #   fst is the value of Fst that is to be simulated between two populations
  #   r.actual is the inbreeding coefficient of the inbred individuals
  #   c is the proportion of the population inbred to degree r.actual
  #   r.consider is the value of the inbreeding coefficient being considered for the      
  #     analysis of the simulated dataset 
  #   max.alleles places an uppermost limit on the number of alleles considered
  #   f.resolution is the resolution of the Fst parameter
  #   c.resolution is the resolution on the c parameter
  #   iteration is the number of iterations of the simulation run through in order
  #     to arrive at the specified, simulated Fst
  #  Example of use:
  #  > simulate(N=200, num.loc=12, fst=0.05, r.actual=0.05, c=0.5, r.consider=0.05, max.alleles=100, f.resolution=100, c.resolution=100, iteration=10000)
  #  It is necessary to specify number of alleles per locus (e.g., for 4 loci):
  #  num.alleles = c(8,8,8,8)
#
#  *****************************************************************
#  MAX.LIKELIHOOD
#  *****************************************************************
max.likelihood = function(data, max.alleles, resolution){
# Input file defined and NAs converted to 0
infile = read.table(data)
infile [is.na(infile)] <- 0
# Population size is defined by the number of lines
N = nrow(infile)
# Number of loci is defined by number of columns
num.loc = ncol(infile) / 2
# Initialise variables
all.count = matrix(0, num.loc, max.alleles)
all.freq = matrix(0, num.loc, max.alleles)
homozygous.frequency = matrix(0, num.loc, max.alleles)
count.alleles = matrix(0, max.alleles)
likelihood.individual = matrix(0, N)
ln.likelihood = matrix(0, resolution)
e.likelihood = matrix(0, resolution)
hom.total = matrix(0, num.loc)
#
#set up input data
#
# Read diploid genotypes as sets of two alleles
allele.1 = matrix(0, N, (2 * num.loc))
allele.2 = matrix(0, N, (2 * num.loc))
for(i in 1:N){
  for(j in seq(1, by = 2, (2 * num.loc - 1))){
    allele.1[i,j] = infile[i,j]
  }
}
for(i in 1:N){
  for(j in seq(2, by = 2, (2 * num.loc))){
    allele.2[i,j] = infile[i,j]
  }
}
allele.1.seq = seq(1, 2 * num.loc, by = 2)
allele.2.seq = seq(2, 2 * num.loc, by = 2)
for(i in 1:N){
  for (j in 1:num.loc){
    allele.1[i, j] = allele.1[i, allele.1.seq[j]]
  }
}
for(i in 1:N){
  for (j in 1:num.loc){
    allele.2[i, j] = allele.2[i, allele.2.seq[j]]
  }
}

# Calculate allele frequencies
for(i in 1:N){
  for (j in 1:num.loc){
    if (allele.1[i, j] >0){
      all.count[j, allele.1[i, j]] = all.count[j, allele.1[i, j]] + 1
    }
  }
}
for(i in 1:N){
  for (j in 1:num.loc){
    if (allele.2[i, j] > 0){
      all.count[j, allele.2[i, j]] = all.count[j, allele.2[i, j]] + 1
    }
  }
}
#
locus.n = matrix(0, num.loc)
for(i in 1:num.loc){
  for(j in 1:max.alleles){
    locus.n[i] = sum(all.count[i, ])
  }
}
for(i in 1:num.loc){
  for(j in 1:max.alleles){
    all.freq[i, j]=all.count[i, j] / locus.n[i]
  }
}
#
# Calculate maximum Fst using Hedrick’s equation (1-Hs)/Hs, where
# Hs is the expected heterozygosity 
homozygous.frequency = matrix(0, num.loc,max.alleles)
for(i in 1:num.loc){
  for(j in 1:max.alleles){
    if(all.freq[i, j] > 0){
      homozygous.frequency[i, j] = homozygous.frequency[i, j] + (all.freq[i, j]) ^ 2
    }
  }
}
fst.max = (sum(homozygous.frequency / num.loc) / (1 - sum(homozygous.frequency) / num.loc))
if(fst.max >= 1){
  fst.max = 1
}
#
# Calculate probability of multilocus genotypes for range of f values
#
# substructure.homozygous calculates the probability of homozygous genotypes
# given a value of f
substructure.homozygous = function(f, all.freq.1){
  all.freq.1 * (f + (1 - f) * all.freq.1)
}
# substructure.heterozygous calculates the probability of homozygous genotypes
# given a value of f
substructure.heterozygous = function(f, all.freq.1, all.freq.2){
  2 * all.freq.1 * all.freq.2 * (1 - f)
}
# Initialise f value
f = 0
# res keeps track of number of number of iterations
res = 1
# Initialise ln.likelihood values (the ln of the probabilities of the 
# individual geneotypes
ln.likelihood = matrix(0, resolution)
# likelihood.locus = the probabilities of the individual genotypes, for
# each individual and each locus
likelihood.locus = matrix(0, N, num.loc)
# The probability of the multicolous dataset is calculated for a range of 
# values of f (0 - fst.max)
repeat{
  likelihood.locus = matrix(0, N, num.loc)
# Take missing genotypes and make likelihood = 1
for(i in 1:N){
  for (j in 1:num.loc){
    if(allele.1[i, j] == 0){
      likelihood.locus[i, j] = 1
    }
  }
}
# Identify homozygous genotypes and calculate probability, given f
  for(i in 1:N){
    for (j in 1:num.loc){
      if(allele.1[i, j] > 0){
        if(allele.1[i, j] == allele.2[i, j]){  
          likelihood.locus[i, j] = substructure.homozygous(f, all.freq[j, allele.1[i, j]])
        }
      }
    }	
  }
# Identify heterozygous genotypes and calculate probability, given f
  for(i in 1:N){
    for (j in 1:num.loc){
      if(allele.1[i, j] > 0){
        if(allele.1[i, j] != allele.2[i, j]){
          likelihood.locus[i, j] = substructure.heterozygous(f, all.freq[j, allele.1[i, j]], all.freq[j, allele.2[i, j]])
        }
      }
    }
  }
# Take log of the product of each individual’s multilocus genotype
  for(i in 1:N){
    likelihood.individual[i] = log(prod(likelihood.locus[i, ]))
  }
# ln.likelihood[res] is the log.probability of the entire dataset, given a value of f
  ln.likelihood[res] = sum(likelihood.individual[, 1])
# To identify the maximum likelihood, the maximum value is 
# subtracted from each of the ln.likelihood values
  max.value = max(ln.likelihood)
  temp = matrix(0, resolution)
  for(i in 1:resolution){
    temp[i]=ln.likelihood[i] - max.value
  }
  for(i in 1:resolution){
    ln.likelihood[i] = temp[i]
  }
  res = res + 1
  f.res = fst.max / resolution
  f = f + f.res
  if(res == (resolution + 1)){
    break
  }
}
#
# The maximum likelihood is then the ln.likelihood value that = 0
# The corresponding value of F is found:
for(i in 1:resolution){
  if(ln.likelihood[i] == 0){
    ML = ((i - 1) / resolution) * fst.max
  }
}
# A plot of the distribution of F; the f-axis goes from 0 - fst.max
f = seq(0, by = f.res, (fst.max - f.res))
f.axis <<- f
#
#
# Likelihood values are converted back from logs
for(i in 1:resolution){
  e.likelihood[i] = exp(ln.likelihood[i])
}
# Area of distribution sums to 1
total = sum(e.likelihood)
e.likelihood = e.likelihood / total
probability <<- e.likelihood
plot(f, e.likelihood, type = "l", xlab = "F",ylab = "Likelihood")
#
#
#
writeLines("Maximum value of Fst = ")
print(fst.max)
writeLines("Maximum Likelihood value of Fst = ")
print(ML)
writeLines("G = ")
G = exp(log(max(e.likelihood))-2)
print(G)
}
# end of function
#  *****************************************************************
#  CONSTRUCT
#  *****************************************************************
construct = function(data, max.alleles,f.resolution,c.resolution,r){
# Input file defined and NAs removed
infile = read.table(data)
infile [is.na(infile)] <- 0
#
# Population size is defined by the number of lines 
N = nrow(infile)
#
# Number of loci is defined by number of columns
num.loc = ncol(infile) / 2
#
# Initialise variables
resolution = 1000 # for estimation of f.max
all.count = matrix(0, num.loc, max.alleles)    
all.freq = matrix(0, num.loc, max.alleles)    
homozygous.frequency = matrix(0, num.loc, max.alleles)     
count.alleles = matrix(0, max.alleles)
likelihood.individual = matrix(0, N)
ln.likelihood = matrix(0, resolution)
e.likelihood = matrix(0, resolution)
hom.total = matrix(0, num.loc)
#
#
# Set up input data
# Read diploid genotypes as sets of two alleles
allele.1 = matrix(0, N, (2 * num.loc))
allele.2 = matrix(0, N, (2 * num.loc))
for(i in 1:N){                                  #  Reads in allele 1 for each locus
  for(j in seq(1, by = 2, (2 * num.loc - 1))){
    allele.1[i, j] = infile[i, j]
  }
}
for(i in 1:N){                                  #  Reads in allele 1 for each locus
  for(j in seq(2, by = 2, (2 * num.loc))){
    allele.2[i, j] = infile[i, j]
  }
}
allele.1.seq = seq(1, 2 * num.loc, by = 2)      #  Order the 1st alleles as odd (1,3,5,etc)
allele.2.seq = seq(2, 2 * num.loc, by = 2)      #  Order the 2nd alleles as even (2,4,6,etc)
for(i in 1:N){
  for (j in 1:num.loc){
    allele.1[i, j] = allele.1[i, allele.1.seq[j]]
  }
}
for(i in 1:N){
  for (j in 1:num.loc){
    allele.2[i, j] = allele.2[i, allele.2.seq[j]]
  }
}
#
#
#
# Calculate allele frequencies
for(i in 1:N){
  for (j in 1:num.loc){
    if (allele.1[i, j] > 0){
      all.count[j, allele.1[i, j]] = all.count[j, allele.1[i, j]] + 1
    }
  }
}
for(i in 1:N){
  for (j in 1:num.loc){
    if (allele.2[i, j] > 0){
      all.count[j, allele.2[i, j]] = all.count[j, allele.2[i, j]] + 1
    }
  }
}
locus.n = matrix(0, num.loc)
for(i in 1:num.loc){
  for(j in 1:max.alleles){
    locus.n[i] = sum(all.count[i, ])
  }
}
for(i in 1:num.loc){
  for(j in 1:max.alleles){
    all.freq[i, j] = all.count[i, j] / locus.n[i]
  }
}
#
# Calculate maximum Fst for two populations according to Hedrick et als
# (1-Hs)/Hs, where Hs is the expected heterozygosity
homozygous.frequency=0
for(i in 1:num.loc){
  for(j in 1:max.alleles){
    if(all.freq[i, j] > 0){
      homozygous.frequency = homozygous.frequency + (all.freq[i, j]) ^ 2
    }
  }
}
homozygous.frequency = homozygous.frequency / num.loc
#
#
#
fst.max = homozygous.frequency / (1 - homozygous.frequency) 
if(fst.max > 1){fst.max = 1}
writeLines("Maximum value of Fst = ")
print(fst.max)
#
#
#
# Calculate probability of multilocus genotypes for range of f values
#
#
# Specify the functions that calculate the probability of genotype frequencies
# For a homozygous genotype in a substructured population where f=fst:
substructure.homozygous = function(f, all.freq.1){
  all.freq.1 * (f + (1 - f) * all.freq.1)
}
# For a heterozygous genotype in a substructured population where f=fst:
substructure.heterozygous = function(f, all.freq.1, all.freq.2){
  2 * all.freq.1 * all.freq.2 * (1 - f)
} 
# For a homozygous genotype in a consanguineous and/or substructured population where f=fst:
consanguinity.homozygous = function(f, r, all.freq.1){
  all.freq.1 * (r + (1 - r) * (f + (1 - f) * all.freq.1))
}
# For a heterozygous genotype in a consanguineous and /or substructured population where f=fst:
consanguinity.heterozygous = function(f, r, all.freq.1, all.freq.2){
  2 * all.freq.1 * all.freq.2 * (1 - r) * (1 - f)
}
#
# Calculate maximum likelihood Fst
f = 0      #initilaise f (fst) as zero
res = 1    #initialise the resolution step (1 - 1000)
# Initialise the likelihood variables
ln.likelihood = matrix(0, resolution)
likelihood.locus = matrix(0, N, num.loc)
likelihood.locus.res = matrix(0, resolution, N, num.loc)
# Iterate through the algorithm for increasing values of f
repeat{
  likelihood.locus = matrix(0, N, num.loc)
# Take missing genotypes and make likelihood = 1
for(i in 1:N){
  for (j in 1:num.loc){
    if(allele.1[i, j] == 0){
      likelihood.locus[i, j] = 1
    }
  }
}
# Identifies homozygotes and calculates probability of genotype, given value of f
  for(i in 1:N){     
    for (j in 1:num.loc){
      if(allele.1[i, j] > 0){
        if(allele.1[i, j] == allele.2[i, j]){
          likelihood.locus[i, j] = substructure.homozygous(f, all.freq[j, allele.1[i, j]])
        }
      }
    }
  }
# Identifies heterozygotes and calculates probability of genotype, given value of f
  for(i in 1:N){
    for (j in 1:num.loc){
      if(allele.1[i, j] > 0){
        if(allele.1[i, j] != allele.2[i, j]){
          likelihood.locus[i, j] = substructure.heterozygous(f, all.freq[j, allele.1[i, j]],all.freq[j, allele.2[i, j]])
        }
      }
    }
  }
# Convert probability of multilocus genotypes to log of probabilities
  for(i in 1:N){
    likelihood.individual[i] = log(prod(likelihood.locus[i, ]))
  }
# Calculate the probability of the total sample, given value of f
  ln.likelihood[res] = sum(likelihood.individual[, 1]) 
# Identify maximum likelihood value of f by subtracting maximum value
# of ln.likelihood from each ln.likelihood
  max.value = max(ln.likelihood)
  temp = matrix(0, resolution)
  for(i in 1:resolution){
    temp[i] = ln.likelihood[i] - max.value
  }
  for(i in 1:resolution){
    ln.likelihood[i] = temp[i]
  }
  res = res + 1
  f.res = fst.max / resolution
  f = f + f.res
  if(res == (resolution + 1)){
    break
  }
}
#
#
# Identify maximum likelihood value of f, which should correspond to 
# the ln.likelihood that = 0
for(i in 1:resolution){
  if(ln.likelihood[i] == 0){
    ML = ((i - 1) / resolution) * fst.max
  }
}
writeLines("Maximum Likelihood value of F = ")
print(ML)
f = seq(0, by = f.res, (fst.max - f.res))
# convert log likelihood values
for(i in 1:resolution){
  e.likelihood[i] = exp(ln.likelihood[i])
}
total=sum(e.likelihood)
#
#
#
# Analyse existing data set: estimating joint likelihood of c & f
#
# Initialise variables
ln.likelihood.c = matrix(0, c.resolution)
ln.likelihood.s = matrix(0, f.resolution)
likelihood.individual.s = matrix(0, N)
likelihood.individual.c = matrix(0, N)
ln.likelihood.individual = matrix(0, N) 
ln.likelihood.cs = matrix(0, f.resolution, c.resolution)
e.likelihood.cs = matrix(0, f.resolution, c.resolution)
#  Likelihood of an individual’s multilocus genotype considering substructure:
likelihood.locus.s.res = matrix(0, f.resolution, N, num.loc)  
#  Likelihood of an individual’s multilocus genotype considering consanguinity:
likelihood.locus.c.res = matrix(0, c.resolution, N, num.loc)  
#
#
#
f.count = 1
f = 0					
#  Iterate through calculation with f ranging from 0 to fst.max
repeat{
  c.count = 1			
# Iterate through c from 0 - 1
  c = 0	
  repeat{
# Initialise the probability of locus specific genotype given f
    likelihood.locus.s = matrix(0, N, num.loc) 
# Initialise the probability of locus specific genotype given f and r
    likelihood.locus.c = matrix(0, N, num.loc)     
# Calculate probability of homozygous genotype given f, r, and allele frequencies	
    for(i in 1:N){                	   
      for (j in 1:num.loc){
        if(allele.1[i, j] > 0){
          if(allele.1[i, j] == allele.2[i, j]){	                 
            likelihood.locus.s[i, j] = substructure.homozygous(f, all.freq[j, allele.1[i, j]])
            likelihood.locus.c[i, j] = consanguinity.homozygous(f, r, all.freq[j, allele.1[i, j]])
          }
        }
      }
    }
# Calculate probability of heterozygous genotype given f, r, and allele frequencies
    for(i in 1:N){
      for (j in 1:num.loc){
        if(allele.1[i, j] > 0){
          if(allele.1[i, j] != allele.2[i, j]){
            likelihood.locus.s[i, j] = substructure.heterozygous(f, all.freq[j, allele.1[i, j]], all.freq[j, allele.2[i, j]])
            likelihood.locus.c[i, j] = consanguinity.heterozygous(f, r, all.freq[j, allele.1[i, j]],all.freq[j, allele.2[i, j]])
          }
        }
      }
    }
# Take product of individual locus genotype probabilities
# and log the values. Set limit to the value the probabilities
# can take, here being 10^(-100), which here is made equivalent 
# to zero.	
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
#
# Log the likelihood of individual’s multilocus genotype given both 
# scenarios (c) = probability of consanguinity
# and (1-c) = probability of substructure
      ln.likelihood.individual[i] = log(c * likelihood.individual.c[i] + (1 - c)*likelihood.individual.s[i]) 
    }
# Log likelihood of whole dataset given values of c and f
#							
    ln.likelihood.cs[f.count, c.count] = sum(ln.likelihood.individual)
    c.count = c.count + 1
    c = c + 1/c.resolution
    if(c.count == (c.resolution + 1)){	
      break
    }
  }
  f.count = f.count + 1      
  f.values = fst.max / f.resolution #f.values are for the plot scale
  f = f + f.values
  if(f.count == (f.resolution + 1)){
    break
  }
}
#
#
# Scale the values such that the maximum likelihood value is set to zero
max.value = max(ln.likelihood.cs)
temp = matrix(0, f.resolution, c.resolution)
for(i in 1:f.resolution){   #subtract the maximum value to set it at zero
  for(j in 1: c.resolution){	
    temp[i, j] = ln.likelihood.cs[i, j] - max.value
  }
}
for(i in 1: f.resolution){
  for(j in 1: c.resolution){	
    ln.likelihood.cs[i, j] = temp[i, j]
  }	
}
# Identify maximum values of Fst and C

for(i in 1: f.resolution){ 
  for(j in 1:c.resolution){
    if(ln.likelihood.cs[i, j] == 0){
      ML.F = ((i - 1) / f.resolution) * fst.max
      ML.C = (j / c.resolution) - 1/c.resolution
    }
  }
}
#
#
#
writeLines("Maximum Likelihood value:")
writeLines("Fst =")
print(ML.F)
writeLines("C = ")
print(ML.C)
#
for(i in 1: f.resolution){	#  Take exponents
  for(j in 1: c.resolution){
    e.likelihood.cs[i, j] = exp(ln.likelihood.cs[i, j])
  }
}
#
total = sum(e.likelihood.cs)
for(i in 1: f.resolution){
  for(j in 1: c.resolution){
    e.likelihood.cs[i, j] = e.likelihood.cs[i, j] / total
  }
}
write.table(e.likelihood.cs,file="ConStruct.Outfile.txt")
x=seq(0, fst.max - f.values, by = f.values)
y=seq(0, 0.99, by = 1/c.resolution)
contour(x, y, e.likelihood.cs,xlab = expression("F"[ST]),ylab = expression("c"[g]),nlevels=4)
writeLines("G = ")
G = exp(log(max(e.likelihood.cs))-2)
print(G)
#
# Global variables FST(f.axis) cg(c.axis) and e.likelihood(probability)
f.axis <<- x
c.axis <<- y
probability <<- e.likelihood.cs
}
# end of function
#  *****************************************************************
#  SIMULATE
#  *****************************************************************
simulate = function(N, num.loc, fst, r.actual, c, r.consider, max.alleles, f.resolution, c.resolution, iteration){
  # Arguments:
  #   N is the total sample size
  #   num.loc is the number of loci
  #   fst is the value of Fst that is to be simulated between two populations
  #   r.actual is the inbreeding coefficient of the inbred individuals
  #   c is the proportion of the population inbred to degree r.actual
  #   r.consider is the value of the inbreeding coefficient being considered for the      
  #     analysis of the simulated dataset 
  #   max.alleles places an uppermost limit on the number of alleles considered
  #   f.resolution is the resolution of the Fst parameter
  #   c.resolution is the resolution on the c parameter
  #   iteration is the number of iterations of the simulation run through in order
  #     to arrive at the specified, simulated Fst
  #  Example of use:
  #  > simulate(200, 12, 0.05, 0.05, 0.5, 0.05, 100, 100, 100, 10000)
  # It is necessary to specify number of alleles per locus
# ***********************************
# specify number of alleles per locus
# ***********************************
  num.alleles = c(8,8,8,8,8,8,8,8,8,8,8,8)
#
# Check to see if specified number of alleles corresponds with the specified number of 
# loci (num.loc):
if(length(num.alleles) != num.loc){
  stop("number of alleles per locus does not match number of loci specified.")
}

# Generate random allele frequencies for each locus for two pops
# x.old.1, for example are the initial allele frequencies in population 1
# x.new.1, for example are the updated allele frequencies in population 1
# Once the required Fst value has been found, these frequencies are referred to 
# as allele.freq.1 and allele.freq.2
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
  f.limit = f.resolution + 1
  c.limit = c.resolution + 1
  sum.allele = array(0, dim = c(num.loc))
  fst.old = array(0, dim = c(num.loc))
  fst.new = array(0, dim = c(num.loc))
  likelihood.individual.c = array(0, dim = c(N))
  likelihood.individual.s = array(0, dim = c(N))
  ln.likelihood.individual = array(0, dim = c(N))
  ln.likelihood.cs = matrix(0, f.resolution, c.resolution)
  e.likelihood.cs = matrix(0, f.resolution, c.resolution)

  num = 0
#
#
# Initial allele frequencies are uniformly distributed

  for(i in 1:num.loc){

    for(j in 1:num.alleles[i]){

      x.old.1[i, j] = 1 / num.alleles[i]

      x.old.2[i, j] = 1 / num.alleles[i]

    }

  }
#
# 
# Specifies maximum number of alleles given by num.alleles	
for(i in 1:num.loc){
  limit[i] = num.alleles[i]
#
#  Generate random allele frequencies until they sum to 1
  while(sum(x.old.1[i, ]) >= 1){
    for(j in 1:limit[i]){
      x.old.1[i, j] = runif(1, 0, 1)  
    }
  }
  x.old.1[i, limit[i]] = 1 - sum(x.old.1[i, 1:limit[i] - 1])
#

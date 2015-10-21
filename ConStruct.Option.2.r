# *******************************
# R Script for ConStruct.Option.2
# By ADJ Overall
# University of Brighton
# Last Update October 2015
# *******************************
#
#
# To run this script, please copy ConStruct.Option.2.r to working directory
# containing your input files. Then type the following lines into R:
# >source("ConStruct.Option.2.r")
#	
# This option reads in a data file of microsatellite diploid genotypes
# Each individual’s multicolous genotype is on a new line.
#
# You will be prompted to input the magnitude of r (degree of consanguinity)
# that you are investigating (e.g., for offspring of 1st cousins r = 0.0625)
#
# Prompt for input file (.txt)
infile = read.table(file.choose())
#
#
#
# Population size is defined by the number of lines 
N = nrow(infile)
#
#
#
# Number of loci is defined by number of columns
num.loc = ncol(infile) / 2
# Maximum number of alleles is defined
max.alleles = 1000
# Resolution of parameter values estimates (here set at 0.01)
resolution = 100
#
#
#
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
for(i in 1:N){                                  #  Reads in allele 1 for each locus  for(j in seq(2, by = 2, (2 * num.loc))){
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
#
#
# Calculate maximum Fst for two populations
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
#
#
f = 0      #initilaise f (fst) as zero
res = 1    #initialise the resolution step (1 - 100)
# Initialise the likelihood variables
ln.likelihood = matrix(0, resolution)
likelihood.locus = matrix(0, N, num.loc)
likelihood.locus.res = matrix(0, resolution, N, num.loc)
# Iterate through the algorithm for increasing values of f
repeat{
  likelihood.locus = matrix(0, N, num.loc)
  for(i in 1:N){     
    for (j in 1:num.loc){
      if(allele.1[i, j] > 0){
        if(allele.1[i, j] == allele.2[i, j]){
          likelihood.locus[i, j] = substructure.homozygous(f, all.freq[j, allele.1[i, j]])
        }
      }
    }
  }
  for(i in 1:N){
    for (j in 1:num.loc){
      if(allele.1[i, j] > 0){
        if(allele.1[i, j] != allele.2[i, j]){
          likelihood.locus[i, j] = substructure.heterozygous(f, all.freq[j, allele.1[i, j]],all.freq[j, allele.2[i, j]])
        }
      }
    }
  }
  for(i in 1:N){
    for (j in 1:num.loc){
      if(allele.1[i, j] == 0){
        likelihood.locus[i, j] = 1
      }
    }
  }
  for(i in 1:N){
    likelihood.individual[i] = log(prod(likelihood.locus[i, ]))
  }
  ln.likelihood[res] = sum(likelihood.individual[, 1]) 
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
  if(res == 101){
    break
  }
}
#
#
# Identify maximum likelihood value
for(i in 1:resolution){
  if(ln.likelihood[i] == 0){
    ML = ((i - 1) / resolution) * fst.max
  }
}
writeLines("Maximum Likelihood value of F = ")
print(ML)
f = seq(0, by = f.res, (fst.max - f.res))
for(i in 1:resolution){
  e.likelihood[i] = exp(ln.likelihood[i])
}
total=sum(e.likelihood)
#
#
#
# Analyse existing data set estimating c & f
#
#
writeLines("Type in the value of consanguinity you wish to consider (eg, 0.125; 0.0625 etc)")
r = readLines(n = 1)
r = as.numeric(r)
resolution = 100
ln.likelihood.c = matrix(0, resolution)
ln.likelihood.s = matrix(0, resolution)
likelihood.individual.s = matrix(0, N)
likelihood.individual.c = matrix(0, N)
ln.likelihood.individual = matrix(0, N) 
ln.likelihood.cs = matrix(0, resolution, resolution)
e.likelihood.cs = matrix(0, resolution, resolution)
likelihood.locus.s.res = matrix(0, resolution, N, num.loc)  #  Likelihood of an individual’s multicolour genotype considering substructure
likelihood.locus.c.res = matrix(0, resolution, N, num.loc)  #  Likelihood of an individual’s multicolour genotype considering consanguinity
# writeLines("Analysis in progress %")  #  Option to note progress of analysis
#
#
#
f.resolution = 1
f = 0					#  Iterate through f from 0 to fst.max
# print(f.resolution - 1)		#  Option to note progress of analysis (see additional print option below)
repeat{
  c.resolution = 1			#  Iterate through c from 0 - 1
  c = 0	
  repeat{
    likelihood.locus.s = matrix(0, N, num.loc) #probability of genotype given f
    likelihood.locus.c = matrix(0, N, num.loc) #probability of genotype given f and r
    for(i in 1:N){                	   #calculate probability of homozygous genotype given f, r, and allele frequencies
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
    for(i in 1:N){		# Probability of homozygous genotype for missing values
      for (j in 1:num.loc){
        if(allele.1[i, j] == 0){
          likelihood.locus.s[i, j] = 1
          likelihood.locus.c[i, j] = 1		
        }
      }
    }
    for(i in 1:N){              #  Take product of individual locus genotype probabilities and log the values	
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
# Log likelihood of individual’s multilocus genotype given both scenarios (c probability of consanguinity
# and (1-c) probability of substructure
#
      ln.likelihood.individual[i] = log(c * likelihood.individual.c[i] + (1 - c)*likelihood.individual.s[i]) 
    }
# Log likelihood of whole dataset given values of c and f
#							
    ln.likelihood.cs[f.resolution, c.resolution] = sum(ln.likelihood.individual)
    c.resolution = c.resolution + 1
    c = c + 0.01
    if(c.resolution == 101){	
      break
    }
  }
# print(f.resolution)         #  Option to note progress of analysis
  f.resolution = f.resolution + 1      
  f.values = fst.max / resolution #f.values are for the plot scale
  f = f + f.values
  if(f.resolution == 101){
    break
  }
}
#
#
# Scale the values such that the maximum likelihood value is set to zero
max.value = max(ln.likelihood.cs)
temp = matrix(0, resolution, resolution)
for(i in 1:resolution){   #subtract the maximum value to set it at zero
  for(j in 1: resolution){	
    temp[i, j] = ln.likelihood.cs[i, j] - max.value
  }
}
for(i in 1: resolution){
  for(j in 1: resolution){	
    ln.likelihood.cs[i, j] = temp[i, j]
  }	
}
#
for(i in 1: resolution){ #identify maximum values of Fst and C
  for(j in 1:resolution){
    if(ln.likelihood.cs[i, j] == 0){
      ML.F = ((i - 1) / resolution) * fst.max
      ML.C = (j / resolution) - 0.01
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
for(i in 1: resolution){	#  Take exponents
  for(j in 1: resolution){
    e.likelihood.cs[i, j] = exp(ln.likelihood.cs[i, j])
  }
}
#
total = sum(e.likelihood.cs)
for(i in 1: resolution){
  for(j in 1: resolution){
    e.likelihood.cs[i, j] = e.likelihood.cs[i, j] / total
  }
}
write.table(e.likelihood.cs,file="ConStruct.Outfile.txt")
x=seq(0, fst.max - f.values, by = f.values)
y=seq(0, 0.99, by = 0.01)
contour(x, y, e.likelihood.cs, xlab = "FST", ylab = "% Consanguinity")
writeLines("G = ")
G = exp(log(max(e.likelihood.cs))-2)print(G)


                                                                                                                                                                                                                                                                                                                                 

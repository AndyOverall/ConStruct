# *******************************
# R Script for ConStruct.Option.1
# By ADJ Overall
# University of Brighton
# Last Update October 2015
# *******************************
#
#
# To run this script, please copy ConStruct.Option.1.r to working directory
# containing your input files. Then type the following lines into R:
# >source("ConStruct.Option.1.r")
#	
# This option reads in a data file of microsatellite diploid genotypes
# Each individual’s multilocus genotype is on a new line
# The maximum likelihood value for F (excess homozygosity) is estimated
#
#
#
# Prompt for input file (.txt)
# infile = matrix(0, (2*num.loc),max.alleles)
infile = read.table(file.choose()) 
#
#
# Population size is defined by the number of lines
N = nrow(infile)
#
#
# Number of loci is defined by number of columns
num.loc = ncol(infile) / 2
#
# Resolution of parameter values estimates (here set at 0.01)
resolution = 100
#
# Maximum number of alleles is defined
max.alleles = 1000
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
#set up input data
#
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
# Calculate maximum Fst
homozygous.frequency = matrix(0, num.loc,max.alleles)
for(i in 1:num.loc){
  for(j in 1:max.alleles){
    if(all.freq[i, j] > 0){
      homozygous.frequency[i, j] = homozygous.frequency[i, j] + (all.freq[i, j]) ^ 2
    }
  }
}
fst.max = (sum(homozygous.frequency / num.loc) / (1 - sum(homozygous.frequency) / num.loc))
#
# Calculate probability of multilocus genotypes for range of f values
substructure.homozygous = function(f, all.freq.1){
  all.freq.1 * (f + (1 - f) * all.freq.1)
}
substructure.heterozygous = function(f, all.freq.1, all.freq.2){
  2 * all.freq.1 * all.freq.2 * (1 - f)
}
#
f = 0
resolution = 100
res = 1
ln.likelihood = matrix(0, resolution)
likelihood.locus = matrix(0, N, num.loc)
likelihood.locus.res = matrix(0, resolution, N, num.loc)
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
          likelihood.locus[i, j] = substructure.heterozygous(f, all.freq[j, allele.1[i, j]], all.freq[j, allele.2[i, j]])
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
    temp[i]=ln.likelihood[i] - max.value
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
for(i in 1:resolution){
  if(ln.likelihood[i] == 0){
    ML = ((i - 1) / resolution) * fst.max
  }
}
f = seq(0, by = f.res, (fst.max - f.res))
#
#
#
for(i in 1:resolution){
  e.likelihood[i] = exp(ln.likelihood[i])
}
total = sum(e.likelihood)
plot(f, e.likelihood / total, type = "l", ylab = "Likelihood")
#
#
#
writeLines("Maximum value of Fst = ")
print(fst.max)
writeLines("Maximum Likelihood value of Fst = ")
print(ML)


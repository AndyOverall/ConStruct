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
while(sum(x.old.2[i, ]) >= 1){
    for(j in 1:limit[i]){
      x.old.2[i, j] = runif(1, 0, 1)
    }
  }	
  x.old.2[i, limit[i]] = 1 - sum(x.old.2[i, 1:limit[i] - 1])	
  

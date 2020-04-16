# HOW TO USE: 
# To easily run many replicates, for example 10 replicates with HGT and 10 with duplications, a simple bash for-loop can be piped to GNU parallel by using the following command:
# for i in {1..10}; do echo python hgt.py ${i} 0.00 0.005; echo python hgt.py ${i} 0.005 0.00 ; done | parallel

import sys,math,os,random,string,errno,copy
import numpy as np
import time as tm

N=2000				# number of individuals
numg=1				# number of gene types, each cell starts with all of them, fully expressed (here 1)
cost_expr=0.01
mut_expr=0.005 		# chance that promoter-strength is mutated
mut_qual=0.00
mut_del=0.005 		# chance to delete a gene
mut_step = 0.2
mutate_up = 0.1	 	# percent of parameter mutations are that are improving (up)
mut_dupl=float(sys.argv[2]) 	# chance that a gene copies
mut_hgt=float(sys.argv[3]) 	# chance that a gene hgt's from someone to you
Tmax=100000
timesavedata=50

intseed = sys.argv[1]
random.seed(intseed)



if mut_dupl > 0.0:
	treatment = 'dup'
else:
	treatment = 'hgt'


path = "%s_%s_genomes" % (treatment,intseed)
if not os.path.exists(path): # Folder that contains genomes (100 first genomes)
	os.mkdir(path)

fitness_contributions = [(1./(20.**i)) for i in range(numg)]


class Cell(object):
	def __init__(self, numg):
		self.genome = list()
		self.fitness = 0.0
		for i in range(numg):
			self.genome.append(Gene(i,1.0,1.0))
			self.time_birth = 0
			self.sumexpr = []       # Total expression of proteins
			self.sumproduct = []    # Product of total expression and quality
			self.cumdups = 0
			self.cumexprmut = 0
			self.cumhgt = 0
			self.cumdel = 0
		self.setfitness()
		self.setexpressions()

	def setfitness(self):
		self.sumexpr = [0. for _ in range(numg) ]
		self.sumproduct = [ 0. for _ in range(numg) ]
		present = [ 0 for _ in range(numg) ]
		for gene in self.genome:
			self.sumexpr[gene.id] += gene.expression
			self.sumproduct[gene.id] += gene.expression*gene.quality
			present[gene.id] = 1
		#if(sum(present)<numg):         # If this statement is enabled, cells that miss a gene, no matter how small the fitness contribution, are not viable
		#	return 0.0
		fitness = 0.0
		for i in range(numg):
			expression = self.sumexpr[i]
			product = self.sumproduct[i]
			fitness -= cost_expr*expression                 # Cost are payed for all the expression
			if(product > 1.0): product = 1.0                # No more than 100% of the fitness contribution can be achieved by expressing a gene
			fitness += fitness_contributions[i]*product
		self.fitness = fitness
		
	def setexpressions(self):
		exprstring = ''
		for gene in self.genome:
			exprstring = exprstring + str(gene.expression) + ','
		exprstring = exprstring[:-1] # remove last comma
		self.expressionpattern = exprstring
	
	def genomesummary(self):
		return str(self.time_birth) + ', ' + str(len(self.genome)) + ', ' + str(self.fitness) + ', ' +  self.expressionpattern + '\n'
		
class Gene(object):
    def __init__(self, id, expression, quality):
		self.id = id
		self.expression = expression
		self.quality = quality

lpop = []
for i in range(N):
	lpop.append(Cell(numg))

time=0
now = int(round(tm.time()*10000))

# Function that saves all the output etc.
def saveoutput():
	
	if(time == 0):
		fdat = open("%s_%s.dat" % (treatment,intseed), "w+")
		fdat.write('Time\tAvgfit\tAvglen\tAvgDup\tAvgHgt\tAvgDel\tAvgExprmut\t')
		for i in range(numg):
			fdat.write('cpn_%d\texpr_%d\tqual_%d\t' % (i,i,i))
		fdat.write('\n')
		fdat.close()
	sumexprs = [ 0. for _ in range(numg) ]
	sumqual = [0. for _ in range(numg) ]
	sumcpn = [ 0. for _ in range(numg) ]
	alllen = [len(cell.genome) for cell in lpop]
	totlen = sum(alllen)
	sumdup = sum([c.cumdups for c in lpop])
	sumhgt = sum([c.cumhgt for c in lpop])
	sumdel = sum([c.cumdel for c in lpop])
	sumexprmut = sum([c.cumexprmut for c in lpop])
	allfit = [cell.fitness for cell in lpop] # (outside of the data-write loop as I used them later for reproduction events)
	totfit = sum(allfit)

	
	num_ancestors = 3 # How many cells are written to file

	tracing = 1		   # How many have been written to file


	fpop = open("%s/t%d_pop.dat" % (path,time), "w+")
	for cell in lpop[1:100]:
		fpop.write(cell.genomesummary())
	fpop.close()
	
	for i in range(numg):
		for cell in lpop:
			for gene in cell.genome:
				if(gene.id == i):
					sumcpn[gene.id]+=1.
					sumexprs[i] += gene.expression
					sumqual[i] += gene.quality
	fdat = open("%s_%s.dat" % (treatment,intseed), "a+")
	fdat.write('%d\t%f\t%f\t%f\t%f\t%f\t%f\t' % (time, totfit/float(N), totlen/float(N), sumdup/float(N), sumhgt/float(N), sumdel/float(N), sumexprmut/float(N)))

	#print 'Average expression of genetypes: '
	for i in range(numg):
		if(sumcpn[i] > 0):
			fdat.write( '%f\t%f\t%f\t' % (sumcpn[i]/float(N), sumexprs[i]/sumcpn[i], sumqual[i]/sumcpn[i]))
		else:
			fdat.write( '0.00\t0.00\t0.00\t')

	fdat.write('\n')
	fdat.close()


# Here is the main loop that includes the moran process, potential events (here, swapping of mutation rates) and saving of output
while time<Tmax:
	# Step 1: output
	if time % timesavedata ==0:
		saveoutput()

	# Step 2: events
	#if(time==50000):
	#	old_dup = mut_dupl
	#	mut_dupl = mut_hgt
	#	mut_hgt = old_dup		
		
	# Step 3: moran process, reproducing N random individuals proportional to fitness
	for i in range(N):
		allfit = [cell.fitness for cell in lpop] # (outside of the data-write loop as I used them later for reproduction events)
		totfit = sum(allfit)
		#Below is the update-step, allowing a random individual to reproduced proportional to its fitness
		rpos = random.uniform(0.,totfit) # number between 0 and sumfit+N (+N ensures fitness isnt only a relative number)
		this=0.
		thepos=-1
		for i,fitness in enumerate(allfit):
			this += fitness
			if this > rpos:
				thepos=i
				break
		#Cell object at index thepos is the cell that reproduces
		pos=random.randrange(N) # where to replicate into (self included)
		newcell = copy.deepcopy(lpop[thepos])
		newcell.time_birth = time
		lpop[pos] = newcell

		mutation = False
		for gene in lpop[pos].genome:
			which = random.uniform(0.,1.) 	# which makes a "roulette wheel" for which mutation happens, proportional to probabilities of each mutations
			if which < mut_expr:
				mutation = True
				if random.uniform(0.,1.) < mutate_up:
					gene.expression = min(1.0,max(0.01,gene.expression + random.uniform(0.0,mut_step))) # Promoter strength not smaller than
				else:
					gene.expression = min(1.0,max(0.01,gene.expression + random.uniform(-mut_step,0.0))) # Promoter strength not smaller than
				lpop[pos].cumexprmut +=1
			elif which < mut_expr+mut_qual:
				mutation = True
				if random.uniform(0.,1.) < mutate_up:
					gene.quality = min(1.0,max(0.0,gene.quality + random.uniform(0.0,mut_step))) # Promoter strength not smaller than
				else:
					gene.quality = min(1.0,max(0.0,gene.quality + random.uniform(-mut_step,0.0))) # Promoter strength not smaller than
			elif which < mut_expr+mut_qual+mut_dupl:
				mutation = True
				newgene = copy.deepcopy(gene)
				lpop[pos].genome.append(newgene)
				lpop[pos].cumdups += 1
			elif which < mut_expr+mut_qual+mut_dupl+mut_hgt:
				mutation = True
				donor = lpop[random.randrange(N)].genome
				if(len(donor)>0):
					donorpos = random.randrange(len(donor))
					newgene = copy.deepcopy(donor[donorpos])
					lpop[pos].genome.append(newgene)
					lpop[pos].cumhgt += 1
			elif which < mut_expr+mut_qual+mut_dupl+mut_hgt+mut_del:
				mutation = True
				lpop[pos].genome.remove(gene)
				lpop[pos].cumdel += 1
		if mutation:
			lpop[pos].setfitness()
			lpop[pos].setexpressions()

	# END OF MORAN PROCESS

	
	# END OF PRINTING OUTPUT
	time+=1

later = int(round(tm.time()*10000))

#print int(later - now)

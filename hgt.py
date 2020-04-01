import sys,math,os,random,string,errno,copy
import numpy as np
import time as tm

N=2000				# number of individuals
numg=1				# number of gene types, each cell starts with all of them, fully expressed
cost_expr=0.01
mut_expr=0.005 		# chance that promoter-strength is mutated
mut_qual=0.00
mut_del=0.005 		# chance to delete a gene
mut_step = 0.2
mutate_up = 0.1	 	# percent of parameter mutations are that are improving (up)
mut_dupl=float(sys.argv[2]) 	# chance that a gene copies
mut_hgt=float(sys.argv[3]) 	# chance that a gene hgt's from someone to you
Tmax=50000
timesavedata=10

#print sys.argv[2]
random.seed(sys.argv[1])
fitness_contributions = [(1./(20.**i)) for i in range(numg)]
#print fitness_contributions
#fitness_contributions = [1 for i in range(numg)]

#print fitness_contributions

ancestors = []

class Cell(object):
	def __init__(self, numg):
		self.genome = list()
		self.fitness = 0.0
		self.parent = None
		for i in range(numg):
			self.genome.append(Gene(i,1.0,1.0))
			self.sumexpr = []       # Total expression of proteins
			self.sumproduct = []    # Product of total expression and quality
			self.cumdups = 0
			self.cumexprmut = 0
			self.cumhgt = 0
			self.cumdel = 0
		self.setfitness()

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

while time<Tmax:
	# START OF MORAN PROCESS, REPEATED N TIMES PER ARBRITRARY UNIT OF TIME
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

		ancestors.append(lpop[pos]) # Store the overwritten cell in the ancestor list.
		newcell.parent = lpop[thepos] # Give the new cell a reference to its parent
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
	# END OF MORAN PROCESS

	# START OF PRINTING OUTPUT
	if time % timesavedata ==0:
		if(time == 0):
			print 'Time\tAvgfit\tAvglen\tAvgDup\tAvgHgt\tAvgDel\tAvgExprmut\t',
			for i in range(numg):
				print 'cpn_%d\texpr_%d\tqual_%d\t' % (i,i,i),
			print ''
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

		for i in range(numg):
			for cell in lpop:
				for gene in cell.genome:
					if(gene.id == i):
						sumcpn[gene.id]+=1.
						sumexprs[i] += gene.expression
						sumqual[i] += gene.quality
		print '%d\t%f\t%f\t%f\t%f\t%f\t%f\t' % (time, totfit/float(N), totlen/float(N), sumdup/float(N), sumhgt/float(N), sumdel/float(N), sumexprmut/float(N)),

		#print 'Average expression of genetypes: '
		for i in range(numg):
			if(sumcpn[i] > 0):
				print '%f\t%f\t%f\t' % (sumcpn[i]/float(N), sumexprs[i]/sumcpn[i], sumqual[i]/sumcpn[i]),
			else:
				print '0.00\t0.00\t0.00\t',

		print ''
	# END OF PRINTING OUTPUT
	time+=1

later = int(round(tm.time()*10000))

#print int(later - now)

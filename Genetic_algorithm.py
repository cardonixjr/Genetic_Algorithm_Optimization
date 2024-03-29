import random as r 						# Library to deal with random numbers
import numpy as np 						# Math functions library 
from plot_3d import *			 		# Library to plot

class Individual():
	def __init__(self, phenotype, mask):
		self.phenotype = phenotype;
		self.mask = mask

		self.genotype = int(phenotype, 2)
		self.fx = 0.0
		self.prob = 0.0
		self.is_elite = False

	def change_prob(self,sum):
		''' 
		Change the crossing chance, where prob is individual fx divided by the sum of all individuals fx.
		'''
		self.prob = self.fx/sum

	def calculate_fx(self,function):
		'''
		Calculates the fx
		'''

		variables = []
		i = 0
		var = ''

		# Break the phenotype in single variables
		for gene in self.phenotype:
			var += gene
			if self.mask[i] == '1':
				variables.append(float(int(var,2)))
				var = ''
			i+=1

		function_caller = "function("
		for var in variables:
			function_caller += str(var) + ","
		
		function_caller = function_caller[:-1] + ")"
		self.fx = eval(function_caller)

	def recalculate_genotype(self):
		self.genotype = int(self.phenotype, 2)


def cross(ind1, ind2, mutation_chance):
	'''
	Function that crosses two individuals phenotypes
	
	Each individual gain half phenotype of other individual
	'''
	new_ind1_phenotype = ""
	new_ind2_phenotype = ""



	for g in range(len(ind1.phenotype)):
		if g%2 :
			new_ind1_phenotype += ind2.phenotype[g]
			new_ind2_phenotype += ind1.phenotype[g]

		else:
			new_ind1_phenotype += ind1.phenotype[g]
			new_ind2_phenotype += ind2.phenotype[g]

	#The elite phenotype don't change 
	if not ind1.is_elite: ind1.phenotype = new_ind1_phenotype
	if not ind2.is_elite: ind2.phenotype = new_ind2_phenotype
	
	# Mutation
	mutation = r.random()
	if mutation < mutation_chance:
		mutated_phenotype_1 = ""
		mutated_phenotype_2 = ""

		mutated_gene = r.randint(0,len(ind1.phenotype)-1)

		# Concatenates the new mutated phenotype
		mutated_phenotype_1 = ind1.phenotype[:mutated_gene]
		mutated_phenotype_2 = ind2.phenotype[:mutated_gene]

		# Invert the mutated gene (0 -> 1 or 1 -> 0)
		mutated_phenotype_1 += str(abs(int(ind1.phenotype[mutated_gene]) - 1))
		mutated_phenotype_2 += str(abs(int(ind2.phenotype[mutated_gene]) - 1))

		# Finish the new mutated phenotype
		mutated_phenotype_1 += ind1.phenotype[mutated_gene+1:]
		mutated_phenotype_2 += ind2.phenotype[mutated_gene+1:]

		# Overwrite the phenotype
		ind1.phenotype = mutated_phenotype_1
		ind2.phenotype = mutated_phenotype_2

def genetic_algorithm_optimization(f, generations, number_of_variables, variables_mask, population_lenght, opt_type = "min"):
	'''
	Function that uses a Genetic Algorithm to optimizate the "f" function.

	f is the function you want to optimize;
	number_of_variables is the number of variables of the "f" function, default is 3 (x,y,z);
	population_lenght is the number os individuals of the population;
	opt_type is the type of optimization. Default is "min", that finds the smallest value of the function;
	
	Other types of optimization are:
	
	"max": find the highest value of the function
	'''
	
	# Creates the initial population
	population = []
	elite = 0
	for pop in range(population_lenght):
		# Creates the phenotype
		phenotype = ''
		for i in range(len(variables_mask)):
			phenotype += str(r.randint(0,1))

		population.append(Individual(phenotype, variables_mask))

	generation = 0
	while generation < generations:
		# Recalculates the individual genotype
		for i in population:
			i.recalculate_genotype()


		# Calculates the function to each individual
		total = 0
		for i in range(population_lenght):
			population[i].calculate_fx(f)
			total += population[i].fx

			# Find the elite individual
			if opt_type == "min":
				if population[i].fx < population[elite].fx:
					population[elite].is_elite = False
					population[i].is_elite = True
					elite = i
	
			if opt_type == "max":
				if population[i].fx > population[elite].fx:
					population[elite].is_elite = False
					population[i].is_elite = True
					elite = i


		# Calculates the crossing probability for each individual
		for i in range(population_lenght):
			population[i].change_prob(total)

		# Cross the individuals
		prob_array = []
		for i in range(population_lenght):
			if opt_type == "min":
				probability = int(round(population[i].prob*100))
			elif opt_type == "max":
				probability = int(round((1-population[i].prob)*100))

			for j in range(probability):
				prob_array.append(i)

		for ind in population:
			choice = r.choice(prob_array)

			# Reroll if the same individual is choiced
			while choice == population.index(ind):
				
				choice = r.choice(prob_array)

			cross(ind, population[choice],0.05)

		generation += 1

	return(population[elite])


if __name__ == "__main__":

	f = lambda x,y:  -((x + 2*y - 7)**2 + (2*x + y - 5)**2)
	#f = lambda x,y: x**2 - y**2

	best = genetic_algorithm_optimization(f, 1000, 2, "00010001", 50, opt_type = "max")
	print("phenotype: " + best.phenotype)
	print("fx: " + str(best.fx))
	my_plot(f, 0, 10, 0, 10, markers = [])
import random

target_organism = '>ti;t{0};target_organism_{0}'
filter_organism = '>ti;f{0};filter_organism_{0}'

def randna(length=100):
	nts = ['A', 'T', 'C', 'G']
	dna = ""
	for i in range(length):
		dna += nts[random.randint(0,3)]
	return dna

fasta = []
with open('srdat/reference.fa') as infile:
	for i, line in enumerate(infile):
		if i > 0:
			fasta.append(line)

with open('refs/target_refs/target_ref_1.fa', 'w') as outfile:	
	for i in range(1, 4):
		outfile.write(target_organism.format(i)+'\n')
		for line in fasta[(i-1)*100:i*100]:
			outfile.write(line)
	for i in range(10, 110):
		outfile.write(target_organism.format(i)+'\n')
		for j in range(100):
			outfile.write(randna(length=100)+'\n')

with open('refs/target_refs/target_ref_2.fa', 'w') as outfile:	
	for i in range(4, 6):
		outfile.write(target_organism.format(i)+'\n')
		for line in fasta[(i-1)*100:i*100]:
			outfile.write(line)
	for i in range(110, 210):
		outfile.write(target_organism.format(i)+'\n')
		for j in range(100):
			outfile.write(randna(length=100)+'\n')

with open('refs/target_refs/target_ref_3.fa', 'w') as outfile:	
	for i in range(6, 10):
		outfile.write(target_organism.format(i)+'\n')
		for line in fasta[(i-1)*100:i*100]:
			outfile.write(line)
	for i in range(210, 310):
		outfile.write(target_organism.format(i)+'\n')
		for j in range(100):
			outfile.write(randna(length=100)+'\n')			

with open('refs/filter_refs/filter_ref_1.fa', 'w') as outfile:	
	for i in range(1, 3):
		outfile.write(filter_organism.format(i)+'\n')
		for line in fasta[(i-1)*100:i*100]:
			outfile.write(line)
	for i in range(10, 20):
		outfile.write(filter_organism.format(i)+'\n')
		for j in range(100):
			outfile.write(randna(length=100)+'\n')

with open('refs/filter_refs/filter_ref_2.fa', 'w') as outfile:	
	for i in range(20, 30):
		outfile.write(filter_organism.format(i)+'\n')
		for j in range(100):
			outfile.write(randna(length=100)+'\n')

with open('refs/filter_refs/filter_ref_3.fa', 'w') as outfile:	
	for i in range(30, 40):
		outfile.write(filter_organism.format(i)+'\n')
		for j in range(100):
			outfile.write(randna(length=100)+'\n')			

#python Comp_Tab.py <input1> <input2> <output> 
#não esquecer de definir e editar a coluna abaixo.


import sys, csv

input_file = sys.argv[1]
input_file2 = sys.argv[2]
input_file3 = sys.argv[3]

file_w = open(input_file3, 'w')

with open(input_file, 'r') as file:
	for line in file:
		col = line.split("\t")
		comp1 = col[0].rstrip()
#		print comp1+"."

		with open (input_file2, 'r') as file2:
			for line2 in file2:
				col2 = line2.split("\t")
				comp2 = col2[1].rstrip()

				if (comp1 == comp2):
					file_w.write(line.rstrip()+"\t"+line2.rstrip()+"\n")
#				print comp2+"."


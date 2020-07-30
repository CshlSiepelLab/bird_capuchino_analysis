import sys
import numpy as np

sample_size = int(sys.argv[1])
number_regions = int(sys.argv[2])
pop_interest = int(sys.argv[3])
input_file = sys.argv[4]
output_file = str(pop_interest) + "_" + input_file

pop_dist = [24, 24]

file = open(output_file, "w") 
file.write("./windowedMSOutput ")
file.write(str(pop_dist[pop_interest]))
file.write(" ")
file.write(str(number_regions))
file.write("\n")
file.write("blah")

def concatenate_list_data(list):
  result= ''
  for element in list:
    result += str(element)
  return result

with open(input_file) as f:
  lines = f.readlines()

i = 4
sample_counter = 1
while i < (sample_size*number_regions + number_regions*4):
  seg_sites = int(lines[i].strip("\n").split(" ")[1])                                 #get number of seg sites
  positions = np.asfarray(lines[i+1].strip("\n").split(" ")[1:(seg_sites+1)])         #get positions
  mat = lines[(i+2):(i+2+sample_size)]
  D = np.zeros(shape=(sample_size, len(positions)))
  for j in range(0,len(mat)):                                                         #construct geno matrix
    D[j] = map(int,list(mat[j].strip("\n")))
  count = 0
  pop_index = 0
  for pop in pop_dist:
    if pop_index==pop_interest:
      temp = D
      temp = temp[count:(count+pop),:]
      seg = np.sum(temp, axis=0)
      index = np.where((seg>0) & (seg<np.shape(temp)[0]))[0]
      mat = temp[:,index]
      pos_pop = positions[index]
      file.write("\n\n//\n")
      file.write("segsites: ")
      file.write(str(np.shape(mat)[1]))
      file.write("\n")
      file.write("positions:")
      for iter in range(len(pos_pop)):
        file.write(" ")
        file.write(str(pos_pop[iter]))
      for iter in range(len(mat)):
        file.write("\n")
        file.write(concatenate_list_data(map(int,mat[iter])))
    pop_index = pop_index + 1
    count = count + pop
  sample_counter = sample_counter + 1
  i = i + sample_size + 4
file.close()

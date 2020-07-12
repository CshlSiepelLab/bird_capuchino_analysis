import sys
import numpy as np

sample_size = 24
sample_number = 1
rep = int(sys.argv[1])
window_1 = int(sys.argv[2])
window_2 = int(sys.argv[3])
type = int(sys.argv[4])
dir = sys.argv[5]
class_ = sys.argv[6]
input_file_1 = class_ + "/" + str(type) + "_" + dir + str(rep) + ".txt_" + str(window_1) + ".msWin"
input_file_2 = class_ + "/" + str(type) + "_" + dir + str(rep) + ".txt_" + str(window_2) + ".msWin"
output_file = class_ + "/LD/" + str(type) + "_" + dir + str(rep) + "_" + str(window_2) + ".txt"
print(input_file_1)
print(input_file_2)

#compute linkage disequilibrium (LD) distribution
def computeLD(window_1, window_2, mat_1, mat_2, pos_pop_1, pos_pop_2):
  seg_sites1 = len(pos_pop_1)
  seg_sites2 = len(pos_pop_2)
  D1 = mat_1
  D2 = mat_2
  i = 0
  if window_1==window_2:
    LDdist = np.zeros(seg_sites1*(seg_sites1-1)/2)
    for k in range(0,seg_sites1):
      for u in range(k+1,seg_sites2):
        l1 = D1[:,k]
        l2 = D2[:,u]
        r = computeLD_freq(l1,l2)
        LDdist[i] = r
        i = i + 1
  else:
    LDdist = np.zeros(seg_sites1*seg_sites2)
    print(seg_sites1)
    #print(seg_sites2)
    for k in range(seg_sites1):
      #print(k)
      for u in range(seg_sites2):
        l1 = D1[:,k]
        l2 = D2[:,u]
        r = computeLD_freq(l1,l2)
        LDdist[i] = r
        i = i + 1
  return(LDdist)

def computeLD_freq(l1,l2):
  freq_p1 = np.count_nonzero(l1)
  freq_p2 = len(l1) - freq_p1
  freq_q1 = np.count_nonzero(l2)
  freq_q2 = len(l2) - freq_q1
  freq_p1 = float(freq_p1) / len(l1)
  freq_p2 = float(freq_p2) / len(l1)
  freq_q1 = float(freq_q1) / len(l2)
  freq_q2 = float(freq_q2) / len(l2)
  freq_p1q1 = float(len(np.intersect1d( np.where(l1==1)[0],  np.where(l2==1)[0]))) / len(l1)
  freq_p1q2 = float(len(np.intersect1d( np.where(l1==1)[0],  np.where(l2==0)[0]))) / len(l1)
  freq_p2q1 = float(len(np.intersect1d( np.where(l1==0)[0],  np.where(l2==1)[0]))) / len(l1)
  freq_p2q2 = float(len(np.intersect1d( np.where(l1==0)[0],  np.where(l2==0)[0]))) / len(l1)
  D_ = freq_p1q1*freq_p2q2 - freq_p1q2*freq_p2q1
  r = D_ / np.sqrt((freq_p1*freq_p2*freq_q1*freq_q2))
  r = np.power(r,2)
  return(r)

def process_data(lines, i, sample_size):
  seg_sites = int(lines[i].strip("\n").split(" ")[1])                                 #get number of seg sites
  positions = np.asfarray(lines[i+1].strip("\n").split(" ")[1:(seg_sites+1)])         #get positions
  mat = lines[(i+2):(i+2+sample_size)]
  D = np.zeros(shape=(sample_size, len(positions)))
  for j in range(0,len(mat)):                                                         #construct geno matrix
    D[j] = map(int,list(mat[j].strip("\n")))
  return([D, seg_sites, positions])

with open(input_file_1) as f:
  lines_1 = f.readlines()

with open(input_file_2) as f:
  lines_2 = f.readlines()

i = 4
sample_counter = 1
LD = []
while i < (sample_size*sample_number + sample_number*4):
  if sample_counter == sample_number:
    [D_1, seg_sites_1, positions_1] = process_data(lines_1, i, sample_size)
    [D_2, seg_sites_2, positions_2] = process_data(lines_2, i, sample_size)
    if np.shape(D_1)[1]<=1 or np.shape(D_2)[1]<=1:
      LD = []
    else:
      LD = computeLD(window_1, window_2, D_1, D_2, positions_1, positions_2)
  sample_counter = sample_counter + 1
  i = i + sample_size + 4
np.savetxt(output_file, LD, fmt = '%1.6f', delimiter=' ')

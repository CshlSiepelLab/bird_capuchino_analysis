import sys
import numpy as np
import pickle
from pandas import DataFrame as df

dir = sys.argv[1]

def get_replicates(features):
  replicate = []
  for i in range(0, np.shape(features)[0]):
    index1 = (len(np.where(features[i,f1]<=maxf1)[0]) + len(np.where(features[i,f2]<=maxf2)[0]) + len(np.where(features[i,f3]<=maxf3)[0]) + len(np.where(features[i,f4]<=maxf4)[0]) + len(np.where(features[i,f5]<=maxf5)[0]) + len(np.where(features[i,f6]<=maxf6)[0]) + len(np.where(features[i,f10]<=maxf10)[0]) + len(np.where(features[i,f11]<=maxf11)[0]))
    index2 = (len(np.where(features[i,f1]>=minf1)[0]) + len(np.where(features[i,f2]>=minf2)[0]) + len(np.where(features[i,f3]>=minf3)[0]) + len(np.where(features[i,f4]>=minf4)[0]) + len(np.where(features[i,f5]>=minf5)[0]) + len(np.where(features[i,f6]>=minf6)[0]) + len(np.where(features[i,f10]>=minf10)[0]) + len(np.where(features[i,f11]>=minf11)[0]))
    if index1==80 & index2==80:
      replicate.append(i) 
  return(replicate)

files = ['unprocessed_neutral.pkl', 'unprocessed_softp2.pkl', 'unprocessed_softp3.pkl', 'unprocessed_hardf2.pkl', 'unprocessed_hardf3.pkl', 'unprocessed_partialp2.pkl', 'unprocessed_partialp3.pkl', 'unprocessed_ancestral.pkl', 'unprocessed_parallel.pkl', dir]

for input_file in files:
  print(input_file)
  with open(input_file, 'rb') as f:
    if input_file == files[0]:
      [neutral, sneutral, gneutral] = pickle.load(f)
    if input_file == files[1]:
      [softp2, ssoftp2, gsoftp2] = pickle.load(f)
    if input_file == files[2]:
      [softp3, ssoftp3, gsoftp3] = pickle.load(f)
    if input_file == files[3]:
      [hardp2, shardp2, ghardp2] = pickle.load(f)
    if input_file == files[4]:
      [hardp3, shardp3, ghardp3] = pickle.load(f)
    if input_file == files[5]:
      [partialp2, spartialp2, gpartialp2] = pickle.load(f)
    if input_file == files[6]:
      [partialp3, spartialp3, gpartialp3] = pickle.load(f)
    if input_file == files[7]:
      [ancestral, sancestral, gancestral] = pickle.load(f)
    if input_file == files[8]:
      [parallel, sparallel, gparallel] = pickle.load(f)
    if input_file == files[9]:
      [mat, smat, gmat] = pickle.load(f)

#110 + 80 = 190
print(np.shape(neutral))
print(np.shape(softp2))
print(np.shape(softp3))
print(np.shape(hardp2))
print(np.shape(hardp3))
print(np.shape(partialp2))
print(np.shape(partialp3))
print(np.shape(ancestral))
print(np.shape(parallel))
print(np.shape(mat))

#110 = 11*5*2 normalized
print(np.shape(sneutral))
print(np.shape(ssoftp2))
print(np.shape(ssoftp3))
print(np.shape(shardp2))
print(np.shape(shardp3))
print(np.shape(spartialp2))
print(np.shape(spartialp3))
print(np.shape(sancestral))
print(np.shape(sparallel))
print(np.shape(smat))

#55 = 11*5 normalized
print(np.shape(gneutral))
print(np.shape(gsoftp2))
print(np.shape(gsoftp3))
print(np.shape(ghardp2))
print(np.shape(ghardp3))
print(np.shape(gpartialp2))
print(np.shape(gpartialp3))
print(np.shape(gancestral))
print(np.shape(gparallel))
print(np.shape(gmat))

#2 populations / 5 windows / 11 features
f1 = np.array([0, 11*1, 11*2, 11*3, 11*4, 11*5, 11*6, 11*7, 11*8, 11*9])
f2 = np.add(f1,1)
f3 = np.add(f1, 2)
f4 = np.add(f1, 3)
f5 = np.add(f1, 4)
f6 = np.add(f1, 5)
f7 = np.add(f1, 6)
f8 = np.add(f1, 7)
f9 = np.add(f1, 8)
f10 = np.add(f1, 9)
f11 = np.add(f1, 10)

maxf1 = max(np.amax(mat[:,f1], axis=0))
maxf2 = max(np.amax(mat[:,f2], axis=0))
maxf3 = max(np.amax(mat[:,f3], axis=0))
maxf4 = max(np.amax(mat[:,f4], axis=0))
maxf5 = max(np.amax(mat[:,f5], axis=0))
maxf6 = max(np.amax(mat[:,f6], axis=0))
maxf7 = max(np.amax(mat[:,f7], axis=0))
maxf8 = max(np.amax(mat[:,f8], axis=0))
maxf9 = max(np.amax(mat[:,f9], axis=0))
maxf10 = max(np.amax(mat[:,f10], axis=0))
maxf11 = max(np.amax(mat[:,f11], axis=0))

minf1 = min(np.amin(mat[:,f1], axis=0))
minf2 = min(np.amin(mat[:,f2], axis=0))
minf3 = min(np.amin(mat[:,f3], axis=0))
minf4 = min(np.amin(mat[:,f4], axis=0))
minf5 = min(np.amin(mat[:,f5], axis=0))
minf6 = min(np.amin(mat[:,f6], axis=0))
minf7 = min(np.amin(mat[:,f7], axis=0))
minf8 = min(np.amin(mat[:,f8], axis=0))
minf9 = min(np.amin(mat[:,f9], axis=0))
minf10 = min(np.amin(mat[:,f10], axis=0))
minf11 = min(np.amin(mat[:,f11], axis=0))

replicate1 = get_replicates(neutral)
neutral = neutral[replicate1,:]
sneutral = sneutral[replicate1,:]
gneutral = gneutral[replicate1,:]

replicate2 = get_replicates(softp2)
softp2 = softp2[replicate2,:]
ssoftp2 = ssoftp2[replicate2,:]
gsoftp2 = gsoftp2[replicate2,:]

replicate3 = get_replicates(softp3)
softp3 = softp3[replicate3,:]
ssoftp3 = ssoftp3[replicate3,:]
gsoftp3 = gsoftp3[replicate3,:]

replicate4 = get_replicates(hardp2)
hardp2 = hardp2[replicate4,:]
shardp2 = shardp2[replicate4,:]
ghardp2 = ghardp2[replicate4,:]

replicate5 = get_replicates(hardp3)
hardp3 = hardp3[replicate5,:]
shardp3 = shardp3[replicate5,:]
ghardp3 = ghardp3[replicate5,:]

replicate6 = get_replicates(partialp2)
partialp2 = partialp2[replicate6,:]
spartialp2 = spartialp2[replicate6,:]
gpartialp2 = gpartialp2[replicate6,:]

replicate7 = get_replicates(partialp3)
partialp3 = partialp3[replicate7,:]
spartialp3 = spartialp3[replicate7,:]
gpartialp3 = gpartialp3[replicate7,:]

replicate8 = get_replicates(ancestral)
ancestral = ancestral[replicate8,:]
sancestral = sancestral[replicate8,:]
gancestral = gancestral[replicate8,:]

replicate9 = get_replicates(parallel)
parallel = parallel[replicate9,:]
sparallel = sparallel[replicate9,:]
gparallel = gparallel[replicate9,:]

sneutral = np.hstack((sneutral, gneutral))
ssoftp2 = np.hstack((ssoftp2, gsoftp2))
ssoftp3 = np.hstack((ssoftp3, gsoftp3))
shardp2 = np.hstack((shardp2, ghardp2))
shardp3 = np.hstack((shardp3, ghardp3))
spartialp2 = np.hstack((spartialp2, gpartialp2))
spartialp3 = np.hstack((spartialp3, gpartialp3))
sancestral = np.hstack((sancestral, gancestral))
sparallel = np.hstack((sparallel, gparallel))

sneutral = np.hstack((sneutral, neutral[:,110:190]))
ssoftp2 = np.hstack((ssoftp2, softp2[:,110:190]))
ssoftp3 = np.hstack((ssoftp3, softp3[:,110:190]))
shardp2 = np.hstack((shardp2, hardp2[:,110:190]))
shardp3 = np.hstack((shardp3, hardp3[:,110:190]))
spartialp2 = np.hstack((spartialp2, partialp2[:,110:190]))
spartialp3 = np.hstack((spartialp3, partialp3[:,110:190]))
sancestral = np.hstack((sancestral, ancestral[:,110:190]))
sparallel = np.hstack((sparallel, parallel[:,110:190]))

print(np.shape(sneutral))
print(np.shape(ssoftp2))
print(np.shape(ssoftp3))
print(np.shape(shardp2))
print(np.shape(shardp3))
print(np.shape(spartialp2))
print(np.shape(spartialp3))
print(np.shape(sancestral))
print(np.shape(sparallel))

smat = np.hstack((smat, gmat))
smat = np.hstack((smat, mat[:,110:190]))

output_file = 'processed_' + dir
with open(output_file, 'wb') as f:
  pickle.dump([sneutral, ssoftp2, ssoftp3, shardp2, shardp3, spartialp2, spartialp3, sancestral, sparallel, smat], f)

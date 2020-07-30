import sys
import numpy as np
from scipy.stats import kurtosis, skew
import pickle
from pandas import DataFrame as df

dir = sys.argv[1]
_id = int(sys.argv[2])

num_populations = 2
num_windows = 5
features = 11
features_LD = 8

num_features = features*num_windows*num_populations + features_LD*num_windows*num_populations #+ features_LD*num_windows
mat = np.empty((0,num_features), float)
ids = [_id]

for id in ids:
  temp = []
  temp_LD = []
  for j in range(0,num_populations):
    for window in range(0, num_windows):
      filename = dir + "/features/" + "feature_" + str(window) + "_" + str(j) + "_" + str(id) + ".txt.stats"
      data = np.loadtxt(filename, skiprows=1)
      if j==0 and window==0:
        temp = data
      else:
        temp = np.hstack((temp, data))
  num_samples = np.shape(temp)[0]
  numLDfeatures = features_LD*num_windows*num_populations #+ features_LD*num_windows
  mat_LD = np.empty((0,numLDfeatures), float)
  for replicate in range(1,num_samples+1):
    print(replicate)
    tempLD = []
    for j in range(0,num_populations):
      for i in range(0,num_windows):
        filename = dir + "/LD/feature_" + str(j) + "_" + str(i) + "_" + str(id) + "_" + str(replicate) + "_LD.txt"
        data = np.loadtxt(filename, ndmin=1)
        tempLD = np.append(tempLD, [np.mean(data), float(df(data).quantile(.05)), float(df(data).quantile(.1)), float(df(data).quantile(.25)), float(df(data).quantile(.5)), float(df(data).quantile(.75)), float(df(data).quantile(.9)), float(df(data).quantile(.95))])
    #global_stats = []
    for window in range(0, num_windows):
      filename = dir + "/LDG/feature_" + str(window) + "_" + str(id) + "_" + str(replicate) + "_LD.txt"
      #data = np.loadtxt(filename)
      #if len(np.shape(data))==0:
      #  global_stats = np.append(global_stats, [data, data, data, data, data])
      #else:
      #  global_stats = np.append(global_stats, [np.mean(data), float(df(data).quantile(.05)), float(df(data).quantile(.1)), float(df(data).quantile(.9)), float(df(data).quantile(.95))])
    #mat_LD = np.vstack((mat_LD, np.append(tempLD,global_stats)))
    mat_LD = np.vstack((mat_LD, tempLD))
  mat = np.concatenate((mat, np.hstack((temp,mat_LD))), axis=0)
  print(np.shape(mat))

smat = np.empty((0,features*num_windows*num_populations),float)
for id in ids:
  temp = []
  for j in range(0,num_populations):
    filename = dir + "/features/" + "feature_" + str(j) + "_" + str(id) + ".txt"
    data = np.loadtxt(filename, skiprows=11)
    if j==0:
      temp = data
    else:
      temp = np.hstack((temp, data))
  smat = np.vstack((smat, temp))


gmat = np.empty((0,features*num_windows),float)
for id in ids:
  filename = dir + "/features/" + "gfeature_" + str(id) + ".txt"
  data = np.loadtxt(filename, skiprows=11)
  gmat = np.vstack((gmat, data))

#Fst = np.empty((0,num_windows),float)
#for id in ids:
#  filename = "empirical/Fst/" + str(id) + ".txt"
#  data = np.loadtxt(filename)
#  Fst = np.vstack((Fst, data))

print(np.shape(mat))
print(np.shape(smat))
print(np.shape(gmat))
#print(np.shape(Fst))

output_file = 'unprocessed_' + str(_id) + "_" + dir
with open(output_file, 'wb') as f:
  pickle.dump([mat, smat, gmat], f)

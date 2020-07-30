import sys
import numpy as np
import pickle
from pandas import DataFrame as df

num_populations = 2
num_windows = 5
features = 11
features_LD = 8

dir = sys.argv[1]
type = sys.argv[2]
num_samples = int(sys.argv[3])

def sloaddata(subfolder1, class_, num_populations, num_samples, num_features):
  mat = np.empty((0,num_features), float)
  for replicate in range(1,num_samples+1):
    temp = []
    for j in range(0,num_populations):
      filename = subfolder1 + "/" + str(j) + '_' + class_ + str(replicate) + ".txt"
      data = np.loadtxt(filename, skiprows=11)
      temp = np.hstack((temp, data))
    mat = np.vstack((mat, temp))
  return(mat)

def gloaddata(subfolder1, class_, num_samples, num_features):
  mat = np.empty((0,num_features), float)
  for replicate in range(1,num_samples+1):
    filename = subfolder1 + "/g" + class_ + str(replicate) + ".txt"
    data = np.loadtxt(filename, skiprows=11)
    mat = np.vstack((mat, data))
  return(mat)

num_features = features*num_windows*num_populations
Y = sloaddata(dir, type, num_populations, num_samples, num_features)
num_features = features*num_windows
G = gloaddata(dir, type, num_samples, num_features)

def loaddata(subfolder1, class_, num_populations, num_samples, num_features, num_windows):
  replicate_vec = []
  mat = np.empty((0,num_features), float)
  for replicate in range(1,num_samples+1):
    print(replicate)
    temp = []
    tempLD = []
    flag = 1
    for j in range(0,num_populations):
      for window in range(0, num_windows):
        filename = subfolder1 + "/" + str(j) + '_' + class_ + str(replicate) + ".txt_" + str(window) + ".msWin.stats"      
        data = np.loadtxt(filename, skiprows=1)
        temp = np.hstack((temp, data))
        filename = subfolder1 + "/LD/" + str(j) + "_" + class_ + str(replicate) + "_" + str(window) + ".txt"
        data = np.loadtxt(filename, ndmin=1)
        if len(data)>1:
          tempLD = np.append(tempLD, [np.mean(data), float(df(data).quantile(.05)), float(df(data).quantile(.1)), float(df(data).quantile(.25)), float(df(data).quantile(.5)), float(df(data).quantile(.75)), float(df(data).quantile(.9)), float(df(data).quantile(.95))])
        else:
          flag = 0
    if flag == 1:
      global_stats = []
      for window in range(0, num_windows):
        filename = subfolder1 + "/LDG/" + class_ + str(replicate) + "_" + str(window) + ".txt"
        #data = np.loadtxt(filename, ndmin=1)
        #global_stats = np.append(global_stats, [np.mean(data), float(df(data).quantile(.05)), float(df(data).quantile(.1)), float(df(data).quantile(.25)), float(df(data).quantile(.5)), float(df(data).quantile(.75)), float(df(data).quantile(.9)), float(df(data).quantile(.95))])
      #tempLD = np.append(tempLD, global_stats)
      replicate_vec = np.append(replicate_vec, replicate)
      mat = np.vstack((mat, np.hstack((temp,tempLD))))
  return(mat,replicate_vec)

num_features = features*num_windows*num_populations + features_LD*num_windows*num_populations #+ features_LD*num_windows
print(num_features)
X, rep_X = loaddata(dir, type, num_populations, num_samples, num_features, num_windows)

Y = Y[rep_X.astype(int)-1,]
G = G[rep_X.astype(int)-1,]

print(np.shape(Y))
print(np.shape(G))
print(np.shape(X))

output_file = 'unprocessed_' + dir +'.pkl'
with open(output_file, 'wb') as f:
  pickle.dump([X, Y, G], f)

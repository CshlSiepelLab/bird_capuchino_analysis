import sys
import numpy as np
from sklearn.model_selection import cross_val_score
from sklearn import svm
from sklearn.metrics import confusion_matrix
import pickle

type_ = sys.argv[1]
lambda_= float(sys.argv[2])

input_file = 'processed_unprocessed_empirical_' + type_
with open(input_file, 'rb') as f:  # Python 3: open(..., 'rb')
  neutral, softp2, softp3, hardp2, hardp3, partialp2, partialp3, ancestral, parallel, mat = pickle.load(f)

print(np.shape(neutral))
print(np.shape(softp2))
print(np.shape(softp3))
print(np.shape(hardp2))
print(np.shape(hardp3))

training_rep1 = 8000
test_rep1 = training_rep1 + 1000

training_rep2 = 8000
test_rep2 = training_rep2 + 1000

training_rep3 = 8000
test_rep3 = training_rep3 + 1000

training_rep4 = 8000
test_rep4 = training_rep4 + 1000

training_rep5 = 8000
test_rep5 = training_rep5 + 1000

INDEX_ = []

INDEX_1 = range(0,5)			#pi
INDEX_2 = range(5,10)			#ss
INDEX_3 = range(10,15)			#thetaH
INDEX_4 = range(15,20)			#tajD
INDEX_5 = range(20,25)			#fayWuH
INDEX_6 = range(25,30)			#HapCount
INDEX_7 = range(45,50)			#Omega
INDEX_8 = range(50,55)			#ZnS

INDEX_ = np.append(INDEX_, INDEX_1)
INDEX_ = np.append(INDEX_, INDEX_2)
INDEX_ = np.append(INDEX_, INDEX_3)
INDEX_ = np.append(INDEX_, INDEX_4)
INDEX_ = np.append(INDEX_, INDEX_5)
INDEX_ = np.append(INDEX_, INDEX_6)
INDEX_ = np.append(INDEX_, INDEX_7)
INDEX_ = np.append(INDEX_, INDEX_8)

num=55
INDEX_1 = range(0+num,5+num)		#pi
INDEX_2 = range(5+num,10+num)		#ss
INDEX_3 = range(10+num,15+num)		#thetaH
INDEX_4 = range(15+num,20+num)		#tajD
INDEX_5 = range(20+num,25+num)		#fayWuH
INDEX_6 = range(25+num,30+num)		#HapCount
INDEX_7 = range(45+num,50+num)		#Omega
INDEX_8 = range(50+num,55+num)		#ZnS

INDEX_ = np.append(INDEX_, INDEX_1)
INDEX_ = np.append(INDEX_, INDEX_2)
INDEX_ = np.append(INDEX_, INDEX_3)
INDEX_ = np.append(INDEX_, INDEX_4)
INDEX_ = np.append(INDEX_, INDEX_5)
INDEX_ = np.append(INDEX_, INDEX_6)
INDEX_ = np.append(INDEX_, INDEX_7)
INDEX_ = np.append(INDEX_, INDEX_8)

num=110
INDEX_1 = range(0+num,5+num)            #pi
INDEX_2 = range(5+num,10+num)           #ss
INDEX_3 = range(10+num,15+num)          #thetaH
INDEX_4 = range(15+num,20+num)          #tajD
INDEX_5 = range(20+num,25+num)          #fayWuH
INDEX_6 = range(25+num,30+num)          #HapCount
INDEX_7 = range(45+num,50+num)          #Omega
INDEX_8 = range(50+num,55+num)          #ZnS

INDEX_ = np.append(INDEX_, INDEX_1)
INDEX_ = np.append(INDEX_, INDEX_2)
INDEX_ = np.append(INDEX_, INDEX_3)
INDEX_ = np.append(INDEX_, INDEX_4)
INDEX_ = np.append(INDEX_, INDEX_5)
INDEX_ = np.append(INDEX_, INDEX_6)
INDEX_ = np.append(INDEX_, INDEX_7)
INDEX_ = np.append(INDEX_, INDEX_8)

num_LD=8
begin_=165
base = []
for i in range(0,10):
  base = np.append(base, begin_+num_LD*i)

INDEX_3 = []
features_to_use = [5,6,7]
for i in features_to_use:
  INDEX_3 = np.concatenate((INDEX_3, base+i))
INDEX_3 = np.sort(INDEX_3)
INDEX_ = np.append(INDEX_, INDEX_3)

s = 0 
i1 = [INDEX_3[s], INDEX_3[s+3], INDEX_3[s+6], INDEX_3[s+9], INDEX_3[s+12]]
s = 1 
i2 = [INDEX_3[s], INDEX_3[s+3], INDEX_3[s+6], INDEX_3[s+9], INDEX_3[s+12]]
s = 2
i3 = [INDEX_3[s], INDEX_3[s+3], INDEX_3[s+6], INDEX_3[s+9], INDEX_3[s+12]]

s = 15
i4 = [INDEX_3[s], INDEX_3[s+3], INDEX_3[s+6], INDEX_3[s+9], INDEX_3[s+12]]
s = 16
i5 = [INDEX_3[s], INDEX_3[s+3], INDEX_3[s+6], INDEX_3[s+9], INDEX_3[s+12]]
s = 17
i6 = [INDEX_3[s], INDEX_3[s+3], INDEX_3[s+6], INDEX_3[s+9], INDEX_3[s+12]]

i1 = np.int32(i1)
i2 = np.int32(i2)
i3 = np.int32(i3)
i4 = np.int32(i4)
i5 = np.int32(i5)
i6 = np.int32(i6)

for IND in [i1, i2, i3, i4, i5, i6]:
  print(IND) 
  neutral[:,IND] = neutral[:,IND]/np.reshape(np.sum(neutral[:,IND], axis=1), (np.shape(neutral)[0], 1))
  softp2[:,IND] = softp2[:,IND]/np.reshape(np.sum(softp2[:,IND], axis=1), (np.shape(softp2)[0], 1))
  softp3[:,IND] = softp3[:,IND]/np.reshape(np.sum(softp3[:,IND], axis=1), (np.shape(softp3)[0], 1))
  hardp2[:,IND] = hardp2[:,IND]/np.reshape(np.sum(hardp2[:,IND], axis=1), (np.shape(hardp2)[0], 1))
  hardp3[:,IND] = hardp3[:,IND]/np.reshape(np.sum(hardp3[:,IND], axis=1), (np.shape(hardp3)[0], 1))

print(len(INDEX_))
INDEX_ = np.int32(INDEX_)
neutral = neutral[:,INDEX_]
softp2 = softp2[:,INDEX_]
softp3 = softp3[:,INDEX_]
hardp2 = hardp2[:,INDEX_]
hardp3 = hardp3[:,INDEX_]
mat = mat[:,INDEX_]

print(np.shape(neutral))
print(np.shape(softp2))
print(np.shape(softp3))
print(np.shape(hardp2))
print(np.shape(hardp3))

fv = np.vstack((neutral[0:training_rep1,], softp2[0:training_rep2,], softp3[0:training_rep3,], hardp2[0:training_rep4,], hardp3[0:training_rep5,]))

n = np.zeros(int(np.shape(neutral[0:training_rep1,])[0]))
s1 = np.ones(int(np.shape(softp2[0:training_rep2,])[0]))
s2 = np.ones(int(np.shape(softp3[0:training_rep3,])[0])) + 1
h1 = np.ones(int(np.shape(hardp2[0:training_rep4,])[0])) + 2
h2 = np.ones(int(np.shape(hardp3[0:training_rep5,])[0])) + 3
target = np.append(n,s1)
target = np.append(target,s2)
target = np.append(target,h1)
target = np.append(target,h2)

clf = svm.SVC(kernel='linear', C=lambda_, probability=True)
#scores = cross_val_score(clf, fv, target, cv=10)
#print(scores)
print(np.shape(fv))
clf.fit(fv, target)

fv1 = np.vstack((neutral[training_rep1:test_rep1,], softp2[training_rep2:test_rep2,], softp3[training_rep3:test_rep3,], hardp2[training_rep4:test_rep4,], hardp3[training_rep5:test_rep5,]))
n = np.zeros(int(np.shape(neutral[training_rep1:test_rep1,])[0]))
s1 = np.ones(int(np.shape(softp2[training_rep2:test_rep2,])[0]))
s2 = np.ones(int(np.shape(softp3[training_rep3:test_rep3,])[0])) + 1
h1 = np.ones(int(np.shape(hardp2[training_rep4:test_rep4,])[0])) + 2
h2 = np.ones(int(np.shape(hardp3[training_rep5:test_rep5,])[0])) + 3
target = np.append(n,s1)
target = np.append(target,s2)
target = np.append(target,h1)
target = np.append(target,h2)

pred = clf.predict(fv1)
pred1 = clf.predict_proba(fv1)

print(len(pred[pred==target])/float(len(pred)))
print(confusion_matrix(target, pred))

output = "proba_hardsoft" + type_ + ".txt"
np.savetxt(output, pred1, fmt = '%f', delimiter=' ')

import pickle
output_file = 'model_hardsoft' + type_ + '.pkl'
with open(output_file, 'wb') as f:  # Python 3: open(..., 'wb')
  pickle.dump([clf, pred1], f)

names = []
features = np.array(["local_pi", "local_seg", "local_thetaH", "local_tajima", "local_fayWuH", "local_HapCount", "local_Omega", "local_ZnS"], dtype=np.object)
for i in range(1,3):
  for j in features:
    for num in range(1,6):
      names = np.append(names, j+"_"+str(num)+"_"+str(i))

features = np.array(["global_pi", "global_seg", "global_thetaH", "global_tajima", "global_fayWuH", "global_HapCount", "global_Omega", "global_ZnS"], dtype=np.object)
for j in features:
  for num in range(1,6):
    names = np.append(names, j+"_"+str(num))

features = np.array(["local_quantile75", "local_quantile90", "local_quantile95"], dtype=np.object)
for i in range(1,3):
  for num in range(1,6):
    for j in features:
      names = np.append(names, j+"_"+str(num)+"_"+str(i))

print(len(clf.coef_[0]))
print(len(names))
print(len(clf.coef_))
for i in range(0,len(clf.coef_)):
  imp,fnames = zip(*sorted(zip(abs(clf.coef_[i]),names)))
  print(imp[(len(imp)-10):len(imp)])
  print(fnames[(len(imp)-10):len(imp)])

pred1 = clf.predict_proba(mat)
output = "empirical_hardsoft_" + type_ + ".txt"
np.savetxt(output, pred1, fmt = '%f', delimiter=' ')

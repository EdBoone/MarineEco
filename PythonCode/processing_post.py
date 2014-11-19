from pylab import *

f1 = open('posterior_points')
X = []
Y = []
W = []
temp = str()

for line in f1:
    if line[0] == '[':
        temp = line.strip()
    else:
        temp = temp +'  '+ line.strip()
    
    if temp[-1] == ']':
        temp = temp[1:-1]
        temp = array(temp.split(),dtype=float)
        X.append(temp)
        temp = str()
f1.close()

f2 = open('posterior_weights')

for line in f2:
    W.append(float(line))
f2.close()

X = array(X)
W = array(W)

weighted_mean = mean(X.T*W,1)*(len(W)/sum(W))

f3 = open('weighted_mean_142','w+')
for i in weighted_mean:
    f3.write(str(i)+'\n')
f3.close()
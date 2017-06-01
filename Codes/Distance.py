from sklearn import preprocessing
from sklearn.neighbors import DistanceMetric 
import pandas as pd
import numpy as np
from ete3 import Tree
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.pairwise import cosine_distances
import scipy
import pylab
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as sd
import matplotlib.pyplot as plt
import matplotlib.collections as collections

#Functions for flattening the lists and graphics###################################
def flatten(vari,arr1,arr2,arr3):
	for x in arr1:
		vari.append(x)
	for x in arr2:
		vari.append(x)
	for x in arr3:
		vari.append(x)		
	return vari

def condFlatten(vari,arr):
	for x in arr:
		vari.append(x)		
	return vari
	
def triatpos(pos=(0,0), rot=0):
    r = np.array([[-1,-1],[1,-1],[1,1],[-1,-1]])*.5
    rm = [[np.cos(np.deg2rad(rot)), -np.sin(np.deg2rad(rot))],
           [np.sin(np.deg2rad(rot)),np.cos(np.deg2rad(rot)) ] ]
    r = np.dot(rm, r.T).T
    r[:,0] += pos[0]
    r[:,1] += pos[1]
    return r

def triamatrix(a, ax, rot=0, cmap=plt.cm.viridis, **kwargs):
    segs = []
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            segs.append(triatpos((j,i), rot=rot) )
    col = collections.PolyCollection(segs, cmap=cmap, **kwargs)
    col.set_array(a.flatten())
    ax.add_collection(col)
    return col
##################################################################################


#Distance storage arrays##########################################################
GC = [] 
aa = []
nucl = []
dinucl = []
diaa = []
##################################################################################


#Distance lists : convention = gENOME#############################################
g1 = []
g2 = []
g3 = []
g4 = []
g5 = []
##################################################################################


#File I/O pointers
files = ["anal_out.txt","diaa.txt","dinucl.txt","nucl.txt","gc.txt"]
##################################################################################


#Dataframes for individual files
dfGC = pd.DataFrame.from_csv('GC.txt', sep='\t', header=0)
dfAA = pd.DataFrame.from_csv('anal_out.txt', sep='\t', header=0)
dfDiAA = pd.DataFrame.from_csv('diaa.txt', sep='\t', header=0)
dfNU = pd.DataFrame.from_csv('nucl.txt', sep='\t', header=0)
dfDiNU = pd.DataFrame.from_csv('dinucl.txt', sep='\t', header=0)


#"Correcting" dataframe
dfDiAA.index = dfDiAA.index + dfDiAA.Frame.map(str) 
dfAA.index = dfAA.index + dfAA.Frame.map(str) 
##################################################################################


#Normalising the dataframes#######################################################
#Method 1
#mean = dfGC.GC.dropna().mean()
#max = dfGC.GC.dropna().max()
#min = dfGC.GC.dropna().min()
#dfGC['GC'] = dfGC['GC'].apply(lambda x: (x - mean) / (max -min))

#Method 2
#x = dfGC.values #returns a numpy array
#min_max_scaler = preprocessing.MinMaxScaler()
#x_scaled = min_max_scaler.fit_transform(x)
#dfGC = pd.DataFrame(x_scaled)

#Method 3 - https://stackoverflow.com/questions/40197156/python-pandas-best-way-to-normalize-data
dfAA = (dfAA - dfAA.mean())/dfAA.std()
dfDiAA = (dfDiAA - dfDiAA.mean())/dfDiAA.std()
dfGC = (dfGC - dfGC.mean())/dfGC.std()
dfNU = (dfNU - dfNU.mean())/dfNU.std()
dfDiNU = (dfDiNU - dfDiNU.mean())/dfDiNU.std()
##################################################################################


#Creating Meta arrays to append distance information############################## 
#For DiAA
temp = dfDiAA.values.tolist()
g1 = flatten(g1,temp[0],temp[1],temp[2])
g2 = flatten(g2,temp[3],temp[4],temp[5])
g3 = flatten(g3,temp[6],temp[7],temp[8])
g4 = flatten(g4,temp[9],temp[10],temp[11])
g5 = flatten(g5,temp[12],temp[13],temp[14]) 

#For AA
temp = dfAA.values.tolist()
g1 = flatten(g1,temp[0],temp[1],temp[2])
g2 = flatten(g2,temp[3],temp[4],temp[5])
g3 = flatten(g3,temp[6],temp[7],temp[8])
g4 = flatten(g4,temp[9],temp[10],temp[11])
g5 = flatten(g5,temp[12],temp[13],temp[14]) 

#For NU
temp = dfNU.values.tolist()
g1 = condFlatten(g1,temp[0])
g2 = condFlatten(g2,temp[1])
g3 = condFlatten(g3,temp[2])
g4 = condFlatten(g4,temp[3])
g5 = condFlatten(g5,temp[4])

#For diNU
temp = dfDiNU.values.tolist()
g1 = condFlatten(g1,temp[0])
g2 = condFlatten(g2,temp[1])
g3 = condFlatten(g3,temp[2])
g4 = condFlatten(g4,temp[3])
g5 = condFlatten(g5,temp[4])

#For GC
temp = dfDiNU.values.tolist()
g1 = condFlatten(g1,temp[0])
g2 = condFlatten(g2,temp[1])
g3 = condFlatten(g3,temp[2])
g4 = condFlatten(g4,temp[3])
g5 = condFlatten(g5,temp[4])
##################################################################################


#Numpy operations#################################################################
g1 = np.array(g1)
n, bins, patches = plt.hist(g1)
plt.show()
g2 = np.array(g2)
g3 = np.array(g3)
g4 = np.array(g4)
g5 = np.array(g5)
X = (g1,g2,g3,g4,g5)
np.set_printoptions(threshold=np.inf, linewidth=np.inf)
##################################################################################


#Distance matrices################################################################  
dist = DistanceMetric.get_metric('minkowski',p=1)
distance = dist.pairwise(X)
np.savetxt('minkowski_p1.txt', distance, delimiter=',')

with open("minkowski_p1", 'w') as f:
	f.write(np.array2string(distance, separator=', '))
	
dist = DistanceMetric.get_metric('minkowski',p=2)
distance = dist.pairwise(X)

with open("minkowski_p2.txt", 'w') as f:
	f.write(np.array2string(distance, separator=', '))
	
dist = DistanceMetric.get_metric('euclidean')
distance = dist.pairwise(X)

with open("euclidean.txt", 'w') as f:
	f.write(np.array2string(distance, separator=', '))
	
distance = cosine_distances(X)
with open("cosine.txt", 'w') as f:
	f.write(np.array2string(distance, separator=', '))

##################################################################################


#Comparing between euclidean and cosine###########################################
distanceC = cosine_distances(X)
dist = DistanceMetric.get_metric('euclidean')
distanceE = dist.pairwise(X)
##################################################################################

#Plots############################################################################

label_colors = {'a': 'r', 'b': 'g', 'c': 'b', 'd': 'm'}


# Compute and plot first dendrogram.
fig = pylab.figure(figsize=(8,8))
ax1 = fig.add_axes([0.01,0.1,0.2,0.6])
Y = sch.average(sd.squareform(distanceC))
Z1 = sch.dendrogram(Y, orientation='left',labels=["G1", "G2", "G3", "G4", "G5"])
#ax1.set_xticks([])
#ax1.set_yticks([])
xlbls = ax1.get_xmajorticklabels()





# Compute and plot second dendrogram.
ax2 = fig.add_axes([0.25,0.74,0.6,0.2])
Y = sch.average(sd.squareform(distanceE))
Z2 = sch.dendrogram(Y,labels=["G1", "G2", "G3", "G4", "G5"])
#ax2.set_xticks([])
#ax2.set_yticks([])
xlbls = ax2.get_xmajorticklabels()


#Plot distance matrix.
axmatrix = fig.add_axes([0.25,0.1,0.6,0.6])
axmatrix2 = triamatrix(distanceC, axmatrix, rot=0, cmap="Blues")
idx1 = Z1['leaves']
idx2 = Z2['leaves']
distance = distanceE[idx1,:]
distance = distanceE[:,idx2]
im = axmatrix.matshow(distanceE, aspect='auto', origin='lower', cmap="Blues")
axmatrix.set_xticks([])
axmatrix.set_yticks([])

pylab.xticks(rotation=-90, fontsize=8)


axmatrix.set_xticklabels(idx1, minor=False)
axmatrix.xaxis.set_label_position('bottom')
axmatrix.xaxis.tick_bottom()

axmatrix.set_yticklabels(idx2, minor=False)
axmatrix.yaxis.set_label_position('right')
axmatrix.yaxis.tick_right()

#Plot colorbar.
axcolor = fig.add_axes([0.87,0.1,0.02,0.6])
axcolor2 = fig.add_axes([0.93,0.1,0.02,0.6])
pylab.colorbar(im, cax=axcolor)
pylab.colorbar(axmatrix2, cax=axcolor2)
fig.show()
fig.savefig('dend.pdf')



# ax = fig.add_axes([0.3,0.1,0.6,0.6])
# im1 = triamatrix(distanceC, ax, rot=0, cmap="Blues")
# im2 = triamatrix(distanceE, ax, rot=90, cmap="Reds")
# ax.set_xlim(-.5,distanceC.shape[1]-.5)
# ax.set_ylim(-.5,distanceE.shape[0]-.5)
# ax.set_xticks([])
# ax.set_yticks([])
# #ax.colorbar(im1, fig=fig, )
# #ax.colorbar(im2, fig=fig, )

# fig.show()
# fig.savefig('dendrogram.png')










# Plot distance matrix.
#axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
#idx1 = Z1['leaves']
#idx2 = Z2['leaves']
#distance = distance[idx1,:]
#distance = distance[:,idx2]
#im = axmatrix.matshow(distance, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
#axmatrix.set_xticks([])
#axmatrix.set_yticks([])

# Plot colorbar.
#axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
#pylab.colorbar(im, cax=axcolor)
#fig.show()
#fig.savefig('dendrogram.png')
##################################################################################

















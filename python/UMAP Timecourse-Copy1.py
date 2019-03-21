
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')
import os
import pylab
import numpy as np
import umap
import glob


# In[2]:


sns.set(context='paper', style='darkgrid', rc={'figure.facecolor':'white'}, font_scale=2)


# In[95]:


#get file names in directory that end in csv
os.chdir('/Users/s1249052/PhD/flow data/vac69a/t cells only/experiment_210618_files/csv_by_person/all')
os.getcwd()
filenames = glob.glob('*.csv')
print(filenames, len(filenames))


# In[96]:


big_boi = pd.DataFrame()

for n in filenames:
    file = pd.read_csv(n)
    file['id'] = '{}'.format(n)
    big_boi = big_boi.append([file], ignore_index=True, sort=False)


print(len(big_boi))
print(len(big_boi.columns))


# In[97]:


#print(big_boi.columns)
#print(len(big_boi.columns))
print(big_boi.iloc[:, 1:35].head(n=2))


# In[98]:


#######great
reducer = umap.UMAP(n_neighbors=20, min_dist=0.1, n_components=2, metric='euclidean')
get_ipython().run_line_magic('time', 'embedding = reducer.fit_transform(big_boi.iloc[:, 1:35])')
embedding.shape 
big_boi['umap1'] = embedding[:, 0]
big_boi['umap2'] = embedding[:,1]

channels=list(big_boi.columns.values)

for n in channels[0:35]:
    cmap=plt.get_cmap('nipy_spectral')
    fig=plt.figure(figsize=(15,10))
    ax=plt.subplot()
    ax.scatter(big_boi['umap1'], big_boi['umap2'], c=big_boi[n], cmap=cmap, s=1)
    plt.title('{}'.format(n))
    #plt.savefig('{}.pdf'.format(n))
    plt.show()


# In[63]:


# tried but didn't look nice:

#(n_neighbors=15, min_dist=0.1, n_components=2, metric='euclidean')
#(n_neighbors=15, min_dist=0.2, n_components=2, metric='euclidean')
#(n_neighbors=20, min_dist=0.2, n_components=2, metric='euclidean')

# maybe try increasing nearest neighbor iteratively? makes sense that it would have to go up as
# a result of the much higher density in features conserved between people & timepoints

# ask irina about algorithm that finds good metaparameters


# In[99]:


#print(np.unique(big_boi['id'])[0])
#maybs = big_boi.loc[big_boi['id'] == 'C+10_02_UMAP.csv']
#maybs.head(n=20)


# In[103]:


os.chdir('/Users/s1249052/PhD/flow data/vac69a/t cells only/experiment_210618_files/csv_by_person/all/big_one')
for source in np.unique(big_boi['id']):
    new_df = pd.DataFrame()
    new_df = big_boi.loc[big_boi['id'] == source]
    new_df = new_df.drop(columns = ['id'])
    new_df = new_df.drop(columns = ['Unnamed: 0'])
    new_df = new_df.drop(columns = ['Timepoint'])                                   
    new_df.to_csv('{}'.format(source), sep=',') 
    
os.getcwd()


# In[101]:


print(big_boi.columns)


# In[102]:


Baseline_csv = sub_concat[sub_concat['Timepoint'] == 'C-1']
C_10_csv = sub_concat[sub_concat['Timepoint'] == 'C+10']
C_12_csv = sub_concat[sub_concat['Timepoint'] == 'C+12']
DoD_csv = sub_concat[sub_concat['Timepoint'] == 'DoD']
T_6_csv = sub_concat[sub_concat['Timepoint'] == 'T+6']

Baseline_csv = Baseline_csv.drop(columns = ['Timepoint'])
C_10_csv = C_10_csv.drop(columns = ['Timepoint'])
C_12_csv = C_12_csv.drop(columns = ['Timepoint'])
DoD_csv = DoD_csv.drop(columns = ['Timepoint'])
T_6_csv = T_6_csv.drop(columns = ['Timepoint'])


# In[ ]:


Baseline_csv.to_csv('baseline_09_UMAP.csv', sep=',') 
C_10_csv.to_csv('C+10_09_UMAP.csv', sep=',') 
C_12_csv.to_csv('C+12_09_UMAP.csv', sep=',') 
DoD_csv.to_csv('DoD_09_UMAP.csv', sep=',') 
T_6_csv.to_csv('T+6_09_UMAP.csv', sep=',') 

os.getcwd()


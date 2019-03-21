
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


# In[23]:


os.chdir('/Users/s1249052/PhD/flow data/vac69a/t cells only/experiment_210618_files')
os.getcwd()

flow_data = pd.read_csv("vac69a_07_t+6.csv",sep=',')


# In[49]:


filenames


# In[48]:


filenames = glob.glob('*.csv')



# In[ ]:


#organise files into folders according to individual

for i in filenames:
    
    temp=pd.read_csv("*.csv")
    temp=flow_data[["In115Di","Pr141Di", "Nd142Di", "Nd143Di", "Nd144Di", "Nd145Di", "Nd146Di", "Nd148Di", "Sm149Di", "Nd150Di", "Eu151Di", "Eu153Di", "Sm154Di", "Gd155Di", "Gd156Di", "Gd158Di",
    "Tb159Di", "Gd160Di", "Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di", "Ho165Di", "Er166Di", "Er167Di", "Er168Di",
    "Tm169Di", "Yb171Di", "Yb172Di", "Yb173Di", "Yb174Di", "Lu175Di", "Pt198Di", "Bi209Di"]].copy()
    temp.columns.columns = ['115In_CD57', '141Pr_HLA-DR', '142Nd_BCL-2', '143Nd_CD45RA', '144Nd_GZB', '145Nd_CD4', '146Nd_Vd2',
    '148Nd_ICOS', '149Sm_CXCR5', '150Nd_CD95', '151Eu_CD103', '153Eu_Va7.2', '154Sm_TIM-3', '155Gd_PD1',
    '156Gd_CD161', '158Gd_CD27', '159Tb_FoxP3', '160Gd_CTLA4', '161Dy_Tbet', '162Dy_IntegrinB7', '163Dy_CD28', '164Dy_Ki-67',
    '165Ho_CD45RO', '166Er_CD56', '167Er_CCR7', '168Er_CD127', '169Tm_CD38', '171Yb_CD49d', '172Yb_CD25', '173Yb_CD39',
    '174Yb_CLA', '175Lu_Perforin', '198Pt_CD8', '209Bi_CD16']
    temp=temp.apply(np.arcsinh, 'columns')

    reducer = umap.UMAP(n_neighbors=30, min_dist=0.01, n_components=2, metric='euclidean')
    embedding = reducer.fit_transform(cytobank)
    embedding.shape 
    temp['umap1'] = embedding[:, 0]
    temp['umap2'] = embedding[:,1]
    cytobank.to_csv('temp.csv', sep=',')   
    
    
}


# In[24]:


#cytobank=flow_data[["Cd114Di", "In115Di","Pr141Di", "Nd142Di", "Nd143Di", "Nd144Di", "Nd145Di", "Nd146Di", "Sm147Di", "Nd148Di", "Sm149Di", "Nd150Di", "Eu151Di", "Eu153Di", "Sm154Di", "Gd155Di", "Gd156Di", "Gd158Di",
"Tb159Di", "Gd160Di", "Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di", "Ho165Di", "Er166Di", "Er167Di", "Er168Di",
"Tm169Di", "Yb171Di", "Yb172Di", "Yb173Di", "Yb174Di", "Lu175Di", "Yb176Di", "Pt198Di", "Bi209Di"]].copy()

#cytobank.columns = ['114Cd_CD14', '115In_CD57', '141Pr_HLA-DR', '142Nd_BCL-2', '143Nd_CD45RA', '144Nd_GZB', '145Nd_CD4', '146Nd_Vd2',
'147Sm_CD20', '148Nd_ICOS', '149Sm_CXCR5', '150Nd_CD95', '151Eu_CD103', '153Eu_Va7.2', '154Sm_TIM-3', '155Gd_PD1',
'156Gd_CD161', '158Gd_CD27', '159Tb_FoxP3', '160Gd_CTLA4', '161Dy_Tbet', '162Dy_IntegrinB7', '163Dy_CD28', '164Dy_Ki-67',
'165Ho_CD45RO', '166Er_CD56', '167Er_CCR7', '168Er_CD127', '169Tm_CD38', '171Yb_CD49d', '172Yb_CD25', '173Yb_CD39',
'174Yb_CLA', '175Lu_Perforin', '176Yb_CX3CR1', '198Pt_CD8', '209Bi_CD16']


# In[38]:


#copy with the channels that worked/are interesting

cytobank=flow_data[["In115Di","Pr141Di", "Nd142Di", "Nd143Di", "Nd144Di", "Nd145Di", "Nd146Di", "Nd148Di", "Sm149Di", "Nd150Di", "Eu151Di", "Eu153Di", "Sm154Di", "Gd155Di", "Gd156Di", "Gd158Di",
"Tb159Di", "Gd160Di", "Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di", "Ho165Di", "Er166Di", "Er167Di", "Er168Di",
"Tm169Di", "Yb171Di", "Yb172Di", "Yb173Di", "Yb174Di", "Lu175Di", "Pt198Di", "Bi209Di"]].copy()

cytobank.columns = ['115In_CD57', '141Pr_HLA-DR', '142Nd_BCL-2', '143Nd_CD45RA', '144Nd_GZB', '145Nd_CD4', '146Nd_Vd2',
'148Nd_ICOS', '149Sm_CXCR5', '150Nd_CD95', '151Eu_CD103', '153Eu_Va7.2', '154Sm_TIM-3', '155Gd_PD1',
'156Gd_CD161', '158Gd_CD27', '159Tb_FoxP3', '160Gd_CTLA4', '161Dy_Tbet', '162Dy_IntegrinB7', '163Dy_CD28', '164Dy_Ki-67',
'165Ho_CD45RO', '166Er_CD56', '167Er_CCR7', '168Er_CD127', '169Tm_CD38', '171Yb_CD49d', '172Yb_CD25', '173Yb_CD39',
'174Yb_CLA', '175Lu_Perforin', '198Pt_CD8', '209Bi_CD16']


# In[39]:


len(cytobank.columns)


# In[41]:


cytobank=cytobank.apply(np.arcsinh, 'columns')

reducer = umap.UMAP(n_neighbors=30, min_dist=0.01, n_components=2, metric='euclidean')

get_ipython().run_line_magic('time', 'embedding = reducer.fit_transform(cytobank)')
embedding.shape  


# In[42]:


cytobank['umap1'] = embedding[:, 0]
cytobank['umap2'] = embedding[:,1]


# In[43]:


channels=list(cytobank.columns.values)

for n in channels:
    cmap=plt.get_cmap('nipy_spectral')
    fig=plt.figure(figsize=(15,10))
    ax=plt.subplot()
    ax.scatter(cytobank['umap1'], cytobank['umap2'],c=cytobank[n], cmap=cmap, s=1)
    plt.title('{}'.format(n))
    #plt.savefig('{}.pdf'.format(n))
    plt.show()


# In[21]:


os.getcwd()


# In[44]:


cytobank.to_csv('vac69a_03_t+6_umap.csv', sep=',')


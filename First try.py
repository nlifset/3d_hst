
# coding: utf-8

# In[ ]:

ls


# In[ ]:

cd Downloads/


# In[ ]:

ls


# In[6]:

cd aegis_3dhst.v4.1.cats/


# In[7]:

ls


# In[8]:

cd Catalog/


# In[9]:

ls


# In[ ]:

from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt
import pylab


# In[11]:

data = ascii.read("aegis_3dhst.v4.1.cat")
data[0]


# In[12]:

data[1]


# In[13]:

cd


# In[14]:

ls


# In[15]:

cd Downloads/


# In[16]:

ls


# In[17]:

cd aegis_3dhst.v4.1.cats/


# In[18]:

ls


# In[19]:

cd Eazy/


# In[20]:

ls


# In[21]:

data2 = ascii.read("aegis_3dhst.v4.1.zout")


# In[30]:

data2[0]


# In[23]:

cd


# In[24]:

cd Downloads/


# In[25]:

ls


# In[26]:

cd cosmos_3dhst.v4.1.cats/


# In[27]:

ls


# In[28]:

cd Catalog/


# In[21]:

use_phot = data["use_phot"]


# In[22]:

idx, = np.where(use_phot == 1.0)


# In[23]:

len(data[idx])


# In[24]:

data_flagged = data[idx]


# In[25]:

ls


# In[26]:

cd Downloads/


# In[27]:

cd ..


# In[28]:

cd ..


# In[29]:

ls


# In[30]:

cd aegis_3dhst.v4.1.cats/


# In[31]:

ls


# In[32]:

cd Fast/


# In[33]:

ls


# In[34]:

data_fast = ascii.read("aegis_3dhst.v4.1.fout")
data_fast[0]


# In[35]:

data_fast_flagged = data_fast[idx]


# In[36]:

len(data_fast)
len(data_fast_flagged)


# In[37]:

len(data_fast)


# In[ ]:

lmass = data_fast_flagged["lmass"]


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




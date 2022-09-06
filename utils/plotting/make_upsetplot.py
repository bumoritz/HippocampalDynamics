# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 14:21:32 2022
@author: Moritz
"""

#from upsetplot import generate_counts
#from upsetplot import plot
#from matplotlib import pyplot


import numpy as np
import pandas as pd
import os
import upsetplot
from matplotlib import pyplot



path = r'C:\Users\Moritz\Documents\MATLAB\SniffinHippo\temp\upsetplot.npy'



data_npy = np.load(path,allow_pickle=True)
#numCells = np.size(data_npy,0)
#numPredictorGroups = np.size(data_npy,1)

data_bool = data_npy.astype(bool)
data_dict = {"OdourA": data_bool[:,0], "OdourX": data_bool[:,1], "OdourB": data_bool[:,2], "OdourY": data_bool[:,3], "Reward": data_bool[:,4]}

#data_upset = upsetplot.from_indicators(data_bool)
data_upset = upsetplot.from_indicators(data_dict)

upsetplot.UpSet(data_upset, sort_by='degree' , sort_categories_by=None, subset_size='sum', orientation='vertical', show_counts=True, show_percentages=True).plot()
#max_degree=2
pyplot.suptitle('Mixed selectivity')


pyplot.show(block=True)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:49:03 2021

@author: criph
"""

import pandas as pd

#carregamento da base de dados e separação de previsores da classe
base = pd.read_csv("census.csv")
previsores = base.iloc[: ,0:14].values
classe = base.iloc[:,14].values

#Transformação de dados categoricos para numéricos:
from sklearn.preprocessing import LabelEncoder
labelencoder_previsores = LabelEncoder()
#labels = labelencoder_previsores.fit_transform(previsores[:,1])

for i in [1,3,5,6,7,8,9,13] :
      previsores[:,i] = labelencoder_previsores.fit_transform(previsores[:,i])


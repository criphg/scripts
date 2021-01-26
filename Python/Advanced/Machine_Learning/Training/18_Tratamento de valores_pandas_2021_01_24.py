#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 16:35:34 2021

@author: criph
"""

import pandas as pd
base = pd.read_csv('credit_data.csv')
#print(base.describe())

#base.loc[base["age"] < 0]

#apagar a coluna - registros inteiros
#base.drop('age',1, inplace=True)

#apagar apenas os registros com problemas
#base.drop(base[base.age <0].index, inplace = True)

#mudar o nome dos registros
#preencher com a média
#base.mean()
#base['age'][base.age > 0].mean()
#base.loc[base.age < 0, 'age'] = 40.92

#print(base.age[base['age'] < 0])

pd.isnull(base.age)
base[pd.isnull(base.age)]
base.loc[pd.isnull(base.age)]

#separação dos previsores e classe, id não entra pois não fornece dados para classificação
previsores = base.iloc[:, 1:4]
classe = base.iloc[:,4]



#preste atenção no perigo da generalização: ele preenche os dados com a média, mas leva em consideração os valores de idade negativo por exemplo. A média de idade preenchida deveria ser 40.92 e náo 40.8.
#Além disso pela nova versão ou devido a caracteristicas do python3 não pode haver denominação no número de linhas na chamada do dataframe, a não ser que use '.iloc()' antes do objeto dataframe.
#print(previsore[: , 0:3]) #não funciona!!!!
#por isso o código sugerido pelo curso não funciona.
from sklearn.impute import SimpleImputer
imputer = SimpleImputer()
imputer = imputer.fit(previsores)
previsores2 = imputer.transform(previsores)

print(previsores2)

#print(previsores.iloc[0:3, 0:1])
#print(previsores)
#type(previsores)

#Faltou o escalonamento. Será adicionado.
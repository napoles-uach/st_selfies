import streamlit as st
import pandas as pd
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

X0={'[F]':['F',1],'[=O]':['O',2],'[#N]':['N',3] ,'[O]':['O',2] ,'[N]':['N',3] ,'[=N]':['N',3] ,'[C]':['C',4] ,'[=C]':['C',4] ,'[#C]':['C',4] }
X1={'[F]':['F',0],'[=O]':['O',0],'[#N]':['N',0] ,'[O]':['O',1] ,'[N]':['N',2] ,'[=N]':['N',2] ,'[C]':['C',3] ,'[=C]':['C',3] ,'[#C]':['C',3] }
X2={'[F]':['F',0],'[=O]':['=O',0],'[#N]':['=N',0] ,'[O]':['O',1] ,'[N]':['N',2] ,'[=N]':['=N',2] ,'[C]':['C',3] ,'[=C]':['=C',2] ,'[#C]':['=C',2] }
X3={'[F]':['F',0],'[=O]':['=O',0],'[#N]':['#N',0] ,'[O]':['O',1] ,'[N]':['N',2] ,'[=N]':['=N',1] ,'[C]':['C',3] ,'[=C]':['=C',2] ,'[#C]':['#C',1] }
X4={'[F]':['F',0],'[=O]':['=O',0],'[#N]':['#N',0] ,'[O]':['O',1] ,'[N]':['N',2] ,'[=N]':['=N',1] ,'[C]':['C',3] ,'[=C]':['=C',2] ,'[#C]':['#C',1] }

X=[X0,X1,X2,X3,X4]

chain=''


def braketoff(lista):
  new_list=[]
  for el in lista:
    new=el.replace('[','')
    new=new.replace(']','')
    new_list.append(new)

  return new_list


i=1
key=0

SELFIES=[]
dic=X[0]


init=st.sidebar.selectbox('token:'+str(key),dic,key=str(key))

SELFIES.append(init)


i=1
while i!=0:
    key+=1


    dic=X[i]    
    init=st.sidebar.selectbox('token:'+str(key),dic,key=str(key))

    SELFIES.append(init)
    i=int(dic[init][1])

    if i==0:
      final=braketoff(SELFIES)
      st.title('SELFIES 😎 (BETA)')
      st.write(SELFIES)
      for f in final:
        chain+='['+f+']'

      #st.title(chain)

      compound_smiles=sf.decoder(chain) 
      
      m = Chem.MolFromSmiles(compound_smiles)

      Draw.MolToFile(m,'mol.png')
      st.image('mol.png')
  

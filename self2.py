import streamlit as st
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# Diccionarios de tokens válidos según la valencia restante
X = [
    {'[F]': ['F', 1], '[=O]': ['O', 2], '[#N]': ['N', 3], '[O]': ['O', 2],
     '[N]': ['N', 3], '[=N]': ['N', 3], '[C]': ['C', 4], '[=C]': ['C', 4], '[#C]': ['C', 4]},
    
    {'[F]': ['F', 0], '[=O]': ['O', 0], '[#N]': ['N', 0], '[O]': ['O', 1],
     '[N]': ['N', 2], '[=N]': ['N', 2], '[C]': ['C', 3], '[=C]': ['C', 3], '[#C]': ['C', 3]},
    
    {'[F]': ['F', 0], '[=O]': ['=O', 0], '[#N]': ['=N', 0], '[O]': ['O', 1],
     '[N]': ['N', 2], '[=N]': ['=N', 2], '[C]': ['C', 3], '[=C]': ['=C', 2], '[#C]': ['=C', 2]},
    
    {'[F]': ['F', 0], '[=O]': ['=O', 0], '[#N]': ['#N', 0], '[O]': ['O', 1],
     '[N]': ['N', 2], '[=N]': ['=N', 1], '[C]': ['C', 3], '[=C]': ['=C', 2], '[#C]': ['#C', 1]},
    
    {'[F]': ['F', 0], '[=O]': ['=O', 0], '[#N]': ['#N', 0], '[O]': ['O', 1],
     '[N]': ['N', 2], '[=N]': ['=N', 1], '[C]': ['C', 3], '[=C]': ['=C', 2], '[#C]': ['#C', 1]}
]

# Función para quitar corchetes de una lista de tokens
def remove_brackets(tokens):
    return [token.replace('[', '').replace(']', '') for token in tokens]

# Inicialización
chain = ''
i = 1
key = 0
SELFIES = []

# Primer token desde estado inicial
dic = X[0]
selected = st.sidebar.selectbox(f'token: {key}', dic, key=str(key))
SELFIES.append(selected)

# Construcción de la cadena SELFIES
while i != 0:
    key += 1
    dic = X[i]
    selected = st.sidebar.selectbox(f'token: {key}', dic, key=str(key))
    SELFIES.append(selected)
    i = int(dic[selected][1])

    if i == 0:
        # Se termina la cadena, convertir a SMILES y mostrar la molécula
        final_tokens = remove_brackets(SELFIES)
        st.title('st_SELFIES 😎 (BETA)')
        st.write(SELFIES)

        chain = ''.join(f'[{tok}]' for tok in final_tokens)
        smiles = sf.decoder(chain)
        mol = Chem.MolFromSmiles(smiles)

        if mol:
            img = Draw.MolToImage(mol)
            st.image(img)
        else:
            st.error("No se pudo generar una molécula válida.")


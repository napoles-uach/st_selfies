import streamlit as st
import pandas as pd
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# Cargar archivo CSV con la tabla de derivación
df = pd.read_csv('selfies_table.csv')

# Construir diccionarios por estado
X_dict = {
    state: {
        row['symbol']: [row['smiles'], int(row['next_state'])]
        for _, row in group.iterrows()
    }
    for state, group in df.groupby('state')
}

# Lista de estados en orden
states = ['X0', 'X1', 'X2', 'X3', 'X4']

# Función para quitar corchetes
def remove_brackets(tokens):
    return [tok.replace('[', '').replace(']', '') for tok in tokens]

# Inicialización
SELFIES = []
chain = ''
key = 0
state_index = 0

# Primer token desde X0
dic = X_dict[states[state_index]]
selected = st.sidebar.selectbox(f'token: {key}', dic, key=str(key))
SELFIES.append(selected)

# Construcción de cadena SELFIES
next_state = int(dic[selected][1])

while next_state != 0:
    key += 1
    dic = X_dict[states[next_state]]
    selected = st.sidebar.selectbox(f'token: {key}', dic, key=str(key))
    SELFIES.append(selected)
    next_state = int(dic[selected][1])

# Mostrar resultados
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

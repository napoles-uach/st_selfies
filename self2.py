import streamlit as st
import pandas as pd
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# Cargar el archivo corregido
df = pd.read_csv('selfies_table_corrected.csv')

# Construir diccionario FSM
X_dict = {
    state: {
        row['symbol']: [row['smiles'], row['next_state']]
        for _, row in group.iterrows()
    }
    for state, group in df.groupby('state')
}

states = ['X0', 'X1', 'X2', 'X3', 'X4']

def remove_brackets(tokens):
    return [tok.replace('[', '').replace(']', '') for tok in tokens]

SELFIES = []
chain = ''
key = 0
state_index = 0

# Primer token desde X0
dic = X_dict[states[state_index]]
selected = st.sidebar.selectbox(f'token: {key}', dic, key=str(key))
SELFIES.append(selected)
next_state = dic[selected][1]

while next_state != 0:
    key += 1
    dic = X_dict[states[next_state]]
    selected = st.sidebar.selectbox(f'token: {key}', dic, key=str(key))
    SELFIES.append(selected)

    # 🔀 Manejo de ramas
    if 'Branch' in selected:
        Q = int(selected.replace('[Branch', '').replace(']', ''))
        branch_state = states[next_state]
        st.sidebar.markdown(f'**→ Dentro de rama de longitud {Q}**')
        for q in range(Q):
            key += 1
            branch_dic = X_dict[branch_state]
            branch_token = st.sidebar.selectbox(f'branch_token: {key}', branch_dic, key=str(key))
            SELFIES.append(branch_token)
            branch_state = states[branch_dic[branch_token][1]]
        next_state = dic[selected][1]

    # 🔁 Manejo de anillos
    elif 'Ring' in selected:
        # Agrega el símbolo sin cambiar de estado
        next_state = dic[selected][1]  # permanece en el mismo estado

    else:
        next_state = dic[selected][1]

# Visualización
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

import streamlit as st
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# Diccionarios de estados (representan valencias restantes y transiciones)
X0 = {'[F]': ['F', 1], '[=O]': ['O', 2], '[#N]': ['N', 3], '[O]': ['O', 2], '[N]': ['N', 3], '[=N]': ['N', 3], '[C]': ['C', 4], '[=C]': ['C', 4], '[#C]': ['C', 4]}
X1 = {'[F]': ['F', 0], '[=O]': ['O', 0], '[#N]': ['N', 0], '[O]': ['O', 1], '[N]': ['N', 2], '[=N]': ['N', 2], '[C]': ['C', 3], '[=C]': ['C', 3], '[#C]': ['C', 3]}
X2 = {'[F]': ['F', 0], '[=O]': ['=O', 0], '[#N]': ['=N', 0], '[O]': ['O', 1], '[N]': ['N', 2], '[=N]': ['=N', 2], '[C]': ['C', 3], '[=C]': ['=C', 2], '[#C]': ['=C', 2]}
X3 = {'[F]': ['F', 0], '[=O]': ['=O', 0], '[#N]': ['#N', 0], '[O]': ['O', 1], '[N]': ['N', 2], '[=N]': ['=N', 1], '[C]': ['C', 3], '[=C]': ['=C', 2], '[#C]': ['#C', 1]}
X4 = {'[F]': ['F', 0], '[=O]': ['=O', 0], '[#N]': ['#N', 0], '[O]': ['O', 1], '[N]': ['N', 2], '[=N]': ['=N', 1], '[C]': ['C', 3], '[=C]': ['=C', 2], '[#C]': ['#C', 1]}
X = [X0, X1, X2, X3, X4]

# Inicialización del estado
if 'tokens' not in st.session_state:
    st.session_state.tokens = []
if 'current_state' not in st.session_state:
    st.session_state.current_state = 0
if 'step' not in st.session_state:
    st.session_state.step = 0

# Función para limpiar corchetes
def braketoff(lista):
    return [el.replace('[', '').replace(']', '') for el in lista]

st.sidebar.title('Constructor SELFIES')
dic = X[st.session_state.current_state]

# Selección de token actual
selected_token = st.sidebar.selectbox(f'Token {st.session_state.step}', dic, key=str(st.session_state.step))

# Botón para agregar token
if st.sidebar.button('Agregar token'):
    st.session_state.tokens.append(selected_token)
    st.session_state.current_state = int(dic[selected_token][1])
    st.session_state.step += 1

# Botón para reiniciar
if st.sidebar.button('Reiniciar'):
    st.session_state.tokens = []
    st.session_state.current_state = 0
    st.session_state.step = 0

# Mostrar resultados cuando el estado llega a 0
if st.session_state.current_state == 0 and st.session_state.tokens:
    st.title('st_SELFIES 😎 (BETA)')
    st.write("SELFIES tokens:", st.session_state.tokens)
    
    selfies_string = ''.join(['[' + token.replace('[', '').replace(']', '') + ']' for token in st.session_state.tokens])
    st.write("SELFIES:", selfies_string)
    
    smiles = sf.decoder(selfies_string)
    st.write("SMILES:", smiles)

    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
        st.image(img)
    else:
        st.error("Error al generar la molécula.")

import streamlit as st
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# Diccionarios de transición por estado de valencia
X0 = {'[F]': ['F', 1], '[=O]': ['O', 2], '[#N]': ['N', 3], '[O]': ['O', 2],
      '[N]': ['N', 3], '[=N]': ['N', 3], '[C]': ['C', 4], '[=C]': ['C', 4], '[#C]': ['C', 4]}
X1 = {'[F]': ['F', 0], '[=O]': ['O', 0], '[#N]': ['N', 0], '[O]': ['O', 1],
      '[N]': ['N', 2], '[=N]': ['N', 2], '[C]': ['C', 3], '[=C]': ['C', 3], '[#C]': ['C', 3]}
X2 = {'[F]': ['F', 0], '[=O]': ['=O', 0], '[#N]': ['=N', 0], '[O]': ['O', 1],
      '[N]': ['N', 2], '[=N]': ['=N', 2], '[C]': ['C', 3], '[=C]': ['=C', 2], '[#C]': ['=C', 2]}
X3 = {'[F]': ['F', 0], '[=O]': ['=O', 0], '[#N]': ['#N', 0], '[O]': ['O', 1],
      '[N]': ['N', 2], '[=N]': ['=N', 1], '[C]': ['C', 3], '[=C]': ['=C', 2], '[#C]': ['#C', 1]}
X4 = {'[F]': ['F', 0], '[=O]': ['=O', 0], '[#N]': ['#N', 0], '[O]': ['O', 1],
      '[N]': ['N', 2], '[=N]': ['=N', 1], '[C]': ['C', 3], '[=C]': ['=C', 2], '[#C]': ['#C', 1]}

X = [X0, X1, X2, X3, X4]

# Función para quitar corchetes
def braketoff(lista):
    return [el.replace('[', '').replace(']', '') for el in lista]

# Encabezado de la app
st.title('st_SELFIES 😎 (BETA)')

# Inicialización de estado
if 'SELFIES' not in st.session_state:
    st.session_state.SELFIES = []
    st.session_state.state = 0
    st.session_state.token_count = 0

# Mostrar el selector de token correspondiente al estado actual
dic = X[st.session_state.state]
token = st.sidebar.selectbox(f'Selecciona token ({st.session_state.token_count})', options=list(dic.keys()), key=f"token_{st.session_state.token_count}")

# Botón para agregar token
if st.sidebar.button("Agregar token"):
    st.session_state.SELFIES.append(token)
    st.session_state.state = int(dic[token][1])
    st.session_state.token_count += 1

# Mostrar SELFIES parcial
st.subheader("Cadena SELFIES actual")
st.write(st.session_state.SELFIES)

# Si se llegó al estado terminal (0), construir molécula
if st.session_state.state == 0 and st.session_state.SELFIES:
    final = braketoff(st.session_state.SELFIES)
    chain = ''.join(f'[{f}]' for f in final)

    st.subheader("SELFIES final:")
    st.code(chain)

    # Decodificación y visualización
    smiles = sf.decoder(chain)
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        st.subheader("Molécula generada:")
        img = Draw.MolToImage(mol)
        st.image(img)
    else:
        st.error("No se pudo generar la molécula desde el SMILES generado.")

# Botón para reiniciar todo
if st.sidebar.button("Reiniciar"):
    st.session_state.SELFIES = []
    st.session_state.state = 0
    st.session_state.token_count = 0

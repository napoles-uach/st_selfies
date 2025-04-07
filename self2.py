import streamlit as st
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# --- Diccionarios de tokens por estado (valencia restante simulada)
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

# --- Inicializar estado
if "tokens" not in st.session_state:
    st.session_state.tokens = []
if "state" not in st.session_state:
    st.session_state.state = 0

st.title("🔬 Constructor de SELFIES interactivo")

# --- Diccionario actual según estado
current_dict = X[st.session_state.state]

# --- Selección de token
selected = st.selectbox("Selecciona el siguiente token", current_dict.keys())

# --- Botón para agregar token
if st.button("Agregar token"):
    st.session_state.tokens.append(selected)
    st.session_state.state = int(current_dict[selected][1])

# --- Mostrar SELFIES parcial
if st.session_state.tokens:
    st.subheader("SELFIES actual")
    st.write(st.session_state.tokens)

    # Armar cadena completa
    chain = "".join(st.session_state.tokens)

    # Mostrar SMILES y molécula si ya se cerró (estado 0)
    if st.session_state.state == 0:
        smiles = sf.decoder(chain)
        st.success(f"SMILES: `{smiles}`")

        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))
            st.image(img)
        else:
            st.error("No se pudo generar la molécula. 😓")

# --- Botón para reiniciar
if st.button("Reiniciar"):
    st.session_state.tokens = []
    st.session_state.state = 0

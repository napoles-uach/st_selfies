import streamlit as st
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# Agregamos enlaces dobles y triples y ramificaciones
extra_tokens = {
    '[Branch1_1]': ['B', 0],
    '[Branch2_2]': ['B', 0],
    '[Ring1]': ['R', 0],
    '[Ring2]': ['R', 0]
}

# Diccionarios base de valencia (puedes enriquecerlos después)
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

# Agregar extra_tokens a todos los diccionarios
for xi in X:
    xi.update(extra_tokens)

# --- Estado inicial
if "tokens" not in st.session_state:
    st.session_state.tokens = []
if "state" not in st.session_state:
    st.session_state.state = 0

st.title("🧬 SELFIES Builder - Versión extendida")

# Diccionario actual
current_dict = X[st.session_state.state]

# Selección de token
selected = st.selectbox("Selecciona el siguiente token", current_dict.keys())

# Agregar token
if st.button("Agregar token"):
    st.session_state.tokens.append(selected)
    st.session_state.state = int(current_dict[selected][1])

# Mostrar SELFIES actual
if st.session_state.tokens:
    st.subheader("📋 Cadena SELFIES actual")
    st.write(" ".join(st.session_state.tokens))

    chain = "".join(st.session_state.tokens)

    # Si estado = 0, mol cerrada
    if st.session_state.state == 0:
        smiles = sf.decoder(chain)
        st.success(f"SMILES: `{smiles}`")
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))
            st.image(img)
        else:
            st.error("Error al generar la molécula")

# Reiniciar
if st.button("🔄 Reiniciar construcción"):
    st.session_state.tokens = []
    st.session_state.state = 0


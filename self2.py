import streamlit as st
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

st.set_page_config(page_title="SELFIES Builder", layout="centered")

st.title("☕ SELFIES Builder – Versión extendida estilo FSM")

# --- Tabla de valencias extendida ---
valencias = {
    '[F]': 1, '[Cl]': 1, '[Br]': 1, '[I]': 1,
    '[O]': 2, '[=O]': 2,
    '[N]': 3, '[=N]': 3, '[#N]': 3, '[n]': 3,
    '[C]': 4, '[=C]': 4, '[#C]': 4, '[c]': 3,
    '[S]': 2, '[=S]': 2,
    '[P]': 3, '[=P]': 3,
    '[B]': 3, '[=B]': 3,
    '[Branch1_1]': 0, '[Branch2_2]': 0,
    '[Ring1]': 0, '[Ring2]': 0, '[Ring3]': 0,
    '[Expl=Ring1]': 0, '[Expl=Ring2]': 0,
    '[Expl#Ring1]': 0, '[Expl#Ring2]': 0,
    '[nop]': 0
}

# --- Construcción de estados X0 a X4 ---
X = []
for estado in range(5):
    tabla_estado = {}
    for token, valencia in valencias.items():
        nueva_valencia = max(estado + valencia - 4, 0)
        tabla_estado[token] = [token.strip("[]"), nueva_valencia]
    X.append(tabla_estado)

# --- Inicialización de sesión ---
if "tokens" not in st.session_state:
    st.session_state.tokens = []
if "state" not in st.session_state:
    st.session_state.state = 0

# --- Diccionario actual según estado ---
estado_actual = st.session_state.state
diccionario_actual = X[estado_actual]

# No permitir ramas ni ciclos como primer token
if len(st.session_state.tokens) == 0:
    diccionario_actual = {k: v for k, v in diccionario_actual.items() if 'Branch' not in k and 'Ring' not in k}

# --- Selección de token ---
selected = st.selectbox("Selecciona el siguiente token:", sorted(diccionario_actual.keys()))

# --- Agregar token ---
if st.button("Agregar token"):
    st.session_state.tokens.append(selected)
    st.session_state.state = int(diccionario_actual[selected][1])

# --- Mostrar SELFIES actual ---
if st.session_state.tokens:
    st.subheader("📋 Cadena SELFIES actual:")
    chain = "".join(st.session_state.tokens)
    st.code(chain)

    # Intentar decodificar y mostrar estructura
    if st.session_state.state == 0:
        try:
            smiles = sf.decoder(chain)
            st.success(f"SMILES generado: `{smiles}`")
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))
                st.image(img, caption="Estructura generada")
            else:
                st.warning("RDKit no pudo generar una molécula válida.")
        except Exception as e:
            st.error(f"Error al decodificar: {e}")

# --- Botón para reiniciar ---
if st.button("🔄 Reiniciar construcción"):
    st.session_state.tokens = []
    st.session_state.state = 0

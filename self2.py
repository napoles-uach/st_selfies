import streamlit as st
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# Obtener todos los tokens válidos de SELFIES
tokens = sf.get_semantic_robust_alphabet()

# Estimar valencias por token (muy simplificado)
valencia_por_token = {
    '[C]': 4, '[=C]': 4, '[#C]': 4, '[c]': 3,
    '[N]': 3, '[=N]': 3, '[#N]': 3, '[n]': 3,
    '[O]': 2, '[=O]': 2,
    '[F]': 1, '[Cl]': 1, '[Br]': 1, '[I]': 1,
    '[S]': 2, '[=S]': 2, '[P]': 3, '[=P]': 3,
    '[B]': 3, '[=B]': 3,
    '[Branch1_1]': 0, '[Branch2_2]': 0,
    '[Ring1]': 0, '[Ring2]': 0, '[Ring3]': 0,
    '[Expl=Ring1]': 0, '[Expl=Ring2]': 0,
    '[Expl#Ring1]': 0, '[Expl#Ring2]': 0,
    '[nop]': 0, '[=nop]': 0  # no-op tokens
}

# Generar estados X0 a X4 basados en valencia restante
X = [{} for _ in range(5)]
for token in tokens:
    val = valencia_por_token.get(token, 1)  # si no se reconoce, se le da 1
    for state in range(5):
        new_state = max(state + val - 4, 0)
        X[state][token] = [token.strip('[]'), new_state]

# Inicializar estado
if "tokens" not in st.session_state:
    st.session_state.tokens = []
if "state" not in st.session_state:
    st.session_state.state = 0

st.title("🧬 SELFIES Builder - Compatible con moléculas complejas")

# Diccionario según estado actual
current_dict = X[st.session_state.state]

# Si estamos en el primer paso, filtramos ciclos y ramas
if len(st.session_state.tokens) == 0:
    current_dict = {k: v for k, v in current_dict.items() if 'Branch' not in k and 'Ring' not in k}

# Selección de token
selected = st.selectbox("Selecciona el siguiente token", current_dict.keys())

# Agregar token
if st.button("Agregar token"):
    st.session_state.tokens.append(selected)
    st.session_state.state = int(current_dict[selected][1])

# Mostrar SELFIES actual
if st.session_state.tokens:
    st.subheader("📋 SELFIES actual")
    st.code("".join(st.session_state.tokens))

    chain = "".join(st.session_state.tokens)

    # Mostrar molécula si el estado es 0
    if st.session_state.state == 0:
        try:
            smiles = sf.decoder(chain)
            st.success(f"SMILES generado: `{smiles}`")
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))
                st.image(img)
        except:
            st.error("❌ Error al decodificar la cadena SELFIES.")

# Botón de reinicio
if st.button("🔄 Reiniciar"):
    st.session_state.tokens = []
    st.session_state.state = 0

import streamlit as st
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# 1. Generar el alfabeto SELFIES extendido a partir de ejemplos aromáticos
ejemplos_smiles = [
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",  # Cafeína
    "c1ccncc1",                    # Piridina
    "c1ccccc1",                    # Benceno
    "C1=CC=CN=C1",                 # Piridina isómera
]
tokens = sf.get_alphabet_from_selfies([sf.encoder(s) for s in ejemplos_smiles])

# 2. Estimar valencias para los tokens
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
    '[nop]': 0, '[=nop]': 0
}

# 3. Construir diccionarios de estado (X0 a X4)
X = [{} for _ in range(5)]
for token in tokens:
    val = valencia_por_token.get(token, 1)  # Asumimos 1 si no está en el diccionario
    for state in range(5):
        new_state = max(state + val - 4, 0)
        X[state][token] = [token.strip('[]'), new_state]

# 4. Inicializar estado de Streamlit
if "tokens" not in st.session_state:
    st.session_state.tokens = []
if "state" not in st.session_state:
    st.session_state.state = 0

st.title("☕ SELFIES Builder extendido – ¡Listo para cafeína!")

# Diccionario del estado actual
current_dict = X[st.session_state.state]

# Al iniciar, no permitir ciclos ni ramas
if len(st.session_state.tokens) == 0:
    current_dict = {k: v for k, v in current_dict.items() if 'Branch' not in k and 'Ring' not in k}

# Selección de token
selected = st.selectbox("Selecciona el siguiente token", sorted(current_dict.keys()))

# Botón para agregar token
if st.button("Agregar token"):
    st.session_state.tokens.append(selected)
    st.session_state.state = int(current_dict[selected][1])

# Mostrar SELFIES actual
if st.session_state.tokens:
    st.subheader("📋 SELFIES actual")
    chain = "".join(st.session_state.tokens)
    st.code(chain)

    # Mostrar SMILES y estructura si el estado es 0
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
if st.button("🔄 Reiniciar construcción"):
    st.session_state.tokens = []
    st.session_state.state = 0

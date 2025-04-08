import streamlit as st
import pandas as pd
import requests
import os
import selfies as sf
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# --- Config ---
HF_API_TOKEN = os.environ.get("HF_API_TOKEN", "YOUR_HF_API_TOKEN_HERE")
MODEL_NAME = "mistralai/Mistral-7B-Instruct-v0.3"

# --- Load SELFIES transition table ---
df = pd.read_csv("selfies_table_full_extended.csv")

# --- Build FSM dict ---
X_dict = {
    state: {
        row['symbol']: [row['smiles'], int(row['next_state'])]
        for _, row in group.iterrows()
    }
    for state, group in df.groupby('state')
}

states = ['X0', 'X1', 'X2', 'X3', 'X4']

# --- Prompt builder ---
def build_selfies_prompt(current_state, options, current_chain, objective=None):
    base = f"""
Estás construyendo una molécula usando SELFIES. Cada token depende del estado actual de valencia.

Estado actual: {current_state}
SELFIES parcial: {' '.join(current_chain)}

Tokens válidos desde este estado:
{list(options.keys())}
"""
    if objective:
        base += f"\nTu objetivo es: {objective}\n"

    base += "\nElige el siguiente token SELFIES exacto. Devuelve solo el token, como por ejemplo: [C], [=O], [Branch1]."
    return base

# --- Hugging Face API ---
def query_hf_model(prompt):
    url = f"https://api-inference.huggingface.co/models/{MODEL_NAME}"
    headers = {"Authorization": f"Bearer {HF_API_TOKEN}"}
    payload = {"inputs": prompt}
    response = requests.post(url, headers=headers, json=payload)
    response.raise_for_status()
    return response.json()[0]['generated_text']

# --- App ---
st.title("🤖 SELFIES Builder con Agente LLM")

objective = st.text_input("🎯 Objetivo químico del agente (opcional)", "incluir un grupo carbonilo [=O] y un nitrógeno [N]")

if st.button("Construir molécula automáticamente"):
    state = 'X0'
    selfies_chain = []

    while True:
        options = X_dict[state]
        prompt = build_selfies_prompt(state, options, selfies_chain, objective)
        try:
            response_text = query_hf_model(prompt)
        except Exception as e:
            st.error(f"Error al consultar Hugging Face API: {e}")
            break

        # Mostrar lo que "piensa" el modelo
        with st.expander(f"🤖 Pensamiento del modelo en {state}"):
            st.markdown("**Prompt enviado:**")
            st.code(prompt)
            st.markdown("**Respuesta completa:**")
            st.code(response_text)

        # Intentar extraer un token válido
        selected = None
        for token in options:
            if token in response_text:
                selected = token
                break

        if not selected:
            st.warning("No se encontró un token válido en la respuesta. Finalizando.")
            break

        selfies_chain.append(selected)
        next_state = options[selected][1]
        state = f"X{next_state}"

        if next_state == 0:
            break

    st.subheader("Cadena SELFIES construida:")
    st.code(" ".join(selfies_chain))

    # Decodificar y mostrar
    final_string = ''.join(selfies_chain)
    smiles = sf.decoder(final_string)
    st.write(f"**SMILES:** {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
        st.image(img)
    else:
        st.error("No se pudo generar la molécula a partir del SMILES.")

import streamlit as st

st.title('Protein Variant Library Generator')

option = st.selectbox(
    'Choose an organism:',
    ('Human', 'E-Coli')
)

uploaded_file = st.file_uploader("Upload your fasta file here:", type=['fasta'])

temp_range = st.slider('DNA Melting Temperature', min_value=50, max_value=100, value=(60, 90))
oligo_length_range = st.slider('Oligonucleotide Length', min_value=30, max_value=2000, value=(60, 90))
cost_per_base = st.number_input('Cost per Base', min_value=0.01, max_value=1.0, value=0.1, step=0.01)
st.button('Submit')

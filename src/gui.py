import streamlit as st

st.title('Protein Variant Library Generator')

option = st.selectbox(
    'Choose an organism:',
    ('Human', 'E-Coli')
)

uploaded_file = st.file_uploader("Upload your fasta file here:", type=['fasta'])

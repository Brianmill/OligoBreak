import streamlit as st

st.title('Protein Variant Library Generator')

option = st.selectbox(
    'Choose an organism:',
    ('Human', 'E-Coli')
)

uploaded_file = st.file_uploader("Upload your fasta file here:", type=['fasta'])

temp_range = st.slider('Input temperature range', min_value=50, max_value=100, value=(60, 90))
oligo_length_range = st.slider('Input oligo length range', min_value=30, max_value=2000, value=(60, 90))

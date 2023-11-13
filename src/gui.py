import streamlit as st

st.title('Protein Variant Library Generator')

with st.form(key='my_form'):
    organism = st.selectbox(
        'Choose an organism:',
        ('Human', 'E-Coli')
    )
    fasta_file = st.file_uploader("Upload your fasta file here:", type=['fasta'])
    degenerate_base_file = st.file_uploader("Upload your degenerate base file here:", type=['fasta'])
    temp_range = st.slider('DNA Melting Temperature', min_value=50, max_value=100, value=(60, 90))
    oligo_length_range = st.slider('Oligonucleotide Length', min_value=30, max_value=2000, value=(60, 1250))
    cost_per_base = st.number_input('Cost per Base', min_value=0.01, max_value=1.0, value=0.1, step=0.01)
    submit = st.form_submit_button(label='Submit')

if submit:
    successful_submit = True
    if fasta_file is None:
        st.warning("Please upload a fasta file!")
        successful_submit = False
    if degenerate_base_file is None:
        st.warning("Please upload a degenerate base file!")
        successful_submit = False

    if successful_submit:
        st.success('File successfully uploaded and processed.')

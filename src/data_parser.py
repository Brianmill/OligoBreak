def import_organism_decodon_set_files():
    organism_decodon_set_files = Path(__file__).parent / 'data' / 'organism_AA_decodon_sets'
    raise NotImplementedError('Importation of organism decodon set files not yet implemented')


def extract_protein_sequence(protein_sequence_file: str, multiple_amino_file: str):
    raise (NotImplementedError('extract_protein_sequence not yet implemented'))


def store_protein_sequence(protein_sequence_data) -> str:
    raise (NotImplementedError('store_protein_sequence not yet implemented'))


def extract_organism_codon_frequencies(codon_frequencies_file: str):

    organism_file = open(codon_frequencies_file, "r")

    # Adjust this part according to your file format
    # This example assumes each line in the file is a key-value pair separated by a comma
    data = {}
    for line in file:
        key, value = line.strip().split(',')
        data[key] = value

    return data
    raise (NotImplementedError('extract_organism_codon_frequencies not yet implemented'))


def store_organism_codon_frequencies() -> str:
    raise (NotImplementedError('store_organism_codon_frequencies not yet implemented'))


def extract_degenerate_codons_key(degenerate_codons_key_file: str):
    raise (NotImplementedError('extract_organism_codon_frequencies not yet implemented'))


def store_degenerate_codons_key() -> str:
    raise (NotImplementedError('store_organism_codon_frequencies not yet implemented'))



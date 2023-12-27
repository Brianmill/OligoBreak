def import_organism_decodon_set_files():
    organism_decodon_set_files = Path(__file__).parent / 'data' / 'organism_AA_decodon_sets'
    raise NotImplementedError('Importation of organism decodon set files not yet implemented')


def parse_organism_decodon_set_file(file):

    organism_file = open(file, "r")

    # Adjust this part according to your file format
    # This example assumes each line in the file is a key-value pair separated by a comma
    data = {}
    for line in file:
        key, value = line.strip().split(',')
        data[key] = value

    return data
    raise (NotImplementedError('Importation of organism decodon set file not yet implemented'))

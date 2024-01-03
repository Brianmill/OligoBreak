import unittest
import src

TEST_FILES = Path(__file__).parent / 'test_files'
class TestDataParser(unittest.TestCase):

    def test_temp(self):
        result = 15
        self.assertEqual(result, 15)

    def test_protein_sequence_extraction(self):
        protein_sequence_file = str(TEST_FILES / "ferredoxin_oxidoreductase.fasta")
        multiple_amino_file = str(TEST_FILES / "ferredoxin_oxidoreductase_degen_condons.txt")

        extracted_protein_secquence_data = extract_protein_sequence()


if __name__ == '__main__':
    unittest.main()


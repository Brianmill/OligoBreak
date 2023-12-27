import unittest
import src

class TestDataParser(unittest.TestCase):

    def test_temp(self):
        result = 15
        self.assertEqual(result, 15)


if __name__ == '__main__':
    unittest.main()


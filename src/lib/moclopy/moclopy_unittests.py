import unittest
import moclo as moclopy
from synbiolib import codon

table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)

correct= "atgagtgaagtaaacctaaaaggaaatacagatgaattagtgtattatcgacagcaaaccactggaaataaaatcgccaggaagagaatcaaaaaagggaaagaagaagtttattatgttgctgaaacggaagagaagatatggaca"
contains_sapi_for = "gatatggacagctcttcaaataaaaaa"
contains_sapi_rev = "gatatggacagaagagcaaataaaaaa"
contains_long_repeat = "atgagtgaagtaaacctaaaaggaaatacagatgaattagtgtattatcgacagcaaaccactggaaataaaatcgccaggaagagaatcaaaaaagggaaagaagaagtttattatgttgctgaaacggaagagaagatatggacaatgagtgaagtaaacctaaaaggaaatacagatgaattagtgtattatcgacagcaaaccactggaaataTAA"
contains_short_repeat = "atgGGGCGATCGGCGCGGGGCCCAagtgaagtaaacctaaaaggaaatacagatgaattagtgtattatcgacagcaaaccactggaaataaaatcgccaggaagagaatcaaaaaagggaaagaagaagtttattatgttgctgaaacggaagagaagatatggacaGGGCGATCGGCGCGGGGCCCATAA"
contains_homopolymer = "gatatggacagctcttcaaataaaaaaaaaaGACTAGCTAGCA"
contains_highgc = "atgagtgaagtaaacctaaaaggaaatacagatgaattagtgtattatcgacagcaaaccactggaaataaaatcgccaatgGGGCGATCGGCGCGAagtgaagtaaacctaaaaggaaatacagatgaattagtgtattatcgacagcaaaccactggaaataaGGGGAGGGGGGGGTGCGCGGACGAGCGACGAGCGAGCGCGGCAGGCAGGCaatcgccaggaagagaatcaaaaaagggaaagaagaagtttattatgttgctgaaacggaagagaagatatggacaGGGCGATCGGCGCGATAAtcaaaaaagggaaagaagaagtttattatgttgctgaaacggaagagaagatatggacagaagagcaaataaaaaacttttctttagacaaatttggtacgcatataccttacatagaaggtcattatacaatcttaaataattacttctttgatttttggggctattttttaggtgctgaaggaattgcgctctatgctcacctaactcgttatgcatacggcagcaaagacttttgctttcctagtctacaaacaatcgctaaaaaaatggacaagactcctgttacagttagaggctacttgaaactgcttgaaaggtacggttttatttggaaggtaaacgtccgtaataaaaccaaggataacacagaggaatccccgatttttaagattagacgtaaggttcctttgctttcagaagaacttttaaatggaaaccctaatattgaaattccagatgacgaggaagcacatgtaaagaaggctttaaaaaaggaaaaagagggtcttccaaaggttttgaaaaaagagcacgatgaatttgttaaaaaaatgatggatgagtcagaaacaattaatattccagaggccttacaatatgacacaatgtatgaagatatactcagtaaaggagaaattcgaaaagaaatcaaaaaacaaatacctaatcctacaacatcttttgagagtatatcaatgacaactgaagaggaaaaagtcgacagtactttaaaaagcgaaatgcaaaatcgtgtctctaagccttcttttgatacctggtttaaaaacactaagatcaaaattgaaaataaaaattgtttattacttgtaccgagtgaatttgcatttgaatggattaagaaaagatatttagaaacaattaaaacagtccttgaagaagctggatatgttttcgaaaaaatcgaactaagaaaagtgcaataa"
contains_lowgc = "GAGGAGGGGAGCGGAGGAGCGAGGCAAATCAGGCATGATATGATGATGTAGATGATGATGATGATGATGATGATGATGAGTAGATGGATATGATGATGTAGATGATGATGATGATGATGATGATGATGAGTAGATGaaaaaaaaaaaaaaaaaaaaaaaatttttttttttaaaaaaatttaatatGATATGATGATGTAGATGATGATGATGATGATGATGATGATGAGTAGATGGATATGATGATGTAGATGATGATGATGATGATGATGATGATGAGTAGATGGAGGAGGGGAGCGGAGGAGCGAGGCAAATCAGGCAT"
recode_test_A = "TGGATAGGAGAAAAAAAAATCGGAGATCCAACATAA"

repA = "atgagtgaagtaaacctaaaaggaaatacagatgaattagtgtattatcgacagcaaaccactggaaataaaatcgccaggaagagaatcaaaaaagggaaagaagaagtttattatgttgctgaaacggaagagaagatatggacagaagagcaaataaaaaacttttctttagacaaatttggtacgcatataccttacatagaaggtcattatacaatcttaaataattacttctttgatttttggggctattttttaggtgctgaaggaattgcgctctatgctcacctaactcgttatgcatacggcagcaaagacttttgctttcctagtctacaaacaatcgctaaaaaaatggacaagactcctgttacagttagaggctacttgaaactgcttgaaaggtacggttttatttggaaggtaaacgtccgtaataaaaccaaggataacacagaggaatccccgatttttaagattagacgtaaggttcctttgctttcagaagaacttttaaatggaaaccctaatattgaaattccagatgacgaggaagcacatgtaaagaaggctttaaaaaaggaaaaagagggtctcccaaaggttttgaaaaaagagcacgatgaatttgttaaaaaaatgatggatgagtcagaaacaattaatattccagaggccttacaatatgacacaatgtatgaagatatactcagtaaaggagaaattcgaaaagaaatcaaaaaacaaatacctaatcctacaacatcttttgagagtatatcaatgacaactgaagaggaaaaagtcgacagtactttaaaaagcgaaatgcaaaatcgtgtctctaagccttcttttgatacctggtttaaaaacactaagatcaaaattgaaaataaaaattgtttattacttgtaccgagtgaatttgcatttgaatggattaagaaaagatatttagaaacaattaaaacagtccttgaagaagctggatatgttttcgaaaaaatcgaactaagaaaagtgcaataa"

class EnzymeFinder(unittest.TestCase):
    def test_enzyme_finder_for(self):
        self.assertEqual(moclopy.enzyme_finder(contains_sapi_for), (10, 7))
    def test_enzyme_finder_rev(self):
        self.assertEqual(moclopy.enzyme_finder(contains_sapi_rev), (10, 7))
    def test_No_enzyme(self):
        self.assertFalse(moclopy.enzyme_finder(correct))

class RepeatFinder(unittest.TestCase):
    def test_repeat_finder_long(self):
        self.assertEqual(moclopy.repeat_finder(contains_long_repeat), (0, 70))
    def test_repeat_finder_small(self):
        self.assertEqual(moclopy.repeat_finder(contains_short_repeat), (3, 21))
    def test_No_repeat(self):
        self.assertFalse(moclopy.repeat_finder(correct))

class HomopolymerFinder(unittest.TestCase):
    def test_homopolymer(self):
        self.assertEqual(moclopy.homopolymer_finder(contains_homopolymer), (21, 10))
    def test_No_homopolymer(self):
        self.assertFalse(moclopy.homopolymer_finder(correct))

class GCrangeFinder(unittest.TestCase):
    def test_GCrange_high(self):
        self.assertEqual(moclopy.gc_range_finder(contains_highgc), (165, 50, "AT"))
    def test_GCrange_low(self):
        self.assertEqual(moclopy.gc_range_finder(contains_lowgc), (136, 50, "GC"))
    def test_No_GCrange(self):
        self.assertFalse(moclopy.gc_range_finder(correct))

class DictionaryGen(unittest.TestCase):
    def test_dictionarygen(self):
        self.assertEqual(moclopy.generate_codon_dict(recode_test_A), {0: ('W', ['TGG']),
  1: ('I', ['ATA']),
  2: ('G', ['GGA']),
  3: ('E', ['GAA']),
  4: ('K', ['AAA']),
  5: ('K', ['AAA']),
  6: ('I', ['ATC']),
  7: ('G', ['GGA']),
  8: ('D', ['GAT']),
  9: ('P', ['CCA']),
  10: ('T', ['ACA']),
  11: ('*', ['TAA'])})

class RecodeSequence(unittest.TestCase):
    def test_recode_GC(self): 
        self.assertEqual(moclopy.replace_codons(table,recode_test_A,moclopy.generate_codon_dict(recode_test_A),(0,6),GC_bias="GC"), ('TGGATCGGAGAAAAAAAAATCGGAGATCCAACATAA',
 {0: ('W', ['TGG']),
  1: ('I', ['ATA', 'ATC']),
  2: ('G', ['GGA']),
  3: ('E', ['GAA']),
  4: ('K', ['AAA']),
  5: ('K', ['AAA']),
  6: ('I', ['ATC']),
  7: ('G', ['GGA']),
  8: ('D', ['GAT']),
  9: ('P', ['CCA']),
  10: ('T', ['ACA']),
  11: ('*', ['TAA'])}))
    def test_recode_AT(self):
        self.assertEqual(moclopy.replace_codons(table,recode_test_A,moclopy.generate_codon_dict(recode_test_A),(0,6),GC_bias="AT"), ('TGGATTGGAGAAAAAAAAATCGGAGATCCAACATAA',
 {0: ('W', ['TGG']),
  1: ('I', ['ATA', 'ATT']),
  2: ('G', ['GGA']),
  3: ('E', ['GAA']),
  4: ('K', ['AAA']),
  5: ('K', ['AAA']),
  6: ('I', ['ATC']),
  7: ('G', ['GGA']),
  8: ('D', ['GAT']),
  9: ('P', ['CCA']),
  10: ('T', ['ACA']),
  11: ('*', ['TAA'])}))


if __name__ == '__main__':
    unittest.main()

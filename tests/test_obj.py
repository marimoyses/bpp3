import unittest
import filecmp
import subprocess
from bpp3.bpp3obj import *
from bpp3.helper import *

class MyTest(unittest.TestCase):
	def test_validate(self):
		self.assertEqual(validate("Chr11"), "11")
	def test_constant(self):
		self.assertEqual(primer3Location,'/PHShome/jtg24/PrimerDesign/primer3-2.3.6/src/primer3_core')
	def test_getsequence(self):
		self.assertEqual(getsequence(2,219273660,100),'aaactcaaacctgaccttctgatccacccgccttggcctcccaaagcaccgggataacNggtgtgagccaccgcacctggccAAAAATATTTTTTTAATG')
	def test_goodinput(self):
		file="./testdata/goodinput.txt"
		testobj=BpP3obj(file,9999)
		testobj.sanitate()
		self.assertEqual(testobj.data["PESR_RD_105711"]["START_A"],39273536)
	def test_p3input(self):
		file="./testdata/goodinput.txt"
		testobj=BpP3obj(file,9999)
		testobj.sanitate()
		testobj.getP3input("PESR_RD_757712",75,75,1)
		self.assertEqual(testobj.data["PESR_RD_757712"]['p3input'],'tctaccatcacacccggctaattttttgtatttttagtagagacagggtttcaccatgttggccaggctggtctcggggtttcaccatgttggccaggctggtctcgaactcccaacctcaggtgacctgcctgccttggcctcccaaag')
	def test_badstart(self):
		with self.assertRaises(ValueError):
			file="./testdata/badinput.txt"
			testobj=BpP3obj(file,9999)
			testobj.sanitate()
	def test_printprimer(self):
		file="./testdata/goodinput.txt"
		testobj=BpP3obj(file,9999)
		testobj.sanitate()
		testobj.runP3(1)
		testobj.printprimers("./testdata/test.out")
		subprocess.call("sort ./testdata/test.out > ./testdata/test.sort",shell=True)
		self.assertTrue(filecmp.cmp("./testdata/test.sort","./testdata/test.ref"))
	# def test_qc(self):
		# file="./testdata/goodinput.txt"
		# testobj=BpP3obj(file,9999)
        # testobj.startBLATserver()
		# testobj.sanitate()
		# testobj.runP3(1)
		# testobj.QC()
        # testobj.printbest3primer("")
		# subprocess.call("sort ./testdata/test.out > ./testdata/test.sort",shell=True)
		# self.assertTrue(filecmp.cmp("./testdata/test.sort","./testdata/test.ref"))

suite = unittest.TestLoader().loadTestsFromModule(MyTest())
unittest.TextTestRunner().run(suite)

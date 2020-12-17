all: ttbar_pos ttbar_neg mumu qqbar bbbar gmZ WW ZZ ZWW ZZZ singleTop
.PHONY: all

ttbar_pos:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_ttbar_pos_ecm365.root

ttbar_neg:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_ttbar_neg_ecm365.root

mumu:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_mumu_ecm365.root

qqbar:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_qqbar_ecm365.root

bbbar:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_bbbar_ecm365.root

gmZ:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_gmZ_ecm365.root

WW:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_WW_ecm365.root

ZZ:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_ZZ_ecm365.root

ZWW:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_ZWW_ecm365.root

ZZZ:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_ZZZ_ecm365.root

singleTop:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_singleTop_ecm365.root


all: ttbar_semi mumu qqbar bbbar gmZ WW ZZ ZWW ZZZ singleTop
.PHONY: all

all1M: ttbar_pos1M ttbar_neg1M mumu1M qqbar1M bbbar1M gmZ1M WW1M ZZ1M ZWW1M ZZZ1M singleTop1M
.PHONY: all1M

ttbar_semi:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py ../ttbar/p8_ee_ttbar_semi_ecm365.root

mumu:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py ../ttbar/p8_ee_mumu_ecm365.root

qqbar:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py ../ttbar/p8_ee_qqbar_ecm365.root

bbbar:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py ../ttbar/p8_ee_bbbar_ecm365.root

gmZ:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py ../ttbar/p8_ee_gmZ_ecm365.root

WW:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py ../ttbar/p8_ee_WW_ecm365.root

ZZ:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py ../ttbar/p8_ee_ZZ_ecm365.root

ZWW:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py ../ttbar/p8_ee_ZWW_ecm365.root

ZZZ:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py ../ttbar/p8_ee_ZZZ_ecm365.root

singleTop:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py ../ttbar/p8_ee_singleTop_ecm365.root



ttbar_pos1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_ttbar_pos_ecm365.root

ttbar_neg1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_ttbar_neg_ecm365.root

mumu1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_mumu_ecm365.root

qqbar1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_qqbar_ecm365.root

bbbar1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_bbbar_ecm365.root

gmZ1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_gmZ_ecm365.root

WW1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_WW_ecm365.root

ZZ1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_ZZ_ecm365.root

ZWW1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_ZWW_ecm365.root

ZZZ1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_ZZZ_ecm365.root

singleTop1M:
	python FCCeeAnalyses/ttbar/ttbarstreamer.py /eos/user/j/jutornda/ttbar/p8_ee_singleTop_ecm365.root


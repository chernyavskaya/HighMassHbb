import math
import ROOT


files = []
files.append(ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/HighMassHbb/plotterOutput/v21/BTagCSV_analysis_v21_TriggerFinal.root"))
files.append(ROOT.TFile.Open("/afs/cern.ch/work/n/nchernya/HighMassHbb/HighMass/plotter/Top_analysis_v21_TriggerFinal.root"))
num_masses = 11
masses = [ '375','450','525','600','675','750','825','900','975','1050']
for i in range(0,len(masses)):
	name = 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/HighMassHbb/plotterOutput/v21/VBFSpin0ToBBbar_W_1p0_M_' + masses[i] + '_analysis_v21_TriggerFinal.root'
#	files.append(ROOT.TFile.Open(name))

out_hists =[]

for n_files,o in enumerate(files):
	print n_files
	hname = "hMbb_met_fsr"
	if (n_files==0) :	hist =  o.Get(hname).Clone("data_"+hname)
	if (n_files==1) :	hist =  o.Get(hname).Clone("top_"+hname)
#	else :  hist =  o.Get(hname).Clone("signal_"+hname+"_M_"+masses[n_files-1])
	out_hists.append(hist)
	hname = "hMbb_met_fsr_bg"
	if (n_files==0) :	hist =  o.Get(hname).Clone("data_"+hname)
	if (n_files==1) :	hist =  o.Get(hname).Clone("top_"+hname)
#	else :  hist =  o.Get(hname).Clone("signal_"+hname+"_M_"+masses[n_files-1])
	out_hists.append(hist)
	hname = "hMbb_met_fsr_sg"
	if (n_files==1) :	hist =  o.Get(hname).Clone("top_"+hname)
	if (n_files==0) :	hist =  o.Get(hname).Clone("data_"+hname)
#	else :  hist =  o.Get(hname).Clone("signal_"+hname+"_M_"+masses[n_files-1])
	out_hists.append(hist)

fout = ROOT.TFile("inputs.root","recreate")
for hist in out_hists:
	hist.Write()


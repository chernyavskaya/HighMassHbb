export WORKDIR=`pwd`
cd $WORKDIR

g++ plotter_highmass.C -g -o plot `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA  

declare -A file_names

cp plot /mnt/t3nfs01/data01/shome/nchernya/HighMass/plotter/
cp  batch.sh /mnt/t3nfs01/data01/shome/nchernya/HighMass/plotter/	

path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/HighMassHbb/main_mva/
#path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store/user/nchernya/Hbb/v21/work/
#path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/HighMassHbb/workflow/
file_names=(
#["VBFSpin0ToBBbar_W_1p0_M_750"]=VBFSpin0ToBBbar_W_1p0_M_750_TuneCUEP8M1_13TeV_pythia8
##["VBFSpin0ToBBbar_W_1p0_M_300"]=VBFSpin0ToBBbar_W_1p0_M_300_TuneCUEP8M1_13TeV_pythia8
##["VBFSpin0ToBBbar_W_1p0_M_375"]=VBFSpin0ToBBbar_W_1p0_M_375_TuneCUEP8M1_13TeV_pythia8
#["VBFSpin0ToBBbar_W_1p0_M_450"]=VBFSpin0ToBBbar_W_1p0_M_450_TuneCUEP8M1_13TeV_pythia8
#["VBFSpin0ToBBbar_W_1p0_M_525"]=VBFSpin0ToBBbar_W_1p0_M_525_TuneCUEP8M1_13TeV_pythia8
#["VBFSpin0ToBBbar_W_1p0_M_600"]=VBFSpin0ToBBbar_W_1p0_M_600_TuneCUEP8M1_13TeV_pythia8
#["VBFSpin0ToBBbar_W_1p0_M_675"]=VBFSpin0ToBBbar_W_1p0_M_675_TuneCUEP8M1_13TeV_pythia8
#["VBFSpin0ToBBbar_W_1p0_M_825"]=VBFSpin0ToBBbar_W_1p0_M_825_TuneCUEP8M1_13TeV_pythia8
#["VBFSpin0ToBBbar_W_1p0_M_900"]=VBFSpin0ToBBbar_W_1p0_M_900_TuneCUEP8M1_13TeV_pythia8
#["VBFSpin0ToBBbar_W_1p0_M_975"]=VBFSpin0ToBBbar_W_1p0_M_975_TuneCUEP8M1_13TeV_pythia8
#["VBFSpin0ToBBbar_W_1p0_M_1050"]=VBFSpin0ToBBbar_W_1p0_M_1050_TuneCUEP8M1_13TeV_pythia8
["BTagCSV"]=BTagCSV
#["QCD_HT100to200"]=QCD_HT100to200
#["QCD_HT200to300"]=QCD_HT200to300
#["QCD_HT300to500"]=QCD_HT300to500
#["QCD_HT500to700"]=QCD_HT500to700
#["QCD_HT700to1000"]=QCD_HT700to1000
#["QCD_HT1000to1500"]=QCD_HT1000to1500
#["QCD_HT1500to2000"]=QCD_HT1500to2000
#["QCD_HT2000toInf"]=QCD_HT2000toInf
#["ST_tW"]=ST_tW_topantitop_5f_inclusiveDecays_13TeV-powheg
#["ST_t-channel_top_4f_inclusiveDecays"]=ST_t-channel_top_4f_inclusiveDecays
#["ST_t-channel_antitop_4f_inclusiveDecays"]=ST_t-channel_antitop_4f_inclusiveDecays
#["ST_s-channel"]=ST_s-channel_4f_leptonDecays_13TeV-amcatnlo
#["DYJetsToLL"]=DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnlo_all
#["WJetsToLNu"]=WJetsToLNu_TuneCUETP8M1_13TeV-amcatnlo
#["WJetsToQQ"]=WJetsToQQ_HT180_13TeV-madgraph
#["DYJetsToQQ"]=DYJetsToQQ
####["TT"]=TT_TuneCUETP8M1_13TeV-amcatnlo
######["TT_madgraph"]=TTJets_TuneCUETP8M1_13TeV-madgraph
#["TT_powheg"]=TT_powheg
#["WW"]=WW
#["ZZ"]=ZZ
#["WZ"]=WZ
)
prefix='main_mva_v21_'
postfix=TriggerFinal
v=v21
ROOT=.root
region=(analysis controlTop controlDY)
applyTrigWeight=1
trigWeightNom=(nom up down)
output_dir=$TMPDIR


for key in ${!file_names[@]}; do
	current_region=0
	while [ $current_region -lt 1  ] 
#	while [ $current_region -lt  3 ] 
	do
		data=0
		if [ $key == BTagCSV ]
		then 
	 		data=1
		fi
		if [ $key == BTagCSV ] || [ $applyTrigWeight -eq 0 ]
		then
			f=$path$prefix${file_names[${key}]}_$v.root
			qsub -q short.q batch.sh  $f ${key} ${region[$current_region]} $data  0 nom $v $postfix
		fi 
		if [ $key != BTagCSV ] && [ $applyTrigWeight -eq 1 ]
		then
			current_trigWeight=0
			while [ $current_trigWeight -lt  1 ] 
	#		while [ $current_trigWeight -lt  3 ] 
			do
				f=$path$prefix${file_names[${key}]}_$v.root
				qsub -q short.q batch.sh  $f ${key} ${region[$current_region]} $data  $applyTrigWeight nom $v $postfix
				current_trigWeight=$(( $current_trigWeight + 1 ))
			done
		fi
		current_region=$(( $current_region + 1 ))
#		break
	done
#	break
done

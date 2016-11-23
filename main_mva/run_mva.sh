export WORKDIR=`pwd`
cd $WORKDIR

g++ run_create_main_tmva_all.C -g -o run_all `root-config --cflags --glibs` 

max_samples_num=12 
path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/HighMassHbb/workflow/
#input_dir=( BTagCSV VBFHToBB_M-125_13TeV_powheg-ext)
#input_dir=( BTagCSV_single_v21 VBFSpin0ToBBbar_W_1p0_M_300_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_375_TuneCUEP8M1_13TeV_pythia8_v21  VBFSpin0ToBBbar_W_1p0_M_450_TuneCUEP8M1_13TeV_pythia8_v21  VBFSpin0ToBBbar_W_1p0_M_525_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_600_TuneCUEP8M1_13TeV_pythia8_v21  VBFSpin0ToBBbar_W_1p0_M_675_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_750_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_825_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_900_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_975_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_1050_TuneCUEP8M1_13TeV_pythia8_v21   )
input_dir=( BTagCSV_v21  VBFSpin0ToBBbar_W_1p0_M_375_TuneCUEP8M1_13TeV_pythia8_v21  VBFSpin0ToBBbar_W_1p0_M_450_TuneCUEP8M1_13TeV_pythia8_v21  VBFSpin0ToBBbar_W_1p0_M_525_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_600_TuneCUEP8M1_13TeV_pythia8_v21  VBFSpin0ToBBbar_W_1p0_M_675_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_750_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_825_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_900_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_975_TuneCUEP8M1_13TeV_pythia8_v21 VBFSpin0ToBBbar_W_1p0_M_1050_TuneCUEP8M1_13TeV_pythia8_v21   )
#data=(1 0 0 0 0 0 0 0 0 0 0 0 )
data=(1  0 0 0 0 0 0 0 0 0 0 )
ROOT=.root


current_sample=1
max_samples_num=11
while [ $current_sample -lt $max_samples_num ]
do	
	./run_all $path${input_dir[ $current_sample ]}$ROOT ${input_dir[ $current_sample ]} ${data[ $current_sample ]}
	current_sample=$(( $current_sample + 1 ))
done

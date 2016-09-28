universe = vanilla
Executable=makeplotellipes.sh
Should_Transfer_Files = YES
Transfer_Input_Files = /eos/uscms/store/user/tahuang/SLHC23_patch1_2023Muon_gen_sim_Pt2_50_1M/GEMCSCAna_ctau0_Pt2_50_0925_v2/160926_194628/0000/combined/out_ana_prompt.root,/eos/uscms/store/user/tahuang/DarkSUSY_MH-125_MGammaD-20000_ctau1000_14TeV_madgraph-pythia6-tauola/GEMCSCAna_DarkSUSY_ctau1000_0925_v2/160926_194711/0000/combined/out_ana_displaced.root 
WhenToTransferOutput  = ON_EXIT_OR_EVICT
Output = condorout/condor_out_$(Process)_$(Cluster).stdout
Error  = condorout/condor_out_$(Process)_$(Cluster).stderr
Log    = condorout/condor_out_$(Process)_$(Cluster).log
Notification    = Error
Arguments = Profile-Pt15-BGPt5-0922 Hybrid-Pt15-BGPt5-0922
use_x509userproxy = true
Queue 1

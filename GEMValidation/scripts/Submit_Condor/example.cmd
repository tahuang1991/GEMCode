universe = vanilla
Executable=makeplotellipes.sh
Should_Transfer_Files = YES
WhenToTransferOutput  = ON_EXIT_OR_EVICT
Transfer_Input_Files = /eos/uscms/store/user/tahuang/SLHC23_patch1_2023Muon_gen_sim_Pt2_50_1M/GEMCSCAna_ctau0_Pt2_50_0925_v2/160926_194628/0000/combined/input.root
Output = condor_out_$(Process)_$(Cluster).stdout
Error  = condor_out_$(Process)_$(Cluster).stderr
Log    = condor_out_$(Process)_$(Cluster).log
Notification    = Error
Arguments = FNAME OUTPUTFILE
Queue 1

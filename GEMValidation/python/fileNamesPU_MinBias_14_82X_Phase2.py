import FWCore.ParameterSet.Config as cms
import os

def getFiles(my_dir):
    theInputFiles = []
    if not os.path.isdir(my_dir):
        print "ERROR: This is not a valid directory: ", my_dir
    for file in os.listdir(my_dir):
        if file.endswith('root'):
            if os.path.getsize(my_dir + file) > 20000000:
                #theInputFiles.extend(['file:' + my_dir + file])
	    	theInputFiles.extend(['/store/group/lpcgem/'+my_dir[29:]+file])
    return theInputFiles

fileNamesPU = cms.untracked.vstring()

fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v1/"))
#print "fileNamePU ",fileNamesPU
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v2/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v3/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v4/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v5/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v6/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v7/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v8/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v9/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v10/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170228/GENSIM_10M_v11/"))

fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170404_10M_v1/GENSIM_10M_v1/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170404_10M_v1/GENSIM_10M_v10/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170404_10M_v1/GENSIM_10M_v2/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170404_10M_v1/GENSIM_10M_v3/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170404_10M_v1/GENSIM_10M_v4/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170404_10M_v1/GENSIM_10M_v5/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170404_10M_v1/GENSIM_10M_v6/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170404_10M_v1/GENSIM_10M_v7/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170404_10M_v1/GENSIM_10M_v8/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170404_10M_v1/GENSIM_10M_v9/"))

fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170407_10M_v1/GENSIM_10M_v1/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170407_10M_v1/GENSIM_10M_v10/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170407_10M_v1/GENSIM_10M_v2/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170407_10M_v1/GENSIM_10M_v3/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170407_10M_v1/GENSIM_10M_v4/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170407_10M_v1/GENSIM_10M_v5/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170407_10M_v1/GENSIM_10M_v6/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170407_10M_v1/GENSIM_10M_v7/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170407_10M_v1/GENSIM_10M_v8/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170407_10M_v1/GENSIM_10M_v9/"))

fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v0/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v1/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v10/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v11/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v12/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v13/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v14/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v15/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v16/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v2/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v3/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v4/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v5/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v6/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v7/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v8/"))
fileNamesPU.extend(getFiles("/eos/uscms/store/user/lpcgem/ME0TDRStudies/MinBias_20170406_10M_v1_Tao/GENSIM_10M_v9/"))

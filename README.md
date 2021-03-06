##Instructions of calculating csv Scale Fators

###Step 0: get a GitHub account. 
If you already have one, go to step 1.

###Step 1: Check out the ttH BEAN code
following instructions here: https://github.com/cms-ttH/BEAN

###Step 2: Check out treeMaking code
Following instructions here: https://github.com/cms-ttH/ttHMultileptonAnalysis

###Step 3: Making trees

 1. tree making script:  ttHMultileptonAnalysis/TemplateMakers/bin/csvSF.C
 2. file lists for the skims of all samples: ttHMultileptonAnalysis/listsForSkims2012_53x_v2_hadoop/*list

// you need to change the path of the skim files accordingly in *list

#### $ To run a test job

> cd ttHMultileptonAnalysis/TemplateMakers/test/

> csvSF ssCondor.py wz test NA 1 10

// csvSF ssCondor.py sample lable JES_syst(NA, JESUp, or JESDown) jobNumber totalJobs

// This runs a test job on sample zz, year 2012_53x, job_label test, JES_syst nominal, job number 1 out of 10 total jobs. The command runs locally, the output goes in batchBEAN/zz_test/  


#### $ To run the full set of jobs, do the following (runs on batch queue!): 

> ./submit_condor_jobs.py csvSF label -j JES_syst

// this submits one group of jobs for each list file in csvSF_lists.txt. The number of jobs in the  group is one for each line in the list file.

// The output goes in batchBEAN/${sample}_${label}/ 

#### $ hadd output root files 
 * datasets

> hadd -f TwoMuon.root DoubleMu_Run2012**{*label}/*root

> hadd -f TwoEle.root DoubleElectron_Run2012**{*label}/*root

> hadd -f MuonEle.root MuEG_Run2012**{*label}/*root

 * MC

> hadd -f allMC.root  ttbar*{label}/*root singlet*{label}/*root wjets*{label}/*root zjets*{label}/*root ww*{label}/*root wz*{label}/*root zz*{label}/*root 

###Step 4: Calculate csvSFs

1. #### Check out csv SF calculation code: 

 > cd $CMSSW_BASE/src

 > git clone https://github.com/cms-ttH/csvReweighting.git

 > cd csvReweighting/

 > mkdir dataMCRootFiles

 (mv the root files from hadd above to dataMCRootFiles)

2. #### version0 SF calculation
 *  light flavor scale factor

 > root -l lightflavorCSVSF.C

 *  heavy flavor scale factor

 > root -l heavyflavorCSVSF.C

  the scale factors version 0 will be stored at **csv_rwt_hf_v0.root** and **csv_rwt_lf_v0.root**

3. #### Do multiple Iterations  (tricky part)

* 3a: replace SF files (x --> iteration version number)

    > mv csv_rwt_lf_v{x}.root BEAN/BEANmaker/data/csv_rwt_lf_IT.root

    > mv csv_rwt_hf_v{x}.root BEAN/BEANmaker/data/csv_rwt_hf_IT.root

* 3b: remake the MC trees following Step 3, new csv weigths (**csvWgtlf, csvWgthf**) based on SFs version{x} will be saved

 //don't need to rerun data trees since there will be no change

* 3c: recalculate SFs version{x+1} following Steps 4 using SFs version{x}, 

 //-- change the input in lightflavorCSVSF.C/heavyflavorCSVSF.C accordingly

 * input 1 --> the version number

 * input 2 --> the input MC root file 

 // the scale factors version x+1 will be stored at **csv_rwt_hf_v{x+1}.root** and **csv_rwt_lf_v{x+1}**.root


###Step5: SF fitting and extrapolation

1. #### Create HistoFiles directory and put all inputs there 

 > mkdir HistoFiles

 > cp $CMSSW_BASE/src/BEAN/BEANmaker/data/mc_btag_efficiency_8TeV_53x.root HistoFiles/

 > cp csv_rwt_*f_IT_v2*.root HistoFiles/

2. #### Run scripts to create Statistics scale factors and wait for that to finish

 > root -b -q fit_csvSF_HF_statErr.C+ > & ! out_fit_csvSF_HF_statErr.log &

 > root -b -q fit_csvSF_LF_statErr.C+ > & ! out_fit_csvSF_LF_statErr.log &

3. #### Run scripts to do the fitting and making image files

 > mkdir Images

 > root -b -q fit_csvSF_HF.C+

 > root -b -q fit_csvSF_LF.C+

4. #### Remove unnecessary histograms from files to reduce size, bloat

 > root -b -q drop_histograms.C+'("HistoFiles/csv_rwt_hf_final_IT_v2.root","csv_rwt_hf_final_IT_v2.root")'

 > root -b -q drop_histograms.C+'("HistoFiles/csv_rwt_lf_final_IT_v2.root","csv_rwt_lf_final_IT_v2.root")'



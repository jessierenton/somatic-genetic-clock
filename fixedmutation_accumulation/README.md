# Fixed mutation accumulation
Data scripts, data and plotting scripts for Figure 3 and Supplementary Figures 2-4, 6, 11-12.

*Some data files used in supplementary figures are excluded from the repository due to size contraints but are available on request.

## Figure 3 and Supplementary Figures 11-12

To produce figures: `julia plot_mainfig.jl`

### Data and scripts:

Fig. 3a: 
  - directory: `data/fulldist_data/`
  - data: `data_branching/`
  - script: `qsub run_moran_branching.sh` 
  - submits `julia --project=$ENV -t 1 run_moran_fulldist.jl $INPUT_ARGS 120 1 20 withoutreplacement data_branching` 
  -  input args: `input.txt`

Fig. 3b and Supplementary Fig. 11a & 12:
  - directory: `data/averaged_data/`
  - data: `data_branch_eelgrass_new/`
  - script: `run_moran_branch_eelgrass_new.sh`
  - submits: `julia --project=$ENV -t 1 run_moran2.jl $INPUT_ARGS 30 0.5 20 data_branch_eelgrass_new`
  - input args: `input_branch_eelgrass.txt`
  
Fig. 3c and Supplementary Fig. 11b: 
  - directory: `data/experimental_sampling_data/`
  - data: `data/`
  - script: `qsub run.sh`
  - submits `julia --project=$ENV -t 1 run_experiment_samples.jl $INPUT_ARGS 30 $JOB_ID data` 
  - input args: `input_experimental_samples.txt`

Fig. 3c and Supplementary Fig. 11b (stars)
  - directory: `data/averaged_data/`
  - data: `data_branch_eelgrass_long/X/` (X=0:4)
  - script: `qsub run_moran_branch_eelgrass_long.sh X`
  - submits: `julia --project=$ENV -t 1 run_moran2.jl $INPUT_ARGS 200 10 20 data_branch_eelgrass_long/$X`
  - input args: `input_branch_eelgrass.txt`

## Supplementary Figure 2-3

To produce figures: `julia plot_supp_compare_r.jl`

### Data and scripts:

- directory: `data/fulldist_data/`
- data: `data_long/` (*not included in repository due to size)
- script: `run_moran_long.sh`
- submits: `julia --project=$ENV -t 1 run_moran_fulldist.jl $INPUT_ARGS 480 10 100 split data_long`
- input args: `input.txt`

## Supplementary Figure 4

To produce figure: `julia plot_supp_grads.jl`

### Data and scripts:

Module splitting:
- directory: `data/fulldist_data/` 
- data: `data_long/` (*not included in repository due to size)
- see Supp Fig 2-3

Module branching:
- directory: `data/fulldist_data/`
- data: `data_long_branching/` (*not included in repository due to size)
- script: `run_moran_long_branching.sh`
- submits: `julia --project=$ENV -t 1 run_moran_fulldist.jl $INPUT_ARGS 480 10 20 withoutreplacement_nomutations data_long_branching`
- input args: `input.txt`

Fitted linear regression in:
`data/fulldist_data/linear_fit_fulldist_long_split.csv`
`data/fulldist_data/linear_fit_fulldist_long_branch.csv`


## Supplementary Figure 6

To produce figure: `julia plot_supp_trajectories.jl`

### Data and scripts

- directory `data/fulldist_data/`
- data: `data_tenyrs/` (*not included in repository due to size)
- script: `run_moran_tenyrs.sh`
- submits: `julia --project=$ENV -t 1 run_moran_fulldist.jl $INPUT_ARGS 10 0.1 20 split data_tenyrs`
- input args: `input.txt`

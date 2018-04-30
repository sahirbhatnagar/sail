# README for Simulation Files

This file describes where the simulation results are located, and how some of the figures were
generated.

The code in
"/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail"
is upto date with master on GitHub


1. /mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/mcgill_simulations_genedata.R
has the simulation study for sail only with paramIndex=1. This was
run as a job array on HYDRA using the bash script mcgill-simulations.sh. These results were saved to
/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail_lambda_branch/mcgillsims/thesis_p1000_1a/fit*

These results were used to generate the plots that show the estimated vs. true smooth functions (both main effects and interactions).
The code used to produce those plots is in
/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/mcgill_simulations_analysis.R


2. The full simulation study was done in the file:
/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/main.R
This was run in parallel on an rstudio node, and in tmux.

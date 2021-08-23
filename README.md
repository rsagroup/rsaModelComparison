# RSA model comparison

This repository contains the Matlab code and data to conduct the simulations to compare different method for the comparison of representational models. The simulations are presented in 2 papers: 

* Diedrichsen, J., & Kriegeskorte, N. (2017). Representational models: A common framework for understanding encoding, pattern-component, and representational-similarity analysis. *PLoS Comput Biol.* 
* Diedrichsen, J., Berlot, E., Mur, M., Schütt , H. S., Shahbazi, M., Kriegeskorte, N. (2021). Comparing representational geometries using whitend unbiased distance-matrix similarity. *NBDT*.

## Installation and requirements 

Running the code requires two additional Matlab toolboxes from our group: 

* Dataframe toolbox [https://github.com/jdiedrichsen/dataframe/releases]
* rsatoolbox-matlab [https://github.com/rsagroup/rsatoolbox_matlab] - **develop branch**
* PCM toolbox [https://github.com/jdiedrichsen/pcm_toolbox]

## Conducting simulations 

```rsa_testModelCompare('modelCompare',....);```

is the main method to run and sace simulation. The function simulates multivariate data from 2 or more specific models (specified by the `model` file) and then fits them with these models, using different model comparision `methods`. The signal-to-noise ratio `Omega` can be varied for the simultion and the results are written to `outfile`. For example to conduct the stimulations for experiment 1, using different RSA and PCM models, you can call: 
```
RSA_methods={'spearman','pearson','fixed','loglikPCM'};
rsa_testModelCompare('modelCompare','model','Model_fiveFinger.mat','numSim',1000,'outfile','sim_rsa_Exp1.mat','methods',RSA_methods,'Omega',[0:0.1:0.8]);
```


## Reproducing Figures from the papers

To reproduce the Figures, you first have to conduct the underlying simulations (see commented lines in case 'modelCompare'. Once these files are available in the root folder, the Figures can be reproduced using: 

For the Diedrichsen & Kriegeskorte (2017) paper, the relevant cases for reproducing the Figures are: 
* Fig. 4: `rsa_testModelCompare('Figure_models');`
* Fig. 5: `rsa_testModelCompare('Figure_numReg');`
* Fig. 6: `rsa_testModelCompare('Figure_regularisation');`
* Fig. 7: `rsa_testModelCompare('Figure_pcm_rsa_encode');`
* Fig. 8: `rsa_testModelCompare('Figure_bestMethods');`

For Diedrichsen et al. (2021), the relevant functions are for: 
* Fig 1: `rsa_covarianceDistPaper('Figure1');`
* Fig 3: `rsa_covarianceDistPaper('Figure_variancebias');`
* Fig 4: `rsa_covarianceDistPaper('Figure_crossval_noncrossval');`
* Fig 5: `rsa_covarianceDistPaper('Figure_covariances');`
* Fig 6: `rsa_covarianceDistPaper('Figure_rsa_weight');`
* Fig 7: `rsa_covarianceDistPaper('Figure_numberOfK');`
* Fig R6: `rsa_covarianceDistPaper('Figure_R6_SpatialCorr');`
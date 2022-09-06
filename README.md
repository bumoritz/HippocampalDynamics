# HippocampalDynamics

This repository contains analysis code for my main PhD project, which is about hippocampal dynamics for associating events across time. It's pretty much work-in-progress. A streamlined and annotated version of the code as well as the data will be made available upon publication of the paper.


# Project overview

Episodic memory is the ability to remember events in their spatial and temporal context. Forming new episodic memories requires associating events across time and is known to depend on the hippocampus. We investigate whether the hippocampus supports this ability by generating neural activity sequences which tile the temporal gap between events to be linked. 


### Part 1 - Hippocampal dynamics during learning, retrieval and generalisation
We record the activity of hundreds of neurons in the hippocampus while mice learn and retrieve associations between temporally separated events. To study how the hippocampus generalizes task rules for learning new stimulus associations, we then let expert mice perform the same conceptual task with different stimuli they were never exposed to before.

<p align="center">
  <img src="/Img/ProjectOverview_Part1.png" width="700">
</p>


### Part 2 - Binding memories onto artificial memory traces
An emerging view is that binding event representations onto hippocampal sequence templates lays the basis for acquiring new episodic memories. We directly test this hypothesis by stimulating an arbitrary hippocampal sequence while mice have new experiences. We analyze whether this imprinted sequence might acquire meaning to the brain by starting to encode the new experience.

<p align="center">
  <img src="/Img/ProjectOverview_Part2.png" width="700">
</p>


# Code structure

For each master script, specify which modules to run in the `ops` structure and set parameters in the `/utils/defaults/get_p.m` file.

* `Data2repo_Master.m`: This master script controls the pre-processing of data from an individual experimental session (e.g. suite2p output, synchronisation data, behaviour data, camera data) and saves the relevant information to a repository.
* `Repo2repo_Master.m`: This master script controls further pre-processing of the data in the repository (e.g. tracking the same neurons across days, re-identifying photostimulation targets).
* `Analyses_Master.m`: This master script controls all analyses for individual experimental sessions (e.g. neuronal encoding analysis, decoding analysis, population vector analysis, low-dimensional analysis, sequence analysis, inhibition analysis, response analysis, behavioural analysis, learning analysis)
* `Summary_Master.m`: This master script loads in the analysis outputs from all experimental sessions and controls dataset-wide analyses. 


# Example analyses


### Which behavioural features does each individual neuron encode? -> GLM

`/analyses/nemAnalysis/nemAnalysis.m`




### What is the neural population as a whole encoding? -> Bayesian decoder

`/analyses/decodingAnalysis/decodingAnalysis.m`

<p align="center">
  <img src="/Img/ExampleGLM.png" width="1400">
</p>

This figure summarises the model setup for a GLM fitting task events and behavioural predictors to the neural activity traces of each individual neuron. The depicted predictor image is an excerpt (6 trials) from the design matrix, containing discrete and differently delayed perdictors (odours A, X, B, Y; reward; lick) and continuous behavioural predictors (sniffing, pupil, velocity, acceleration, video motion energy). The model is fit on 80% of the data with 5-fold cross-validation using elastic net regularisation and tested on the 20% of hold-out trials.


### Do the neurons targeted for photostimulation indeed become activated? -> Response analysis

`/analyses/responseAnalysis/responseAnalysis.m`

<p align="center">
  <img src="/Img/ExampleResponseAnalysis.png" width="1400">
</p>




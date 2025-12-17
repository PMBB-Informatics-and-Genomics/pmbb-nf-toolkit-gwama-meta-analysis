
Documentation for GWAMA Meta-Analysis
=====================================

# Module Overview


GWAMA meta-analysis uses the advanced statistical software, GWAMA, for performing population-adjusted meta-analysis of GWAS summary statistics
- [Tool Paper Link for Reference](https://link.springer.com/article/10.1186/1471-2105-11-288)
- [Tool Documentation Link for Reference](https://genomics.ut.ee/en/tools)
- [Example Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/gwama_meta.config)
- [Example nextflow.config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config)

## Software Requirements


* [Nextflow version 24.04.3](https://www.nextflow.io/docs/latest/cli.html)

<<<<<<< HEAD
=======
[Example nextflow.config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config)
## Cloning Github Repository


* Command: `git clone https://github.com/PMBB-Informatics-and-Genomics/geno_pheno_workbench.git`

* Navigate to relevant workflow directory...
## Software Requirements


* [Nextflow version 23.04.1.5866](https://www.nextflow.io/docs/latest/cli.html)

>>>>>>> a5ba8f90b2a7163dc8a0d6179b32b0c1497ed225
* [Singularity 3.8.3](https://sylabs.io/docs/) OR [Docker 4.30.0](https://docs.docker.com/)
## Commands for Running the Workflow


* Singularity Command: `singularity build gwama_meta.sif docker://pennbiobank/gwama_meta:latest`

* Docker Command: `docker pull pennbiobank/gwama_meta:latest`

* Pull from Google Container Registry: `docker pull pennbiobank/gwama_meta:latest`

* Run Command: `nextflow run /path/to/toolkit/module/gwama_meta.nf`

* Common `nextflow run` flags:

    * `-resume` flag picks up workflow where it left off

    * `-stub` performs a dry run, checks channels without executing code

    * `-profile` selects the compute profiles in nextflow.config

    * `-profile standard` uses the Docker image to execute processes

    * `-profile cluster` uses the Singularity container and submits processes to a queue

    * `-profile all_of_us` uses the Docker image on All of Us Workbench

* More info: [Nextflow documentation](https://www.nextflow.io/docs/latest/cli.html)
<<<<<<< HEAD
# Detailed Pipeline Steps

## Part I: Setup


1. Start your own tools directory and go there. You may do this in your project analysis directory, but it often makes sense to clone into a general `tools` location

```sh
# Make a directory to clone the pipeline into
TOOLS_DIR="/path/to/tools/directory"
mkdir $TOOLS_DIR
cd $TOOLS_DIR
```

2. Download the source code by cloning from git

```sh
git clone None
cd $TOOLS_DIR/pmbb-nf-toolkit-gwama-meta
```

3. Build the singularity image
    - you may call the image whatever you like, and store it wherever you like. Just make sure you specify the name in `nextflow.conf`
    - this does NOT have to be done for every saige-based analysis, but it is good practice to re-build every so often as we update regularly.


```sh
cd $TOOLS_DIR/pmbb-nf-toolkit-gwama-meta
singularity build gwama_meta.sif docker://pennbiobank/gwama_meta:latest
```
## Part II: Configure your run


1. Make a separate analysis/run/working directory.
    - The quickest way to get started, is to run the analysis in the folder the pipeline is run. However, subsequent analyses will over-write results from previous analyses.
    - ❗This step is optional, but We Highly recommend making a `tools` directory separate from your `run` directory. We recommend storing the `nextflow.conf` in here as it shouldn't change between runs.


```sh
WDIR="/path/to/analysis/run1"
mkdir -p $WDIR
cd $WDIR
```

2. Fill out the `nextflow.config` file for your system.
    - See [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html) for information on how to configure this file. An example can be found on our GitHub: [Nextflow Config](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/blob/main/Example_Configs/nextflow.config).
    - ❗IMPORTANTLY, you must configure a user-defined profile for your run environments (local, docker, saige, cluster, etc.). If multiple profiles are specified, run with a specific profile using `nextflow run -profile $MY_PROFILE`.
    - For singularity, The profile's attribute `process.container` should be set to `'/path/to/gwama_meta.sif'` (replace `/path/to` with the location where you built the image above). See [Nextflow Executor Information](https://www.nextflow.io/docs/latest/executor.html) for more details.
    - ⚠️As this file remains mostly unchanged for your system, We recommend storing this file in the `tools/pipeline` directory and passing it to the pipeline with `-c /path/to/nextflow.config`.


3. Create a pipeline-specific `.config` file specifying your run parameters and input files. See Below for workflow-specific parameters and what they mean.
    - Everything in here can be configured in `nextflow.config`, however we find it easier to separate the system-level profiles from the individual run parameters.
    - Examples can be found in our Pipeline-Specific [Example Config Files](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs).
    - you can compartamentalize your config file as much as you like by passing
    - There are 2 ways to specify the config file during a run:

        - with the `-c` option on the command line: `nextflow run -c GWAMA_META/gwama_meta.config`
        - in the `nextflow.config`: at the top of the file add: `includeConfig GWAMA_META/gwama_meta.config`

## Part III: Run your analysis
=======
# Input Files for GWAMA_Meta-Analysis


* GWAMA Executable

    * Path to the GWAMA executable

    * Type: Executable

    * Format: exe

* GWAS Summary Statistics

    * Files with GWAS summary stats to be meta-analyzed. 

    * Type: Summary Statistics

    * Format: txt.gz

    * File Header:


    ```
    chromosome      base_pair_location      variant_id      other_allele    effect_allele   effect_allele_count     effect_allele_frequency missing_rate    beta    standard_error      t_statistic     variance        p_value N
    1       722858  chr1:722858:C:T C       T       59.0157 0.0010625       0       0.364206        0.178679        11.4078 31.3223 0.04151698      27772
    1       764898  chr1:764898:G:A G       A       54.5843 0.000982722     0       0.395419        0.189646        10.9943 27.8043 0.03706621      27772
    1       766367  chr1:766367:G:A G       A       57.851  0.00104153      0       0.360632        0.181643        10.9302 30.3085 0.04710096      27772
    1       767638  chr1:767638:T:G T       G       20.2745 0.000365017     0       0.530985        0.229512        10.0803 18.9841 0.02069291      27772
    
    ```
# Output Files for GWAMA_Meta-Analysis


* Meta-Analysis Top Hits

    * Meta-analysis top hits all in one place. If annotate is true, these will also have RSID and nearest gene

    * Type: Summary Table

    * Format: csv

    * File Header:


    ```
    variant_id,CHR,POS,reference_allele,other_allele,eaf,beta,se,beta_95L,beta_95U,z,p-value,_-log10_p-value,q_statistic,q_p-value,i2,n_studies,n_samples,effects,PHENO,ANALYSIS,OR,OR_se,OR_95L,OR_95U,Gene,RSID
    chr10:112974337:A:G,10,112974337,G,A,0.179407,,,,,4.572918,4.88e-06,5.31143,2.193479,0.138596,0.544103,2.0,-9.0,++,T2D,AFR_EUR,1.118071,0.026645,1.065847,1.172854,TCF7L2,rs11196174
    chr10:112976855:T:C,10,112976855,C,T,0.17661,,,,,4.576326,4.8e-06,5.318488,1.65513,0.198262,0.395818,2.0,-9.0,++,T2D,AFR_EUR,1.118606,0.02675,1.066176,1.173613,TCF7L2,rs11196175
    chr10:112986821:C:G,10,112986821,G,C,0.546976,,,,,5.571755,2.58e-08,7.587792,9.51063,0.023219,0.684563,4.0,-9.0,++++,T2D,ALL,1.11404,0.021188,1.072513,1.157176,TCF7L2,rs4073980
    chr10:112986821:C:G,10,112986821,G,C,0.605747,,,,,5.352806,8.85e-08,7.052915,7.074541,0.007819,0.858648,2.0,-9.0,++,T2D,AFR_EUR,1.110165,0.021265,1.068484,1.153471,TCF7L2,rs4073980
    
    ```

* Meta-Analysis Summary Statistics
>>>>>>> a5ba8f90b2a7163dc8a0d6179b32b0c1497ed225

    * Meta-analysis summary stats

<<<<<<< HEAD
❗We HIGHLY recommend doing a STUB run to test the analysis using the `-stub` flag. This is a dry run to make sure your environment, parameters, and input_files are specified and formatted correctly.❗We also HIGHLY recommend doing a TEST run with the included test data in `$TOOLS_DIR/pmbb-nf-toolkit-gwama-meta/test_data`we have several pre-configured analyses runs with input data and fully-specified config files.

```sh
# run an exwas stub
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-gwama-meta/gwama_meta.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c GWAMA_META/gwama_meta.config \
   -stub

# run an exwas for real
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-gwama-meta/gwama_meta.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c GWAMA_META/gwama_meta.config

# resume an exwas run if it was interrupted or ran into an error
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-gwama-meta/gwama_meta.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c GWAMA_META/gwama_meta.config \
   -resume
```
# Pipeline Parameters
=======
    * Type: Summary Statistics

    * Format: txt.gz

    * File Header:


    ```
    variant_id      CHR     POS     reference_allele        other_allele    eaf     beta    se      beta_95L        beta_95U        z       p-value _-log10_p-value q_statisticq_p-value        i2      n_studies       n_samples       effects
    chr1:100000723:G:A      1       100000723       A       G       0.151861        0.004207        0.011077        -0.017504       0.025918        0.379812        0.7041  0.152366    2.172915        0.33741 0.079577        3.0     38693.0 ?+-+
    chr1:10000113:C:T       1       10000113        T       C       0.100108        -0.043836       0.024287        -0.091439       0.003767        -1.804897       0.071099   1.148135 0.0     1.0             1.0     10275.0 ???-
    chr1:100001396:G:C      1       100001396       C       G       0.129179        -0.014217       0.021629        -0.056609       0.028176        -0.657298       0.511007   0.291573 0.0     1.0             1.0     10275.0 ???-
    chr1:100002416:C:T      1       100002416       T       C       0.064698        0.034357        0.029261        -0.022994       0.091709        1.174164        0.240309   0.61923  0.0     1.0             1.0     10275.0 ???+
    
    ```

        * Parallel By: Analysis, Phenotype
# Parameters for GWAMA_Meta-Analysis

## GWAMA


* `gwama_path` (Type: File Path)

    * A path to the GWAMA executable. If using the container, use /app/GWAMA

    * Corresponding Input File: GWAMA Executable

        * Path to the GWAMA executable

        * Type: Executable

        * Format: exe
>>>>>>> a5ba8f90b2a7163dc8a0d6179b32b0c1497ed225

## Input Files for GWAMA_Meta-Analysis


* GWAMA Executable

    * Path to the GWAMA executable

    * Type: Executable

    * Format: exe

* GWAS Summary Statistics

    * Files with GWAS summary stats to be meta-analyzed. 

    * Type: Summary Statistics

    * Format: txt.gz

    * File Header:


<<<<<<< HEAD
    ```
    chromosome      base_pair_location      variant_id      other_allele    effect_allele   effect_allele_count     effect_allele_frequency missing_rate    beta    standard_error      t_statistic     variance        p_value N
    1       722858  chr1:722858:C:T C       T       59.0157 0.0010625       0       0.364206        0.178679        11.4078 31.3223 0.04151698      27772
    1       764898  chr1:764898:G:A G       A       54.5843 0.000982722     0       0.395419        0.189646        10.9943 27.8043 0.03706621      27772
    1       766367  chr1:766367:G:A G       A       57.851  0.00104153      0       0.360632        0.181643        10.9302 30.3085 0.04710096      27772
    1       767638  chr1:767638:T:G T       G       20.2745 0.000365017     0       0.530985        0.229512        10.0803 18.9841 0.02069291      27772
    
    ```
## Output Files for GWAMA_Meta-Analysis
=======
    * Whether or not to annotate results with the RSIDs and nearest genes for plotting and summary files.
## Pre-Processing
>>>>>>> a5ba8f90b2a7163dc8a0d6179b32b0c1497ed225


* `input_col_names` (Type: Map (Dictionary))

    * A Groovy config map where the keys are the required GWAMA columns and the values are the corresponding column names from your input files. The required columns for GWAMA are [MARKERNAME, EA, NEA, N, EAF, CHR, POS] whereas for the effect sizes you can have any of the following sets [OR, SE], [OR, OR_95L, OR_95U], [BETA, SE]

* `sumstats_suffix` (Type: String)

    * Summary stats are expected to be organized in the form ${launchDir}/{cohort}/Sumstats/{pheno}{suffix}. The combinations of cohort and phenotype will be constructed from the lists provided, but you will need to provide the suffix which should be the same for all input files.

    * Corresponding Input File: GWAS Summary Statistics

<<<<<<< HEAD
    * File Header:
=======
        * Files with GWAS summary stats to be meta-analyzed. 

        * Type: Summary Statistics

        * Format: txt.gz

        * File Header:
>>>>>>> a5ba8f90b2a7163dc8a0d6179b32b0c1497ed225


        ```
        chromosome      base_pair_location      variant_id      other_allele    effect_allele   effect_allele_count     effect_allele_frequency missing_rate    beta    standard_error      t_statistic     variance        p_value N
        1       722858  chr1:722858:C:T C       T       59.0157 0.0010625       0       0.364206        0.178679        11.4078 31.3223 0.04151698      27772
        1       764898  chr1:764898:G:A G       A       54.5843 0.000982722     0       0.395419        0.189646        10.9943 27.8043 0.03706621      27772
        1       766367  chr1:766367:G:A G       A       57.851  0.00104153      0       0.360632        0.181643        10.9302 30.3085 0.04710096      27772
        1       767638  chr1:767638:T:G T       G       20.2745 0.000365017     0       0.530985        0.229512        10.0803 18.9841 0.02069291      27772
        
        ```
## Workflow


* `my_python` (Type: File Path)

    * Path to the python executable to be used for python scripts - often it comes from the docker/singularity container (/opt/conda/bin/python)

* `analyses` (Type: Map (Dictionary))

<<<<<<< HEAD
    * File Header:
=======
    * Map of lists where keys are meta-analysis group nicknames and lists are groups of cohorts to include in that meta-analysis. This allows for multiple combinations of meta-analyses, for example all cohorts of one sex/ancestry, leave-one-biobank-out.

* `bin_pheno_list` (Type: List)

    * Binary phenotype list
# Configuration and Advanced Workflow Files

## Example Config File Contents (From Path)

>>>>>>> a5ba8f90b2a7163dc8a0d6179b32b0c1497ed225

```
params {
    analyses = [
        'AFR_EUR': ['AFR_ALL', 'EUR_ALL'],
        'ALL': ['AFR_ALL', 'EUR_ALL', 'EAS_ALL', 'AMR_ALL', 'SAS_ALL'],
        'ALL_M': ['AFR_M', 'EUR_M', 'EAS_M', 'AMR_M', 'SAS_M'],
        'ALL_F': ['AFR_F', 'EUR_F', 'EAS_F', 'AMR_F', 'SAS_F'],
        'Leave_EUR_Out': ['AFR_ALL', 'EAS_ALL', 'AMR_ALL', 'SAS_ALL']
    ]

    
<<<<<<< HEAD
    ```

        * Parallel By: Analysis, Phenotype
## Other Parameters for GWAMA_Meta-Analysis

### GWAMA


* `min_meta_input_EAF` (Type: Float)

    * The minimum effect allele frequency to use when munging the input summary stats. 
### Post-Processing


* `biofilter_close_dist` (Type: Float)

    * The distance in bp for something to be considered “close” vs “far” with respect to nearest gene annotation. Value is often 5E4
### Pre-Processing


* `input_col_names` (Type: Map (Dictionary))

    * A Groovy config map where the keys are the required GWAMA columns and the values are the corresponding column names from your input files. The required columns for GWAMA are [MARKERNAME, EA, NEA, N, EAF, CHR, POS] whereas for the effect sizes you can have any of the following sets [OR, SE], [OR, OR_95L, OR_95U], [BETA, SE]

* `sumstats_suffix` (Type: String)

    * Summary stats are expected to be organized in the form ${launchDir}/{cohort}/Sumstats/{pheno}{suffix}. The combinations of cohort and phenotype will be constructed from the lists provided, but you will need to provide the suffix which should be the same for all input files.

    * Corresponding Input File: GWAS Summary Statistics

        * Files with GWAS summary stats to be meta-analyzed. 

        * Type: Summary Statistics

        * Format: txt.gz

        * File Header:


        ```
        chromosome      base_pair_location      variant_id      other_allele    effect_allele   effect_allele_count     effect_allele_frequency missing_rate    beta    standard_error      t_statistic     variance        p_value N
        1       722858  chr1:722858:C:T C       T       59.0157 0.0010625       0       0.364206        0.178679        11.4078 31.3223 0.04151698      27772
        1       764898  chr1:764898:G:A G       A       54.5843 0.000982722     0       0.395419        0.189646        10.9943 27.8043 0.03706621      27772
        1       766367  chr1:766367:G:A G       A       57.851  0.00104153      0       0.360632        0.181643        10.9302 30.3085 0.04710096      27772
        1       767638  chr1:767638:T:G T       G       20.2745 0.000365017     0       0.530985        0.229512        10.0803 18.9841 0.02069291      27772
        
        ```
### Workflow


* `my_python` (Type: File Path)

    * Path to the python executable to be used for python scripts - often it comes from the docker/singularity container (/opt/conda/bin/python)

* `bin_pheno_list` (Type: List)

    * Binary phenotype list

* `quant_pheno_list` (Type: List)

    * Quantitative phenotype list
# Configuration and Advanced Workflow Files

## Example Config File Contents (From Path)


```
params {
    analyses = [
        'AFR_EUR': ['AFR_ALL', 'EUR_ALL'],
        'ALL': ['AFR_ALL', 'EUR_ALL', 'EAS_ALL', 'AMR_ALL', 'SAS_ALL'],
        'ALL_M': ['AFR_M', 'EUR_M', 'EAS_M', 'AMR_M', 'SAS_M'],
        'ALL_F': ['AFR_F', 'EUR_F', 'EAS_F', 'AMR_F', 'SAS_F'],
        'Leave_EUR_Out': ['AFR_ALL', 'EAS_ALL', 'AMR_ALL', 'SAS_ALL']
    ]

    
=======
>>>>>>> a5ba8f90b2a7163dc8a0d6179b32b0c1497ed225
    // Executables for python and GWAMA
    my_python = '/opt/conda/bin/python'
    gwama_path = '/app/GWAMA'

    // Lists of phenotypes
    bin_pheno_list =  ['T2D', 'AAA']
    quant_pheno_list = ['LDL_median', 'BMI_median']

    // Min Allele Frequency for Meta-Analysis
    min_meta_input_EAF = 0.05

    // Pre- and Post-Processing Params
    sumstats_suffix = '.saige.gz'
    p_cutoff_summarize = 0.00001

    annotate = true

    // The following arguments go with annotate=true and will be used by the biofilter_wrapper sub-workflow
    biofilter_build = '38' // can be 19 or 38
    biofilter_loki = '/path/to/data/loki.db'
    biofilter_script = '/app/biofilter.py' // Must be an executable python file
    biofilter_close_dist = 5E4

    // Column names to map:
    // Keys = GWAMA Column
    // Values = Input Column
    input_col_names = [
        'MARKERNAME' : 'variant_id',
        'EA' : 'effect_allele',
        'NEA' : 'other_allele',
        'OR' : 'odds_ratio',
        'OR_95L' : 'odds_ratio_ci_95L',
        'OR_95U' : 'odds_ratio_ci_95U',
        'BETA' : 'beta',
        'SE' : 'standard_error',
        'N' : 'N',
        'N_CASE': 'n_cases',
        'N_CTRL': 'n_controls',
        'EAF' : 'effect_allele_frequency',
        'CHR' : 'chromosome',
        'POS' : 'base_pair_location'
    ]
}
```
<<<<<<< HEAD
## Current `nextflow.config` contents


```
includeConfig 'gwama_meta.config'

profiles {
    non_docker_dev {
        process.executor = awsbatch-or-lsf-or-slurm-etc
    }

    standard {
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.container = 'katiecardone26/gwama_meta:latest'
        docker.enabled = true
    }

    cluster {
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.queue = 'epistasis_normal'
        process.memory = '15GB'
    	process.container = 'gwama_meta.sif'
        singularity.enabled = true
        singularity.runOptions = '-B /root/,/directory/,/names/'
    }

    all_of_us {
        // CHANGE EVERY TIME! These are specific for each user, see docs
        google.lifeSciences.serviceAccountEmail = service@email.gservicaaccount.com
        workDir = /path/to/workdir/ // can be gs://
        google.project = terra project id

        // These should not be changed unless you are an advanced user
        process.container = 'gcr.io/verma-pmbb-codeworks-psom-bf87/gwama_meta:latest' // GCR SAIGE docker container (static)

        // these are AoU, GCR parameters that should NOT be changed
        process.memory = '15GB' // minimum memory per process (static)
        process.executor = awsbatch-or-lsf-or-slurm-etc
        google.zone = "us-central1-a" // AoU uses central time zone (static)
        google.location = "us-central1"
        google.lifeSciences.debug = true 
        google.lifeSciences.network = "network"
        google.lifeSciences.subnetwork = "subnetwork"
        google.lifeSciences.usePrivateAddress = false
        google.lifeSciences.copyImage = "gcr.io/google.com/cloudsdktool/cloud-sdk:alpine"
        google.enableRequesterPaysBuckets = true
        // google.lifeSciences.bootDiskSize = "20.GB" // probably don't need this
    }
}

params {
    skip_postprocessing_errors = false
}

process {
    withLabel: safe_to_skip {
        errorStrategy=params.skip_postprocessing_errors ? 'ignore' : 'terminate'
    }
}



```
## Advanced Nextflow Users: Take/Emit Info

=======
## Current Dockerfile for Container/Image


```docker
FROM continuumio/miniconda3
WORKDIR /app

# biofilter version argument
ARG BIOFILTER_VERSION=2.4.3

RUN apt-get update \    
    # install packages needed to install GWAMA, biofilter, and NEAT-plots
    && apt-get install -y --no-install-recommends libz-dev g++ gcc git wget tar unzip make \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # install GWAMA
    && wget https://www.geenivaramu.ee/tools/GWAMA_v2.2.2.zip \
    && unzip GWAMA_v2.2.2.zip \
    && make \
    # install python packages needed for pipeline
    && conda install -y -n base -c conda-forge wget libtiff conda-build scipy pandas seaborn matplotlib numpy apsw sqlite \
    && conda clean --all --yes \
    # install NEAT-plots
    && git clone https://github.com/PMBB-Informatics-and-Genomics/NEAT-Plots.git \
    && mv NEAT-Plots/manhattan-plot/ /app/ \
    && conda develop /app/manhattan-plot/ \
    # install biofilter
    && wget https://github.com/RitchieLab/biofilter/releases/download/Biofilter-${BIOFILTER_VERSION}/biofilter-${BIOFILTER_VERSION}.tar.gz -O biofilter.tar.gz \
    && tar -zxvf biofilter.tar.gz --strip-components=1 -C /app \
    && /opt/conda/bin/python setup.py install \
    # make biofilter executable
    && chmod a+rx /app/biofilter.py \
    # remove biofilter tarball and NEAT-plots directory
    && rm -R biofilter.tar.gz NEAT-Plots

USER root

```
## Current `nextflow.config` contents


```
includeConfig 'gwama_meta.config'

profiles {
    non_docker_dev {
        process.executor = awsbatch-or-lsf-or-slurm-etc
    }

    standard {
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.container = 'katiecardone26/gwama_meta:latest'
        docker.enabled = true
    }

    cluster {
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.queue = 'epistasis_normal'
        process.memory = '15GB'
    	process.container = 'gwama_meta.sif'
        singularity.enabled = true
        singularity.runOptions = '-B /root/,/directory/,/names/'
    }

    all_of_us {
        // CHANGE EVERY TIME! These are specific for each user, see docs
        google.lifeSciences.serviceAccountEmail = service@email.gservicaaccount.com
        workDir = /path/to/workdir/ // can be gs://
        google.project = terra project id

        // These should not be changed unless you are an advanced user
        process.container = 'gcr.io/verma-pmbb-codeworks-psom-bf87/gwama_meta:latest' // GCR SAIGE docker container (static)

        // these are AoU, GCR parameters that should NOT be changed
        process.memory = '15GB' // minimum memory per process (static)
        process.executor = awsbatch-or-lsf-or-slurm-etc
        google.zone = "us-central1-a" // AoU uses central time zone (static)
        google.location = "us-central1"
        google.lifeSciences.debug = true 
        google.lifeSciences.network = "network"
        google.lifeSciences.subnetwork = "subnetwork"
        google.lifeSciences.usePrivateAddress = false
        google.lifeSciences.copyImage = "gcr.io/google.com/cloudsdktool/cloud-sdk:alpine"
        google.enableRequesterPaysBuckets = true
        // google.lifeSciences.bootDiskSize = "20.GB" // probably don't need this
    }
}

params {
    skip_postprocessing_errors = false
}

process {
    withLabel: safe_to_skip {
        errorStrategy=params.skip_postprocessing_errors ? 'ignore' : 'terminate'
    }
}



```
## Advanced Nextflow Users: Take/Emit Info

>>>>>>> a5ba8f90b2a7163dc8a0d6179b32b0c1497ed225
### Input Channel (take) Description


a Channel of 3-part tuples with (cohort, phenotype, summary stats file path). If your summary stats are NOT organized as ${launchDir}/{cohort}/Sumstats/{pheno}{suffix}, then you can set up your own input channel to correspond to your input files and pass it to the GWAMA_META named workflow.
### Output Channel (emit) Description
<<<<<<< HEAD
=======


a Channel of 3-part tuples with (analysis, phenotype, summary stats) to potentially be passed to other workflows. Think of “analysis” as a population; it’s a set of cohorts used in this meta-analysis
# Detailed Pipeline Steps


from pathlib import Path
>>>>>>> a5ba8f90b2a7163dc8a0d6179b32b0c1497ed225

detailed_steps_file = Path("Markdowns/Pipeline_Detailed_Steps.md")

# Write the detailed steps content to a separate file
detailed_steps_file

# Detailed Steps for Runnning One of our Pipelines

Note: test data were obtained from the [SAIGE github repo](https://github.com/saigegit/SAIGE).

## Part I: Setup
1. Start your own tools directory and go there. You may do this in your project analysis directory, but it often makes sense to clone into a general `tools` location

```sh
# Make a directory to clone the pipeline into
TOOLS_DIR="/path/to/tools/directory"
mkdir $TOOLS_DIR
cd $TOOLS_DIR
```

2. Download the source code by cloning from git

```sh
git clone https://github.com/PMBB-Informatics-and-Genomics/pmbb-nf-toolkit-saige-family.git
cd ${TOOLS_DIR}/pmbb-nf-toolkit-saige-family/
```

3. Build the `saige.sif` singularity image
- you may call the image whatever you like, and store it wherever you like. Just make sure you specify the name in `nextflow.conf`
- this does NOT have to be done for every saige-based analysis, but it is good practice to re-build every so often as we update regularly. 

```sh
cd ${TOOLS_DIR}/pmbb-nf-toolkit-saige-family/
singularity build saige.sif docker://pennbiobank/saige:latest
```

## Part II: Configure your run

1. Make a separate analysis/run/working directory.
   - The quickest way to get started, is to run the analysis in the folder the pipeline is run. However, subsequent analyses will over-write results from previous analyses. 
   - ❗This step is optional, but We Highly recommend making a  `tools` directory separate from your `run` directory. The only items that need to be in the run directory are the `nextflow.conf` file and the `${workflow}.conf` file.

```sh
WDIR="/path/to/analysis/run1"
mkdir -p 
cd $WDIR
```

2. Fill out the `nextflow.config` file for your system.
    - See [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html) for information on how to configure this file. An example can be found on our GitHub: [Nextflow Config](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/Example_Configs/nextflow.config).
    - ❗IMPORTANTLY, you must configure a user-defined profile for your run environments (local, docker, saige, cluster, etc.). If multiple profiles are specified, run with a specific profile using `nextflow run -profile ${MY_PROFILE}`.
    - For singularity, The profile's attribute `process.container` should be set to `'/path/to/saige.sif'` (replace `/path/to` with the location where you built the image above). See [Nextflow Executor Information](https://www.nextflow.io/docs/latest/executor.html) for more details.
    - ⚠️As this file remains mostly unchanged for your system, We recommend storing this file in the `tools/pipeline` directory and symlinking it to your run directory.

3. Create a pipeline-specific `.config` file specifying your run parameters and input files. See Below for workflow-specific parameters and what they mean.
   - Everything in here can be configured in `nextflow.config`, however we find it easier to separate the system-level profiles from the individual run parameters. 
   - Examples can be found in our Pipeline-Specific [Example Config Files](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/Example_Configs/).
   - you can compartamentalize your config file as much as you like by passing 
   - There are 2 ways to specify the config file during a run:
      - with the `-c` option on the command line: `nextflow run -c /path/to/workflow.conf`
      - in the `nextflow.conf`: at the top of the file add: `includeConfig '/path/to/workflow.conf'` 

## Part III: Run your analysis

- ❗We HIGHLY recommend doing a STUB run to test the analysis using the `-stub` flag. This is a dry run to make sure your environment, parameters, and input_files are specified and formatted correctly. 
- ❗We HIGHLY recommend doing a test run with the included test data in `${TOOLS_DIR}/pmbb-nf-toolkit-saige-family/test_data`
- in the `test_data/` directory for each pipeline, we have several pre-configured analyses runs with input data and fully-specified config files.

```sh
# run an exwas stub
nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -c /path/to/run1/exwas.conf -stub
# run an exwas for real
nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -c /path/to/run1/exwas.conf
# resume an exwas run if it was interrupted or ran into an error
nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -c /path/to/run1/exwas.conf -resume
```

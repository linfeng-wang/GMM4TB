# Genomic_data_analysis
LONDONSCHOOL OF HYGIENE AND TROPICAL MEDICINE

Multi-strain infection model
Models for strain calling and drug resistance profiling

## Model introduction
**Mixed-strain infection in tuberculosis** (TB) is a fascinating and complex phenomenon that challenges our understanding of this ancient disease. Unlike single-strain infections, where a patient harbors a single type of TB bacterium, mixed-strain infections involve the coexistence of multiple distinct strains of the Mycobacterium tuberculosis bacteria within the same host. These strains can vary in their genetic makeup, drug susceptibility, and virulence, potentially leading to diverse clinical outcomes. Studying mixed-strain infections not only sheds light on TB's genetic diversity but also has significant implications for diagnosis, treatment, and control efforts, making it a critical area of research in the fight against this global health threat.

**The Gaussian Mixture Model** (GMM) serves as a powerful tool in disentangling Mixed-Strain Infections (MSI) within tuberculosis cases, allowing for the separation of different lineages and providing valuable insights into drug resistance profiles associated with each lineage. By leveraging the distinct genetic signatures of TB strains, GMM can effectively cluster the heterogeneous mix of strains present in an MSI, helping researchers and clinicians identify individual lineages within the infection.

Furthermore, GMM's ability to unravel these lineages provides a critical advantage in understanding the drug resistance patterns within each lineage. Different TB lineages often exhibit unique resistance profiles due to variations in their genetic makeup. GMM's precision in distinguishing these lineages facilitates the targeted assessment of drug resistance, enabling healthcare professionals to tailor treatment strategies more effectively. This disentanglement of MSI and separation of drug resistance profiles at the lineage level represent a significant advancement in our battle against drug-resistant tuberculosis, offering a promising avenue for improved patient care and public health management.

*Path: /Executable/gmm_model_tbp/main_multi.py
*
required package: 
- tb-profiler(https://github.com/jodyphelan/TBProfiler)
- sckitlearn
- scipy
- uuid
- pathlib
- numpy
- pandas
- plotly
  
Python 3.10.0  

**How to run in terminal:**
Download `Executable/gmm_model_tbp` folder and install above packages. Then run the following command in terminal:

Running using **vcf** as input:

```
python **<Path>**/main_multi.py -vcf <Path_to_vcf.vcf.gz> -g -m -o <output_name> -op <output_path>
```

Running using **json** file as input:

```
python **<Path>**/main_multi.py -json <Path_to_json.json> -g -m -o <output_name> -op <output_path>
```


**-vcf**: path to vcf file

**-g**: enable graphic output

**-m**: enable 2+ mixes MSI detection

**-v**: enable vcf file input (if not enabled, require tb-profiler output file in json format) 

**-json**: path to json file

**-o**: output file name

**-op**: output path


*example output: Executable_eval/results/insilico_mix_new/mixture_results*

## Stand alone prediction method
A stand alone prediction method script has also been produced to predict for any mixed infections or gene reads not only limited for use in TB.
`Genomic_data_analysis/Executable/gmm_model_tbp/stand_alone.py`
Here a input the for of the example file can be used `Executable/gmm_model_tbp/cache/stand_alone_test.txt`
in the form of:

| SNP_Name | Alt allele Fraction |
|---|---|
| rs1001 | 0.25 |
| rs1001 | 0.25 |
| rs1002 | 0.32 |
| rs1002 | 0.32 |
| rs1003 | 0.23 |

No specific column names are needed.

Use example 
python stand_alone.py -i <input.csv> -c <number of mixes expected> -o <ouput file path> -g <graphic output>

Example output file: `Genomic_data_analysis/Executable/gmm_model_tbp/cache/stand_alone_test_out.csv`
| SNP_Name | Fraction | cluster_no |
|---|---|---|
| rs1009 | 0.26 | 1 |
| rs1010 | 0.31 | 1 |
| rs1010 | 0.31 | 1 |
| rs1011 | 0.74 | 2 |
| rs1011 | 0.74 | 2 |

Exampe output graphics:`Executable/gmm_model_tbp/cache/stand_alone_test_out_gmm_fig.png`

## Model function
![model_funciton](img/gmm_process.png)

The Gaussian Mixture Models (GMMs) can accommodate any number of mixtures through the 'multi' option. The results for each sample encompass various parameters, such as the count of mixtures, their attributes (mean and standard deviation (SD)), confidence intervals around the mixtures (±1 SD from the mean), and the proportions in which they are mixed. Each Single Nucleotide Polymorphism (SNP) is allocated within these Gaussian distributions, allowing for deductions about the assigned strain (component) and its associated confidence, typically represented as a probability. This process aids in the identification of individual strains and their corresponding drug resistance patterns.

Profiling drug resistance and lineages from WGS data can be employed for clinical and infection control purposes, as exemplified by the TB-Profiler tool. Nevertheless, existing software tools often fall short when it comes to effectively separating distinct SNPs on different strains within a Mixed Strain Infection (MSI), potentially limiting the precision of profiling. To address this challenge, Gaussian Mixture Models (GMMs) were individually constructed for each sample using Scikit-learn, and subsequently applied to assess the ratios of alternative alleles relative to the total allele count across the SNPs present in a variant calling file (in vcf format).

<h2>Folder directory</h2>
<h3>Artificial_mixture_creation</h3>
Script used to create in silico artifical Multi-strain infection sample order of runing (strain_analysis folder)
    <ol>
    <li>seqtk.sh</li>
    <li>f2m_parallel.sh</li>
    <li>freebayes_parallel.sh</li>
    <li>jsonfile_gen.sh</li>
    </ol>
    <h3> Don't forget to change paths for each file before running</h3>

<h3>in silico aritifical created testing sample</h3>
Some artifically generated mixture sample. Can be used as example file for running model

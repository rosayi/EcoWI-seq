# EcoWI-seq
EcoWI-seq is an enzymatic next-generation sequencing method for genomic mapping of DNA phosphorothioate (PT) modification at single base resolution. This script is developed for analysis of PT modification data generated by EcoWI-seq and for calling of PT modification sites(base positions) in the provided genome. For details, please refer to our preprint:(update later). <br>

Script author: Weiwei Yang (New England Biolabs, wyang@neb.com) <br>

## Before starting <br>
The following programs and python modules are required to be installed: <br>
- samtools (v1.9) <br>
- bedtools (v2.29.2) <br>
- pandas (python module) <br>
- seaborn (python module) <br>
- matplotlib (python module)<br>

## Command-line usage <br>
`python EcoWI_seq.py --input input.bam --genome genome.fasta --read_cutoff 50 --score_cutoff 20 —-window 100`<br>

This script performs modification score calculation to quantify the PT modification level for every base position in the provided genome(s). The calling of modified sites then is processed on the provided cutoff of the modification score (and the others e.g. read_cutoff). <br>

**Input:**<br>
`--input` or `-i`: bam file <br>

`--genome` or `-g`: genome of reference (.fasta or .fa) used for PT modification mapping. Multi-FASTA format is allowed. <br>

**Arguments:** <br>
`--read_cutoff` or `-r`: int, minimum read number for a base to be considered as modified in PT site calling (default=50 if no value is provided) <br>

`--score_cutoff` or `-s`: float, minimum modification score value for a base to be considered as modified in PT site calling (default=0 if no value is provided) <br> 

`--window` or `-w`: int, length of window for calculate modification score, i.e. `-w 50` means the modification score is calculated within +/- 50bp of a base position (default=50 if no value is provided) <br>

**Optional argument:** <br>
`--plotcurve`: plot score_cutoff (x axis) over number of PT modification sites (y axis). This option helps to determine the appropriate score cutoff (visual observation of the deviation from the tangent point). <br>

## Outputs: 
`<input>.modification`: a tab separated file containig 12 fields: <br>
| column | field | description|
| ----------- | ----------- | ---------|
| 1 | seq | genome name |
| 2 | source | bam file source |
|	3	| feature	| motif on the ref sequence: GAAC or GTTC or other |
|	4	| start |	position of the base, 1-based |
|	5	| end	|	position of the base, 1-based |
|	6	| motif_pos	|	position of the corresponding GAAC/GTTC motif or 0 |
|	7	| strand |	+ |
| 8	| 5’ end reads on + strand	|	number of 5’ end reads (reads start at the base position) mapped on + strand |
|	9	| 5’ end reads on - strand	|	number of 5’ end reads (reads end at the base position) mapped on - strand |
|	10	| total reads on + strand	| number of total reads cover this position on + strand |
|	11	| total reads on - strand	|	number of total reads cover this position on - strand |
|	12	| modification score	|	modification score = 5’ end reads / median 5’ end reads within +/- 50bp (or custom value given by the —-window argument) |


`<genome>.allPTsites`: a tab separated file containing the positions of all possible sites and their corresponding motifs (GAAC or GTTC) in the provided genome <br>
e.g. 17  GAAC <br>

`<input>.modification_score.png`: a scatter plot showing the modification scores, every dot represents a base position, plot is grouped by two categories: GAAC/GTTC or other (whether the base is in the GAAC/GTTC motif or not) <br>

`<input>.sites`: the modified sites called by the analysis, containing two tab separated fields:
| column | field | description|
| ----------- | ----------- | ---------|
| 1 | pos | position of the modified base (1-based) |
| 2 | motif | GAAC/GTTC/other |

**Optional output:** <br>
`<input>.curve.png`: plot score_cutoff (x axis) over number of modified sites (y axis). This option helps to determine the appropriate score cutoff (visual observation of the deviation from the tangent point). <br>

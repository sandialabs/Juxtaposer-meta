# Juxtaposer-meta
Juxtaposer-meta aims to detect sharing of mobile DNAs between species in metagenomic samples. Juxtaposer-meta also makes basic improvements to Juxtaposer (https://github.com/sandialabs/juxtaposer), for example recovering far more recombinant reads.

## AUTHORS

- Hannah McClain
- Kelly Williams ([kpwilli@sandia.gov](mailto:kpwilli@sandia.gov))
- Repo maintained by: Ellis Torrance ([eltorra@sandia.gov](mailto:eltorra@sandia.gov))

## CITATIONS

If you use Juxtaposer-meta in your research, please cite:

> Schoeniger JS, Hudson CM, Bent ZW, Sinha A, Williams KP. 2016. *Experimental single-strain mobilomics reveals events that shape pathogen emergence.* *Nucleic Acids Research* **44**:6830-9. doi: [10.1093/nar/gkw601](https://doi.org/10.1093/nar/gkw601)

## SOFTWARE DEPENDENCIES

- `bowtie2` and `bowtie2-build`, in path
- `blastn` and `makeblastdb`, in path
- `hmmscan`, in path
- `prodigal`: binary is included, but if the binary fails on your machine, replace it with a fresh install

Instructions for installing these dependencies in a conda environment are available in [INSTALLATION](#installation) below.

## INSTALLATION

1. Download the GitHub repo from <https://github.com/sandialabs/Juxtaposer-meta> or clone it with git:

```bash
git clone https://github.com/sandialabs/Juxtaposer-meta.git
```

2. Make Juxtaposer-meta files callable system-wide by adding the following to your bash profile. Then restart the terminal session or use `source` to reload your bash profile.

```bash
export PATH=<absolute path to juxtaposer-meta folder>:$PATH
```

3. Install Juxtaposer-meta dependencies as below.

### (OPTIONAL) Creating a Juxtaposer-meta Conda environment

1. Create the conda environment:

```bash
conda create --name Juxtaposer-meta -y
```

2. Activate the conda environment  
   *(Note: you will need to activate the conda environment every time before you run Juxtaposer-meta.)*

```bash
conda activate Juxtaposer-meta
```

3. Install dependencies:

```bash
conda config --add channels bioconda && conda install -y bowtie2 && conda install -y blast && conda install -y hmmer && conda install -y pftools && conda install -y bedtools
```

## INPUT FILES

Juxtaposer-meta requires a (meta)-genome assembly FASTA file and one or more sequencing read files in FASTQ format (acceptable in gzipped form with a `.gz` suffix). We recommend pre-treatment of the reads by quality filtering and primer trimming, for example using `bbduk`.

## RUNNING JUXTAPOSER-META

Aside from the required (meta)-genome assembly FASTA file and one or more read FASTQ files, the user can specify an output directory; otherwise all outputs are written to the current directory.

Additionally, if the user has already run Bowtie, they may provide a path to an alignment file to avoid re-running Bowtie. Other options are shown by running:

```bash
perl PATH/juxtaposer-meta.pl -help
```

There are two main Juxtaposer calls:

### Run Bowtie2 freshly

```bash
perl juxtaposer-meta.pl -r <read files> -g <assembly file>
```

### Parse an existing Bowtie2 file (SAM or BAM) for candidate reads

```bash
perl juxtaposer-meta.pl -r <read files> -g <assembly file> -a <bowtie file>
```

## OUTPUT FILES

Juxtaposer-meta has two primary output files: `juxtas.txt` and `networks.json`.

### `juxtas.txt`

| Field | New label | Old label | Notes |
|---|---|---|---|
| 1 | `dnaLo` | `dnaLo` | Lexically lower contig |
| 2 | `CoordLo` | `CoordLo` | Overlap sequence coordinate |
| 3 | `dnaHi` | `dnaHi` | Lexically higher contig |
| 4 | `CoordHi` | `CoordHi` | Overlap sequence coordinate |
| 5 | `config` | `config` | `++` / `--` / `+-` / `-+` |
| 6 | `call` | `call` | If `dnaLo = dnaHi`: `CircleJunction`, `Scar`, `Palindrome` |
| 7 | `overlap` | `overlap` | Length of shared sequence between Lo and Hi hits |
| 8 | `origin` | `origin` | If `dnaLo = dnaHi`: `+1`, `-1` |
| 9 | `idLo` | `idLo` | Percent identity (BLAST) of Lo hit |
| 10 | `idHi` | `idHi` | Percent identity (BLAST) of Hi hit |
| 11 | `read_termini` | `tnLo` | Designates if each extreme end of the read is included in the juxtaposition: `1,1`, `1,0`, `0,1`, `0,0` |
| 12 | `dna_termini` | `tnOffsetLo` | Designates if the hits reach the terminal end of the DNA (contig): `1,1`, `1,0`, `0,1`, `0,0` |
| 13 | `local_count` | `tnHi` | Total number of juxtapositions in the local group |
| 14 | `local_rank` | `tnOffsetHi` | Rank within the local group, according to `idLo + idHi` |
| 15 | `sampleRead` | `sampleRead` | Recombinant read ID |
| 16 | `overlapSeq` | `overlapSeq` | Sequence within the recombinant read shared by the Lo and Hi hits; length should equal field 7 |
| 17 | `readseq` | `readseq` | Recombinant read sequence |
| 18 | `seqLo` | `seqLo` | Hit sequence at `coordLo` |
| 19 | `seqHi` | `seqHi` | Hit sequence at `coordHi` |
| 20 | `count` | `count` | Number of recombinant reads supporting the juxtaposition |
| 21 | `sites` | `sites` | Number of other rejected juxtapositions supported by the reads |
| 22 | `best_local` | `Mob` | Unique juxtaposition ID of the best local group |
| 23 | `Mob` | `pctMob` | Mobile gene ID and its linkage (if any), if a mobile gene was found near or within a putative network associated with the juxtaposition |
| 24 | `network` | `shortDistMob` | Network serial number; can be used as a key in `networks.json` |

#### Notes on selected fields

- **Field 11: `read_termini`**
  - `1,1`: read ends align to `dnaLo` and `dnaHi`
  - `1,0`: read end aligned on `dnaLo` but not `dnaHi`
  - `0,1`: read end did not align on `dnaLo` but did on `dnaHi`
  - `0,0`: neither read end aligned on `dnaLo` or `dnaHi`

- **Field 12: `dna_termini`**
  - `1,1`: `coordLo` is at a terminus of `dnaLo` and `coordHi` is at a terminus of `dnaHi`
  - `1,0`: `coordLo` is at a terminus of `dnaLo` and `coordHi` is contig-internal
  - `0,1`: `coordLo` is contig-internal and `coordHi` is at a terminus of `dnaHi`
  - `0,0`: `coordLo` and `coordHi` are contig-internal

- **Field 13: `local_count`**  
  A local group is a set of juxtapositions sharing the same contigs, an identical coordinate, and an unshared coordinate within 5 bp. This happens due to read errors or mutations, which may differ from the assembly sequence.

### `networks.json`

The second output file, `networks.json`, provides more information for the networks listed in column 24 of `juxtas.txt`. The keys are identical to those provided in `juxtas.txt` for ease of lookup.

For each network, the following information is provided:

| First key | Second key | Description |
|---|---|---|
| `notes` |  | If any applicable warnings were issued against the network, they are listed here |
|  |  | If an all-broken configuration is detected: |
|  |  | - An `attP` site cannot be established through secondary cycle, and a third cycle was run |
|  |  | - After the 3rd cycle the `attP` prediction was still incomplete |
| `network` | `Network` | Network structure represented by a Python dictionary, where primary juxtapositions are keys to a dictionary of secondary juxtapositions |
|  | `Co-orientations` | Juxtapositions listed in order of network appearance with their co-oriented mobile code |
|  | `Predicted cross over sites` | (Type 1) Predicted genomic `att` sites |
|  | `Predicted recombinant sites` | (Type 2) Predicted recombinant sites after third cycle |
|  | `All broken` | Boolean; `True` if type 2, else `False` |
|  | `Predicted near-terminal` | (Type 1) Predicted genomic `att` sites near the end of a contig |
|  | `Predicted cross over and near-terminal sites` | Union of the two sets |
|  | `Crossover sequence` | Consensus crossover sequence of the network |
|  | `Network-ID` | Unique network ID found in `juxtas.txt`; equal to the key of the network data structure |
| `stats` | `Overlap` | Lists the following information on each parameter: length, sum, average, variance, standard deviation, range, median(s), data, and outliers (by IQR) |
|  | `Coverage (Counts)` |  |
|  | `Local` |  |
|  | `Primary count` |  |
|  | `Secondary count` |  |
|  | `Tertiary count` | Type 2 only |
|  | `Unique contigs` | List of unique contigs and the unique coordinates on those contigs |
|  | `Accounted overlap sequence outliers` | Counts of near-terminal sites with an overlap length that is an outlier |
| `model` | `Decision tree` | Boolean based on the number of primary, secondary, and tertiary (if type 2) connections; unique decision trees for type 1 and type 2 models |
|  | `Prediction` | Boolean based on value of linear regression model, `f(z)` |
|  | `Probability` | Outcome of linear regression model, `f(z)`, expressed as a percentage |
|  | `Model Type` | `1` or `2` |

## WORK FLOW

| Step | Criteria / Threshold | Notes |
|---|---|---|
| Bowtie2 | `--min-score L,-16.6,-0.073` | Reads mapping imperfectly to assembly |
| BLASTN | Penalty: `-5`, Percent identity `>= 95%` | Penalty limits extension past identical region |
| Hit pair filtering | `Overlap >= 0`, unshared sequence ends `>= 8 bp` | Remove nested configurations and short hits |
| Best site selection | Based on length, then percent identity |  |
| Network | A group of juxtapositions sharing the same mobile gene and overlap sequence. Must be one of two network types. | 1. Genomic: adjacent juxtapositions on the same contig on either side of overlap sequence<br>2. Broken: all juxtapositions occur at contig ends |

# RCMAP
Residue Conservation in Multiple Alignment of Proteins

## Installation
### Linux
1. Download and unzip master branch zip file:
    ```
    wget https://github.com/fplewniak/RCMAP/archive/master.zip
    unzip master.zip
    ```  

2. Install the RCMAP Python package and requirements with pip:
 
    <code> pip install RCMAP-master/. </code>
    
### Windows 10
1. Download and unzip master branch zip file.

2. Open the unzipped <code>RCMAP-master</code> directory in a new window

3. Open a Windows PowerShell in the <code>RCMAP-master</code> directory: 
In Windows explorer, type <code>Ctrl-L</code> to select the address bar of 
the <code>RCMAP-master</code> window then type <code>powershell</code> and
validate.

4. In the PowerShell window, install the RCMAP Python package and its 
requirements with pip:

    <code> pip install RCMAP-master/. </code>
    
    Note: python and pip must be in your path
    
## Usage
The RCMAP package provides the <code>evaluate_seq</code> command:
```
usage: evaluate_seq [-h] --file File --seqeval string [string ...] [--positions POSITIONS [POSITIONS ...]]

    required arguments:
      --file File           input multiple protein sequence alignment file containing reference
                            sequences and sequences to evaluate
      --seqeval string [string ...]
                            names of the sequences to evaluate
                            
    optional arguments:
      -h, --help            show help message and exit
      --positions POSITIONS [POSITIONS ...]
                            list of position ranges to examine specified as start and end positions 
                            separated by a ':'. If the start position is not specified, then the
                            range will start at 1. If the end position is not specified, then the
                            range will end at the last position of the alignment. 
                            If no range is specified the whole alignment will be examined.
      --gaps                toggles accounting for gaps when computing information content
      --strict              toggles strict evaluation, an aminoacid is considered compatible only if it has been
                            in at least one of the reference sequences
      --min_info FLOAT      the minimum information content at a given position in the reference alignment required
                            to display the position
```
**Example:** `evaluate_seq --file ArsM_aln.faa --seqeval WP_045226361.1 Q969Z2 --positions 50:70 115:125 200:210 --gaps
`

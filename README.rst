RCMAP
=====

Residue Conservation in Multiple Alignment of Proteins (RCMAP) is a Python package to help manual annotation of protein
sequences by comparing them to a multiple alignment of reference sequences belonging to a functional family.

The RCMAP package provides the shell command ``evaluate_seq`` whose input is a multiple alignment file in FastA format,
containing reference sequences and one or more unknown sequences to annotate. It then displays for each unknown sequence
whether it is consistent at every user-specified position with the aminoacid conservation profile of the reference
sequences.

Installation
------------
Linux
^^^^^
1. Download and unzip master branch zip file::

    wget https://github.com/fplewniak/RCMAP/archive/master.zip
    unzip master.zip

2. Install the RCMAP Python package and requirements with pip::

    pip install RCMAP-master/.

Windows
^^^^^^^
1. Download and unzip master branch zip file.

2. Open the unzipped ``RCMAP-master`` directory in a new window

3. Open a Windows PowerShell in the ``RCMAP-master`` directory. (in Windows explorer, type ``Ctrl-L`` to select the address bar of the ``RCMAP-master`` window then type ``powershell`` and validate)

4. In the PowerShell window, install the RCMAP Python package and its requirements with pip::

    pip install RCMAP-master/.

Note: python and pip must be in your path

Usage
-----
The RCMAP package provides the ``evaluate_seq`` command::

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
          --method              calculation method of the background entropy for the information
          --window              number of positions to calculate the average of information, must be odd
          --window_method       calculation method of the weights of positions to calculate the information of a position,
                                using a window



Examples
^^^^^^^^
Basic example
"""""""""""""
::

    evaluate_seq --file ArsM_aln.faa --seqeval WP_045226361.1 Q969Z2 --positions 200:210 --gaps

    WP_045226361.1 :  200 : S :   True : 4.39 :   {'S'}    {'S': 6}
    WP_045226361.1 :  201 : N :   True : 4.39 :   {'N'}    {'N': 6}
    WP_045226361.1 :  202 : C :   True : 4.39 :   {'C'}    {'C': 6}
    WP_045226361.1 :  203 : - :   True : 4.39 :   set()    {'-': 6}
    WP_045226361.1 :  204 : - :   True : 4.39 :   set()    {'-': 6}
    WP_045226361.1 :  205 : - :   True : 4.39 :   set()    {'-': 6}
    WP_045226361.1 :  206 : - :   True : 4.39 :   set()    {'-': 6}
    WP_045226361.1 :  207 : V :   True : 4.39 :   {'V'}    {'V': 6}
    WP_045226361.1 :  208 : L :   True : 3.47 : Hydrophobic {'C': 4, 'I': 2}
    WP_045226361.1 :  209 : N :   True : 4.39 :   {'N'}    {'N': 6}
    WP_045226361.1 :  210 : L :   True : 4.39 :   {'L'}    {'L': 6}

     WP_045226361.1
     Number of True  11  : Information True  47.4
     Number of False  0  : Information False  0


    Q969Z2 :  200 : S :   True : 4.39 :   {'S'}    {'S': 6}
    Q969Z2 :  201 : D :  False : 4.39 :   {'N'}    {'N': 6}
    Q969Z2 :  202 : I :  False : 4.39 :   {'C'}    {'C': 6}
    Q969Z2 :  203 : P :  False : 4.39 :   set()    {'-': 6}
    Q969Z2 :  204 : F :  False : 4.39 :   set()    {'-': 6}
    Q969Z2 :  205 : G :  False : 4.39 :   set()    {'-': 6}
    Q969Z2 :  206 : K :  False : 4.39 :   set()    {'-': 6}
    Q969Z2 :  207 : K :  False : 4.39 :   {'V'}    {'V': 6}
    Q969Z2 :  208 : F :   True : 3.47 : Hydrophobic {'C': 4, 'I': 2}
    Q969Z2 :  209 : K :  False : 4.39 :   {'N'}    {'N': 6}
    Q969Z2 :  210 : L :   True : 4.39 :   {'L'}    {'L': 6}

     Q969Z2
     Number of True  3  : Information True  12.26
     Number of False  8  : Information False  35.14

All positions between 200 and 210 in the WP_045226361.1 sequence are consistent with the aminoacid observed in the
reference sequences shown in the last two columns. On the other hand, a majority of positions are not compatible with
the reference conservation profile in Q969Z2. Strictly conserved aminoacids at positions 201, 202 207 and 209 are not
conserved in this sequence, and it has an insertion from 203 to 206.


Raw information content accounting for gaps
"""""""""""""""""""""""""""""""""""""""""""
::

    evaluate_seq --file ArsM_aln.faa --seqeval WP_045226361.1 Q969Z2 --positions 50:70 115:125 200:210 --gaps
Displays compatibility at positions from 50 to 70, 115 to 125 and 200 to 210 of sequences WP_045226361.1 and Q969Z2
with the reference alignment in ArsM_aln.faa. Gaps are taken into account when computing information content.

Smoothed information content without gaps
"""""""""""""""""""""""""""""""""""""""""""
::

    evaluate_seq --file ArsM_aln.faa --seqeval WP_045226361.1 Q969Z2 --positions :10 20 200: --window_method hamming --window 5
Displays compatibility at positions from 1 to 10, at 20 and 200 to end of sequences WP_045226361.1 and Q969Z2
with the reference alignment in ArsM_aln.faa. Gaps are not taken into account. Information content
along the aligment is smoothed over a sliding window weighted using the Hamming method.
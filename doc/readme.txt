
To start
  Slr -seqfile alignment.paml -treefile tree_file

E.g. (In Examples/bglobin)
  Slr -seqfile bglobin.paml -treefile bglobin.trees


Options
-------
  Options may be specified in two ways, in a control file or on the
command line. Options given on the command line override those given 
in a control file, which in turn override the default options. The
format for calling Slr is:

  Slr [control_file] [-option value ...]

  The name of a control file can be given as the first argument to the
program, otherwise Slr expects "-option value" pairs, where
"option" is the name of the option to set and "value" what to set it
to. If no control file is given, Slr attempts to read options from
the default control file "slr.ctl".

Examples:
Read options from default control file or to use default
options.
  Slr
Read options from the file "myoptions.ctl"
  Slr myoptions.ctl
Same, but use different tree
  Slr myoptions.ctl -treefile newtree

Control file
------------
  The format of the control file is option, value pairs. The options
and allowed values being the same those that may be specified on the
command line. Each option should be on a separate line.
option1: value1
option2: value2
 etc....
  If a line contains a '#', then all text after the symbol is ignored.

E.g.
  # Control file for beta-globin example
  seqfile: bglobin.paml
  treefile: bglobin.trees
  outfile: out.res
  positive_only: 0



Allowed options
---------------

Default values for each option are given in square brackets.

seqfile [incodon]
  File from which to read alignment of codon sequences. The file
  should be in PAML format.

treefile [intree]
  File from which tree should be read. The tree should be in Nexus
  format

outfile	[slr.res]
  File to which results are written. If the file already exists, it will
  be overwritten.

reoptimize [1]
  Should the branch lengths, omega and kappa be reoptimized?
  0 - no
  1 - yes.

kappa [2.0]
  Value for kappa. If 'reoptimize' is specified, the value
  given will be used as am initial estimate,

omega [0.1]
  Value for omega (dN/dS). If 'reoptimize' is specified, the value
  given will be used as an initial estimate.

codonf [0]
  How codon frequencies are estimated:
    0: F61/F60	Estimates used are the empirical frequencies from the
  data.
    1: F3x4	The frequencies of nucleotides at each codon position
  are estimated from the data and then multiplied together to get the
  frequency of observing a given codon. The frequency of stop codons is
  set to zero, and all other frequencies scaled appropriately.
    2: F1x4	Nucleotide frequencies are estimated from the data
  (not taking into account at which position in the codon it occurs).
  The nucleotide frequencies are multiplied together to get the frequency 
  of observing and then corrected for stop codons.
  
positive_only [0]
  If only positively selected sites are of interest, set this to "1".
  Calculation will be slightly faster, but information about sites under
  purifying selection is lost. 

gencocde [universal]
  Which genetic code to use when determining whether a given mutation
  is synonymous or nonsynonymous. Currently only "universal" and
  "mammalian" mitochondrial are supported.

<B>Model Estimation </B>

A process to estimate amino acid substitution model for a given set of alignments.  <br>
Amino acid replacement matrices play an important role in analyzing protein sequences, especially inferring phylogenies. 
Our program estimates a substitution model that is likely represent of evolutionary process of a given set of alignments the best. 
Using an appropriate substitution model will enhance the accuracy of maximum likelihood tree inference. 
The program is available for Linux operating system.
<br><br><b>Notes:</b><br>
<br>1. The program requires software IQ-TREE (http://www.iqtree.org/) and path to IQTREE folder must be entered at line 42.
<br>2. The data size must large enough to gain a reliable substitution model.
<br><br><b>Commands:</b><br>
<br><i>
python Est.py -d [alignment_directory] -t [number_of_thread] -i [maximum_iteration] -mset [set_of_models] -th [correlation_threshold] 
</i> 
<br>The estimated substitution models will be written to the file <b>[alignment_directory].Mx.PAML</b> where x is the value of the iteration. 
The file contains an exchangeability Matrix and a frequency vector.
<br><i>-d</i>&nbsp;&nbsp;&nbsp;&nbsp;The path to the directory contains all training alignments in Phylip format.
<br><i>-t</i>&nbsp;&nbsp;&nbsp;&nbsp; The number of threads to be divided (optional – defaut: 1).
<br><i>-i</i>&nbsp;&nbsp;&nbsp;&nbsp; The maximum number of iterations can be executed (optional – default: 3).
<br><i>-mset</i>&nbsp;&nbsp;A set of substitution models separated by ‘,’ used as initial models to optimize (optional)
<br>&ensp;&ensp;&ensp;&nbsp;&nbsp; The default model set: FLU,HIVb,HIVw,JTT,WAG,LG.
<br><i>-th</i>&nbsp;&nbsp;&nbsp;&nbsp;	The correlation threshold, so the estimation process will be stopped if the correlation between newly estimated models and the one estimated in the previous iteration greater than the threshold
<br><br><b><i>Example:</b> python Est.py -f Example<i>


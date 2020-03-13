<B>Model Estimation </B>
<br>
A script to estimate an amino acid substitution model from a set of alignments using
the QMaker program (https://www.biorxiv.org/content/10.1101/2020.02.20.958819v1).
<br><br><b>Notes:</b> The script requires the IQ-TREE package (http://www.iqtree.org/). You have to declare the path to your installed IQ-TREE folder at line 34th of the script.
<br><br><b>Commands:</b><br>
<br><i>
python Est.py -d [alignment_directory] -t [number_of_thread] -mset [set_of_models]
</i> 
<br><br>
<b>Where</b>
<br><i>-d</i>&nbsp;&nbsp;&nbsp;&nbsp;&ensp;&ensp;The path to a directory containing all training alignments in Phylip format (required).
<br><i>-t</i>&nbsp;&nbsp;&nbsp;&nbsp;&ensp;&ensp; The number of threads (optional – defaut = 1).
<br><i>-mset</i>&nbsp;&nbsp;A set of amino acid substitution models used to select the best-fit models for alignments (optional, default: FLU,<br>
&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;HIVb, HIVw, JTT, WAG, LG). The models are separated by ‘,’.
<br><br><b><i>Example:</b> python Est.py -f Example</i>

The estimated substitution models will be written to the file <b>[alignment_directory].PAML</b>. The file contains an exchangeability matrix and a frequency vector.

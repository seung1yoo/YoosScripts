# creater : skydj82
# modifier : siyoo

* example working dir : siyoo@titan:/BiO/BioProjects/TBD170324-SNU-Plant-RNAseq-20170529/timeseries
* ln -s ../Rine_Denovo/report/Files/genes.xls
* cat genes.xls | head -n1 | tr '\t' '\n' | grep DEG | grep SELECT > colname.txt

# parameters
## 1. FCFILT (log2)
### log2fc cut 2 ==>  3 ==> 
## 2. p_cut
###
* python mkTimeseriesResult.py genes.xls colname.txt Timeseries

# have to change column index
* python finalTun_timeResult.py

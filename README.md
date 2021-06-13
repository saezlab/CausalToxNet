# Causal-Toxicity-Network (CausalToxNet) Project

The CausalToxNet (CTN) project's repository comprises a snapshot of codes from [CARNIVAL](https://github.com/saezlab/CARNIVAL) with associated input files applied in this project i.e. transcriptional factors weights from [DoRothEA](https://github.com/saezlab/DoRothEA) and pathway activites from [PROGENy](https://github.com/saezlab/progeny). 

In this repository, there are two folders:
- **"CARNIVAL_Raw_Results"** stores the already generated CARNIVAL results which are presented in the manuscript.
- **"CARNIVAL_Scripts"** contains computational scripts to re-generate CARNIVAL results. Within this folder, the pipeline were separated for each of the four datasets from the TG-GATEs and DrugMatrix databases which were analyzed including:
1) **RLR** - rat liver repeated dosing (4,8,15,29 days at low/middle/high doses)
2) **RLS** - rat liver single dosing (3,6,9,24 hours at low/middle/high doses)
3) **RPH** - rat primary hepatocytes (2,8,24 hours at low/middle/high doses)
4) **PHH** - primary human hepatocytes (2,8,24 hours at low/middle/high doses)
5) **DrugMatrix** - rat liver single and repeated dosing (6 hours, 24 hours, 3 days and 7 days at 400 and 1175 mg/kg for carbon tetrachloride) and rat primary hepatocytes (16 hours at 775 mg/kg, 24 hours at 193.7, 387.5 and 775 mg/kg for carbon tetrachloride; 16 hours and 24 hours at 333 mg/kg for monocrotaline)

Within each folder, the customized CARNIVAL scripts for each dataset are stored in the folders started with "CARNIVAL" while the input files for CARNIVAL are stored in the folders "Resources".

In addition, each folder contain a "driver script" where users could re-run the pipeline to generate results as presented in the published article. In the provided script, only a single experimental condition (i.e. a combination of compound, time point and dose  e.g. "acetaminophen, 2 hours, low dose") was assigned in the default setting. Users can change the experimental condition by modifying the counter numbers ('counter_cp' for compound; 'counter_tp' for time point; and 'counter_do' for dose) which could be assigned for parallel computing on computational clusters.

**Important note**: Please change the working directory as well as the file path towards cplex interactive solver in the driver script before running throught the pipeline respectively. For more information on CARNIVAL, please refer to the CARNIVAL GitHub page: [https://github.com/saezlab/CARNIVAL](https://github.com/saezlab/CARNIVAL).

## Author

Panuwat Trairatphisan (panuwat.trairatphisan -at- gmail.com)

## License

Distributed under the GNU GPLv3 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/CARNIVAL/blob/master/LICENSE.txt) or copy at [http://www.gnu.org/licenses/gpl-3.0.html](http://www.gnu.org/licenses/gpl-3.0.html).

## References

[CARNIVAL](https://www.nature.com/articles/s41540-019-0118-z):

> Liu A, Trairatphisan P, Gjerga E, Didangelos A, Barrett J, Saez-Rodriguez J. (2019). From expression footprints to causal pathways: contextualizing large signaling networks with CARNIVAL. *NPJ Systems Biology and Applications*, Issue 5, Nr. 40.

[DoRothEA v2 - Garcia-Alonso et al.](https://www.biorxiv.org/content/early/2018/06/03/337915):

> Garcia-Alonso L, Ibrahim MM, Turei D, Saez-Rodriguez J. (2018). Benchmark and integration of resources for the estimation of human transcription factor activities. *bioRXiv*, https://doi.org/10.1101/337915.

[PROGENy - Schubert et al.](https://www.nature.com/articles/s41467-017-02391-6):

> Schubert M, Klinger B, Klünemann M, Sieber A, Uhlitz F, Sauer S, Garnett MJ, Blüthgen N, Saez-Rodriguez J. (2018). Perturbation-response genes reveal signaling footprints in cancer gene expression. *Nature Communication*, Issue 9, Nr. 20. https://doi.org/10.1038/s41467-017-02391-6.


## Acknowledgement

CARNIVAL has been developed as a computational tool to analyse -omics data within the [TransQST Consortium](https://transqst.org)

"This project has received funding from the Innovative Medicines Initiative 2 Joint Undertaking under grant agreement No 116030. The Joint Undertaking receives support from the European Union's Horizon 2020 research and innovation programme and EFPIA."

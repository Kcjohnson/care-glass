## Longitudinal DNA sequencing analysis in CARE/GLASS datasets

### Overview
The Cellular Analysis of Resistance and Evolution (CARE) and Glioma Longitudinal AnalySiS (GLASS) consortiums both use bulk sequencing to study glioma evolution. The code in this respository was used to perform the pre-processing of whole exome and whole genome sequencing analyses described in both the CARE and GLASS consortium datasets (i.e., same pipelines were used for both datasets). The larger CARE/GLASS dataset is managed internally by the Verhaak Lab using the PostgreSQL database management system. These data are under active curation so future versions will include additional data as well as correct potential errors.

### Data download
The CARE IDH-wildtype data can be downloaded from the `Tables` page [here](https://www.synapse.org/Synapse:syn61979590/tables/) and the `Files` page [here](https://www.synapse.org/Synapse:syn61979590/files/). It is also possible to query the data directly using the the API by using queries. You can read more about that [here](https://docs.synapse.org/articles/tables.html).

### Prior datasets

Prior publications and repositories associated with these workflows are listed below:

Varn, F.S., Johnson, K.C., et al. (2022). Glioma progression is shaped by genetic evolution and microenvironment interactions. Cell (2022). [Paper.](https://pubmed.ncbi.nlm.nih.gov/35649412/) [Code.](https://github.com/fsvarn/GLASSx)

Barthel, F.P., Johnson, K.C., Varn, F.S., Moskalik, A.D., Tanner, G., Kocakavuk, E., Anderson, K.J., Abiola, O., Aldape, K., Alfaro, K.D., et al. (2019). Longitudinal molecular trajectories of diffuse glioma in adults. Nature 576, 112-120. [Paper.](https://www.nature.com/articles/s41586-019-1775-1) [Code.](https://github.com/TheJacksonLaboratory/GLASS)

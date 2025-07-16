# health-rating-index
## Environmental Health Burden in Dublin

This project estimates the population-level health burden from environmental exposures at the small area in Dublin, Ireland — specifically road traffic noise and air pollution (PM2.5, NO₂, and O<sub>3</sub>). Using publicly available spatial datasets and risk models from the World Health Organization ([WHO, 2022](https://www.eionet.europa.eu/etcs/etc-he/products/etc-he-products/etc-he-reports/etc-he-report-2022-10-health-risk-assessment-of-air-pollution-and-the-impact-of-the-new-who-guidelines/@@download/file/ETC%20HE%202022-10_Eionet_report_HRA_FINAL_28-11-2022.pdf)) and existing academic work ([Hanigan et al. 2019](https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-019-0184-x)), this analysis calculates the fraction of mortality burden attributable to noise and air pollution across different health outcomes. It also includes pre-calculated values for access to health 'benefits', including GP and green space accessibility, which are combined with the health burdens callated to create a holistic 'Health Rating Index' (HRI). The HRI is then used as a dependent variable in a 'spatial' random forest model to understand its relationship to social deprivation, non-auto commuting, population density, distance to the nearest primary and secondary road, and various neighbourhood dummies. For more detail on the methodology, please reference the conference paper below.

The core outputs include:
- Population-weighted exposure assessments at small area level
- Attributable fractions or odds ratios by cause
- Transparent, reproducible code using R

This tool is intended to support public health and urban policy decision-making by quantifying the health impacts of environmental exposures in a local context.

🗺️ **Location:** Dublin, Ireland  
🧪 **Tools:** R, tidyverse, sf, data.table  
📊 **Data Sources:** CSO, OSM, [greenR](https://github.com/sachit27/greenR), [Ireland's Open Data Portal](https://data.gov.ie/dataset/family-practice-gp-sites), [Google AirView](https://data.gov.ie/dataset/google-airview-data-dublin-city), [Dublinked](https://data.smartdublin.ie/dataset/noise-maps-from-traffic-sources-in-dublin-city-council), [Pobal](https://data.gov.ie/dataset/pobal-hp-deprivation-index-scores-2022), WHO environmental burden of disease models

Please cite the conference paper when using the prepared data inputs or code:

Credit, K., Damanpreet, K., and Eccles, E. 2025. “Exploring the transport-health-environment nexus through a new 'Health Rating Index' for Dublin, Ireland.” _Proceedings of the 33rd GISRUK Conference_. DOI: https://doi.org/10.5281/zenodo.15183740. 

This work was funded by the Irish Research Council (IRC) under the 2022 New Foundations scheme Strand 1a as a part of the project [Dublin 8 Health + Environment Data Dashboard](https://experience.arcgis.com/experience/04749d06fd0e43d9a58d2e644a4bc71f/). In addition to the authors, important contributions to the project were made by the staff of the Robert Emmet Community Development Project (RECDP), the South Inner City Community Development Association (SICCDA), Prof. Mary Corcoran, Dr. Lidia Manzo, Shayal Kumar, and the various Dublin 8 community members who participated in the project. 

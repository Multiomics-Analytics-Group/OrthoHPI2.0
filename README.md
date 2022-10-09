# OrthoHPI2.0

OrthoHPI 2.0 is a resource that provides homology-derived predictions of host-parasite protein-protein interactions (PPI). This work renews and extends OrthoHPI by integrating new versions of databases, predictors and proteomes of 24 eukaryotic parasites such as apicomplexa and trypanosomatids.

This resource provides tissue and cell type resolution interactions as well as annotation of Biological Processes that allows a comprehensive analysis and comparison of all the parasitic species studied.


<img width="1760" alt="Screenshot 2022-10-08 at 20 25 11" src="https://user-images.githubusercontent.com/1425851/194722121-b6f01c52-57d7-4676-aefe-239a2be3e78a.png">


## Development

OrthoHPI 2.0 has been developed entirely with **Python 3.8.2** and Streamlit (https://docs.streamlit.io/) for visualization and analysis of the results.

The **homology prediction** has been done using:

<a src="http://eggnog5.embl.de/"><img width=100 alt="eggnog" src="https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0/blob/main/images/eggnog.png"> </a>


<a src="https://string-db.org/"><img width=100 alt="string" src="https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0/blob/main/images/string.png"></a>


The **tissue and cellular compartment annotation** used to filter the proteomes come from:

<a src="https://compartments.jensenlab.org/"><img width=100 alt="compartments" src="https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0/blob/main/images/compartments.png"></a>

<a src="https://tissues.jensenlab.org/"><img width=100 alt="tissues" src="https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0/blob/main/images/tissues.png"></a>


The **cell type expression** information has been downloaded from:

<a src="https://www.proteinatlas.org/humanproteome/tissue+cell+type"><img width=200 alt="hpa" src="https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0/blob/main/images/hpa.png"></a>

## Getting started

### Running the website

All the necessary data (predictions and annotations) have been precalculated so you can run the web server:
``` 
$ streamlit run orthohpi_web.py
```

This will open your browser at http://localhost:8501 with the functionality shown in the video:



https://user-images.githubusercontent.com/1425851/194722960-e42a191f-1cec-4c49-a96e-5a03b351677d.mp4


### Installation

If you want to rerun the predictions, you can install and run the pipeline following these instructions.

```
$ pip install -r requirements.txt
```

To obtain the Human-Parasite PPI predictions run:
``` 
$ python main.py
```



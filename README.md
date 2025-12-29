# Unraveling the role of mitochondria dysfunction in early sepsis

![R](https://img.shields.io/badge/R-4.3-276DC3?logo=r)
![Python](https://img.shields.io/badge/Python-3.11-blue?logo=python)
![DESeq2](https://img.shields.io/badge/DESeq2-1.42-green)
![scikit--learn](https://img.shields.io/badge/sklearn-1.3-F7931E?logo=scikit-learn)
![License](https://img.shields.io/badge/License-MIT-green)

This research was conducted by Jamie Scheper (373689)

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

---

# Introduction
Sepsis, a major contributor to in-hospital deaths, is challenging to detect due to its heterogeneous nature, and dysfunctional mitochondria can worsen the condition. Analyzing mitochondria-related genes could aid in early detection and guide treatment.
RNA-Seq and clinical data from 348 septic patients, collected from four emergency rooms (ER) and one intensive care unit (ICU), were used and compared to 44 healthy controls. This study applied supervised and unsupervised algorithms to identify clinically relevant gene signatures.

Mitochondria-related differentially expressed genes (DEGs) were identified by comparing levels of severity based on Sequential Organ Failure Assessment (SOFA) scores. These DEGs helped to establish two unique severity-based endotypes in ER cohorts, which were distinct from the healthy controls and validated on an ICU cohort. This categorization addresses sepsis heterogeneity. Feature selection identified three gene sets that could accurately predict endotype groups. Notably, a logistic regression algorithm with L1 regularization predicted endotypes with 94% and 90% accuracy in ER and ICU cohorts using an eleven-gene set.

The gene signatures and endotypes indicated the essential role of mitochondrial genes in early sepsis detection and as potential novel biomarkers. Future research should aim to develop a multi-modal prediction tool for enhanced patient stratification.

In the flowchart depicted below, we highlight the methodology in more detail. Parts of the research all have their distinct color. Firstly, we focused on establishing the data distribution, imputating NAs, removing redundancies, outliers, and anomolies, and whether normalization was necessary in the form of an exploratory data analysis (EDA) (depicted in the color green). After that, we focused on discovering differentially expressed genes (DEGs) on the severity and mortality variables with DESeq2 (depicted in yellow). The main focus was on establishing an appropriate fold-change threshold. We validated our findings by comparing them with results from edgeR. After that, in purple, we used all the DEGs to form severity-based clusters/endotypes by comparing three different methods: K-means, K-medoids (PAM), and hierarchical clustering. We compared various distance and, where appropriate, linkage methods. Utilizing this, we discovered the existence of two endotypes in the ER cohorts, which we validated with various gene set sizes derived from MAD and the ICU cohort. We developed a machine learning model based on eleven DEGs extracted by comparing multiple feature selection techniques. The machine learning model and these DEGs could accurately predict the endotype and severity status using ER and ICU cohorts (depicted in red).

![Flowchart](misc/figures_report/flowchart.jpg)

# Setup and Usage
We used [Rstudio][r] (version 4.3.1) to clean, restructure, and visualize data and findings. Additionally, we used R primarily for DE and pathway analysis and unsupervised clustering. For that, we used a couple of external packages/libraries:

- [dplyr (version 1.1.3)][dplyr]: data manipulation
- [tidyverse (version 2.0.0)][tidyverse]: data manipulation
- [caret (version 6.0-94)][caret]: removal of near-zero attributes
- [ggplot2 (version 3.4.3)][ggplot2]: data visualization
- [DESeq2 (version 1.42.0)][DESeq2]: DE analysis
- [edgeR (version 4.0.2)][edgeR]: DE analysis
- [ReactomePA (version 1.46.0)][ReactomePA]: pathway analysis
- [enrichR (version 3.2)][enrichR]: pathway analysis
- [SVA (version 3.50.0)][SVA]: normalization (batch correction)
- [dendextend (version 1.17.1)][dendextend]: hierarchical clustering visualization
- [ConsensusClusterPlus (version 1.66.0)][ConsensusClusterPlus]: metric for optimal number of clusters
- [factoextra (version 2.9)][factoextra]: metrics stability clusters and data visualization
- [NbClust (version 3.0.1)][NbClust]: metric for optimal number of clusters
- [ComplexHeatmap (version 2.18.0)][ComplexHeatmap]: data visualization

To install these R packages, please follow the following steps:


Additionally, we performed feature selection and machine learning in [Python][python] (version 3.11.2). An [Jupyter][jupyter] (version 1.0.0) notebook is available under `logs/` wherein we perform all these processes. We used the following packages:

- [Pandas (version 2.1.3)][Pandas]: data manipulation
- [NumPy (version 1.26.2)][NumPy]: mathematics and data manipulation
- [sklearn (version 1.3.2)][sklearn]: machine learning techniques
- [seaborn (version 0.13.0)][seaborn]: data visualization
- [matplotlib (version 3.8.2)][matplotlib]: data visualization
- [scikit-plot (version 0.3.7)][scikit_plot]: data visualization and machine learning
- [mxlextend (version 0.23.0)][mxlextend]: feature selection

To install these Python packages, please follow the following steps:

1) Open Rstudio

Install Rstudio here:
- [Rstudio][r]

2) Install packages

Packages can be install via two ways:

The simple way: install them one by one via the command line:
```shell
install.packages("tidyverse")
```

Via a package manager:

Some packages are only available through `BiocManager` such as `DESeq2`.

Always install `BiocManager` via the regular route.

```shell
install.packages("BiocManager")
```

After that, install your package:

```shell
BiocManager::install("DESeq2")
```

1) Setting up a virtual environment.

Go to your project directory.

```shell
cd path\to\your\project\directory
```

Then, create the virtual environment:

```shell
python -m venv venv
```

Activate the virtual environment:

For Windows:

```shell
.\venv\Scripts\activate
```

For Linux:

```shell
source env_name/bin/activate
```

2) Installation via requirements.txt

Install all the libraries:

```shell
pip install -r misc/requirements.txt
```

3) Loading in the Juypter notebook

Ensure that Jupyter is installed.

```shell
pip install jupyter
```

You can either open the `.ipynb` with your preffered code editor such as Visual Studio Code or via your local web browser.

For the web browser option, execute:

```shell
jupyter notebook
```

# Layout
In the root of this repository, you will find the end report which sums up the findings, the problems we came across and how we challenged those. The subdirectories are set up as follows:

- `data/`: publicly available data, including raw counts and some clinical data points (see `SRATableRun.txt`). Please keep in mind that most clinical data is missing, as this is considered sensitive data.
    - `source/`: contains all the source files and are primarily used in the `EDA.Rmd`.
    - `degs/`: contains the DEG dataset found by DESeq2 and is used in all logs.
    - `counts/`: various subsets of count data, patient location specific (e.g., healthy control counts).
    - `clustering/`: clinical data critical to the endotype classification.

- `logs/`: logs regarding each part of the research. Herein, we explain why certain choices were made. Logs are aimed to be reproducible.
- `pdfs/`: PDFs of logs
- `scripts/`: scripts with helper functions. Are primiarly used in the logs. `helper_functions.R` is for all R-based logs and `train_and_evaluate.py` for Python.
- `misc/`: miscellaneous essential figures, data management plan, poster etc.

In the report, we refer to logs as supplemental material. Here, we mean:

- Supplemental 1: `EDA.Rmd`
- Supplemental 2: `de_and_pathway_analysis.Rmd`
- Supplemental 3: `clustering.Rmd`
- Supplemental 4: `Feature_Selection.ipynb`
- Supplemental 5: `differences_genes.Rmd`

# Contact
If any issue or question remains, please contact us at [j.s.c.scheper@st.hanze.nl](mailto:j.s.c.scheper@st.hanze.nl) or [j.s.c.scheper@umcg.nl](mailto:j.s.c.scheper@umcg.nl).

[R]: https://www.rstudio.com/products/rstudio/download/
[dplyr]: https://dplyr.tidyverse.org/
[tidyverse]: https://www.tidyverse.org/
[caret]: https://cran.r-project.org/web/packages/caret/index.html
[ggplot2]: https://ggplot2.tidyverse.org/
[DESeq2]: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
[edgeR]: https://bioconductor.org/packages/release/bioc/html/edgeR.html
[ReactomePA]: https://bioconductor.org/packages/release/bioc/html/ReactomePA.html
[enrichR]: https://maayanlab.cloud/Enrichr/
[SVA]: https://bioconductor.org/packages/release/bioc/html/sva.html
[dendextend]: https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html
[ConsensusClusterPlus]: https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html
[factoextra]: https://cran.r-project.org/web/packages/factoextra/index.html
[NbClust]: https://cran.r-project.org/web/packages/NbClust/index.html
[ComplexHeatmap]: https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
[python]: https://www.python.org/
[jupyter]: https://jupyter.org/
[Pandas]: https://pandas.pydata.org/
[NumPy]: https://numpy.org/
[sklearn]: https://scikit-learn.org/stable/
[seaborn]: https://seaborn.pydata.org/
[matplotlib]: https://matplotlib.org/
[scikit_plot]: https://pypi.org/project/scikit-plot/
[mxlextend]: https://rasbt.github.io/mlxtend/
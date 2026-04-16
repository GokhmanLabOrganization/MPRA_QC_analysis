FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="df14405700f8890b1518f033680d309f092796cd306fbc3e966730f44838ccc5"
# Conda environment:
#   source: workflow/envs/default.yml
#   prefix: /conda-envs/d2f227aa6d7bbb6926a98fb7d83aa2bc
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - biopython
#     - matplotlib
#     - numpy
#     - pandas
#     - pysam
#     - regex
#     - scikit-learn
#     - scipy
#     - seaborn
#     - click
RUN mkdir -p /conda-envs/d2f227aa6d7bbb6926a98fb7d83aa2bc
COPY workflow/envs/default.yml /conda-envs/d2f227aa6d7bbb6926a98fb7d83aa2bc/environment.yaml

# Conda environment:
#   source: workflow/envs/bcalm.yml
#   prefix: /conda-envs/86d187db96ee09d7a21b9b2009185165
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - r-base
#     - r-bcalm
#     - r-argparse
#     - r-dplyr
#     - r-ggplot2
#     - r-tidyr
#     - r-tibble
RUN mkdir -p /conda-envs/86d187db96ee09d7a21b9b2009185165
COPY workflow/envs/bcalm.yml /conda-envs/86d187db96ee09d7a21b9b2009185165/environment.yaml

# Conda environment:
#   source: workflow/envs/mpralib.yml
#   prefix: /conda-envs/da3ac9efbe61038ecd2a73a053d50729
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - mpralib==0.10.*

RUN mkdir -p /conda-envs/da3ac9efbe61038ecd2a73a053d50729
COPY workflow/envs/mpralib.yml /conda-envs/da3ac9efbe61038ecd2a73a053d50729/environment.yaml


RUN conda env create --prefix /conda-envs/d2f227aa6d7bbb6926a98fb7d83aa2bc --file /conda-envs/d2f227aa6d7bbb6926a98fb7d83aa2bc/environment.yaml && \
conda env create --prefix /conda-envs/86d187db96ee09d7a21b9b2009185165 --file /conda-envs/86d187db96ee09d7a21b9b2009185165/environment.yaml && \
conda env create --prefix /conda-envs/da3ac9efbe61038ecd2a73a053d50729 --file /conda-envs/da3ac9efbe61038ecd2a73a053d50729/environment.yaml && \
conda clean --all -y

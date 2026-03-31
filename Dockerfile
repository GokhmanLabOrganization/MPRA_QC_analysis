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

RUN conda env create --prefix /conda-envs/d2f227aa6d7bbb6926a98fb7d83aa2bc --file /conda-envs/d2f227aa6d7bbb6926a98fb7d83aa2bc/environment.yaml && \
    conda clean --all -y

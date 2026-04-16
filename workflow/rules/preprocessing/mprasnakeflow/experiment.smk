rule preprocessing_mprasnakeflow_experiment_bcalm_elements:
    input:
        reporter_experiment_barcode=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("reporter_experiment_barcode", []),
        labels=config.get("mprasnakeflow", {}).get("experiment", {}).get("labels", []),
        script=getScript("preprocessing/mprasnakeflow/bcalm_elements.R"),
    output:
        output_volcano_plot=report(
            "results/{project}/preprocessing/mprasnakeflow/activity/bcalm/element_volcano_plot.png",
            caption=getReport("preprocessing/bcalm_elements.rst"),
            category="{project}",
            subcategory="Preprocessing",
            labels={
                "analysis": "Preprocessing",
                "type": "Activity",
                "figure": "Volcano",
            },
        ),
        output_density_plot=report(
            "results/{project}/preprocessing/mprasnakeflow/activity/bcalm/element_density_plot.png",
            caption=getReport("preprocessing/bcalm_elements.rst"),
            category="{project}",
            subcategory="Preprocessing",
            labels={
                "analysis": "Preprocessing",
                "type": "Activity",
                "figure": "Density",
            },
        ),
        output_mean_variance_relation=report(
            "results/{project}/preprocessing/mprasnakeflow/activity/bcalm/element_mean_variance_relation.png",
            caption=getReport("preprocessing/bcalm_elements.rst"),
            category="{project}",
            subcategory="Preprocessing",
            labels={
                "analysis": "Preprocessing",
                "type": "Activity",
                "figure": "Mean-Variance Relation",
            },
        ),
        bcalm_result="results/{project}/preprocessing/mprasnakeflow/activity/bcalm/elements.tsv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/bcalm_elements.{project}.log",
    conda:
        getCondaEnv("bcalm.yml")
    container:
        "docker://visze/bcalm:0.100.0"
    params:
        test_label=config.get("mprasnakeflow", {}).get("experiment", {}).get("test_label", ""),
        control_label=config.get("mprasnakeflow", {}).get("experiment", {}).get("control_label", ""),
        percentile=1.0 - (config.get("mprasnakeflow", {}).get("experiment", {}).get("fdr", 0.1) / 2.0),
        normalize=("TRUE" if config.get("mprasnakeflow", {}).get("experiment", {}).get("normalize", True) else "FALSE"),
    shell:
        """
        Rscript {input.script} --counts {input.reporter_experiment_barcode} \
        --labels {input.labels} --test-label {params.test_label} --control-label {params.control_label} \
        --percentile {params.percentile} --normalize {params.normalize} \
        --output-mean-variance-relation-plot {output.output_mean_variance_relation} --output-volcano-plot {output.output_volcano_plot} --output-density-plot {output.output_density_plot} \
        --output {output.bcalm_result} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_activity_df:
    input:
        reporter_experiment_barcode=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("reporter_experiment_barcode", []),
        bcalm_result="results/{project}/preprocessing/mprasnakeflow/activity/bcalm/elements.tsv.gz",
        script=getScript("preprocessing/mprasnakeflow/create_activity_df.py"),
    output:
        activity_df="results/{project}/preprocessing/mprasnakeflow/activity/activity_df.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/activity_df.{project}.log",
    conda:
        getCondaEnv("mpralib.yml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.3--pyhdfd78af_0"
    params:
        fdr=config.get("mprasnakeflow", {}).get("experiment", {}).get("fdr", 0.1),
    shell:
        """
        python {input.script} \
        --reporter-experiment-barcode {input.reporter_experiment_barcode} --bcalm-statistics {input.bcalm_result} \
        --output {output.activity_df} --fdr {params.fdr} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_activity_per_rep:
    input:
        reporter_experiment_barcode=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("reporter_experiment_barcode", []),
        script=getScript("preprocessing/mprasnakeflow/create_activity_per_rep.py"),
    output:
        activity_per_rep_df="results/{project}/preprocessing/mprasnakeflow/activity/activity_per_rep.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/activity_per_rep.{project}.log",
    conda:
        getCondaEnv("mpralib.yml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.3--pyhdfd78af_0"
    shell:
        """
        python {input.script} \
        --reporter-experiment-barcode {input.reporter_experiment_barcode} \
        --output {output.activity_per_rep_df} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_bcalm_comparative:
    input:
        reporter_experiment_barcode=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("reporter_experiment_barcode", []),
        comparative_map=config.get("mprasnakeflow", {}).get("experiment", {}).get("comparative_map", []),
        script=getScript("preprocessing/mprasnakeflow/bcalm_comparative.R"),
    output:
        bcalm_result="results/{project}/preprocessing/mprasnakeflow/activity/bcalm/comparative.tsv.gz",
        output_volcano_plot=report(
            "results/{project}/preprocessing/mprasnakeflow/activity/bcalm/comparative_volcano_plot.png",
            caption=getReport("preprocessing/bcalm_comparative.rst"),
            category="{project}",
            subcategory="Preprocessing",
            labels={
                "analysis": "Preprocessing",
                "type": "Comparative Activity",
                "figure": "Volcano",
            },
        ),
        output_mean_variance_relation=report(
            "results/{project}/preprocessing/mprasnakeflow/activity/bcalm/comparative_mean_variance_relation.png",
            caption=getReport("preprocessing/bcalm_comparative.rst"),
            category="{project}",
            subcategory="Preprocessing",
            labels={
                "analysis": "Preprocessing",
                "type": "Comparative Activity",
                "figure": "Mean-Variance Relation",
            },
        ),
    log:
        "logs/preprocessing/mprasnakeflow/experiment/bcalm_comparative.{project}.log",
    conda:
        getCondaEnv("bcalm.yml")
    container:
        "docker://visze/bcalm:0.100.0"
    params:
        normalize=("TRUE" if config.get("mprasnakeflow", {}).get("experiment", {}).get("normalize", True) else "FALSE"),
    shell:
        """
        Rscript {input.script} --counts {input.reporter_experiment_barcode} \
        --comparative-map {input.comparative_map} --normalize {params.normalize} \
        --output-mean-variance-relation-plot {output.output_mean_variance_relation} --output-volcano-plot {output.output_volcano_plot} \
        --output {output.bcalm_result} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_comparative_df:
    input:
        comparative_map=config.get("mprasnakeflow", {}).get("experiment", {}).get("comparative_map", []),
        bcalm_result="results/{project}/preprocessing/mprasnakeflow/activity/bcalm/comparative.tsv.gz",
        script=getScript("preprocessing/mprasnakeflow/create_comparative_df.py"),
    output:
        comparative_df="results/{project}/preprocessing/mprasnakeflow/activity/comparative_df.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/comparative_df.{project}.log",
    conda:
        getCondaEnv("mpralib.yml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.3--pyhdfd78af_0"
    params:
        fdr=config.get("mprasnakeflow", {}).get("experiment", {}).get("fdr", 0.1),
    shell:
        """
        python {input.script} \
        --bcalm-statistics {input.bcalm_result} --comparative-map {input.comparative_map} \
        --output {output.comparative_df} --fdr {params.fdr} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_allelic_pairs_df:
    input:
        reporter_experiment_barcode=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("reporter_experiment_barcode", []),
        comparative_map=config.get("mprasnakeflow", {}).get("experiment", {}).get("comparative_map", []),
        script=getScript("preprocessing/mprasnakeflow/create_allelic_pairs_df.py"),
    output:
        allelic_pairs_df="results/{project}/preprocessing/mprasnakeflow/activity/allelic_pairs_df.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/allelic_pairs_df.{project}.log",
    conda:
        getCondaEnv("mpralib.yml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.3--pyhdfd78af_0"
    shell:
        """
        python {input.script} \
        --reporter-experiment-barcode {input.reporter_experiment_barcode} --comparative-map {input.comparative_map} \
        --output {output.allelic_pairs_df} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_allelic_pairs_replicates_df:
    input:
        reporter_experiment_barcode=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("reporter_experiment_barcode", []),
        comparative_map=config.get("mprasnakeflow", {}).get("experiment", {}).get("comparative_map", []),
        script=getScript("preprocessing/mprasnakeflow/create_allelic_pairs_replicates_df.py"),
    output:
        allelic_pairs_replicates_df="results/{project}/preprocessing/mprasnakeflow/activity/allelic_pairs_replicates_df.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/allelic_pairs_replicates_df.{project}.log",
    conda:
        getCondaEnv("mpralib.yml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.3--pyhdfd78af_0"
    shell:
        """
        python {input.script} \
        --reporter-experiment-barcode {input.reporter_experiment_barcode} --comparative-map {input.comparative_map} \
        --output {output.allelic_pairs_replicates_df} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_downsampling:
    input:
        reporter_experiment_barcode=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("reporter_experiment_barcode", []),
        script=getScript("preprocessing/mprasnakeflow/activity_downsample.py"),
    output:
        temp(
            expand(
                "results/{{project}}/preprocessing/mprasnakeflow/activity/downsample_tmp/reporter_barcode_{fraction:.1f}.tsv.gz",
                fraction=[0.1 * i for i in range(1, 10)],
            )
        ),
    log:
        "logs/preprocessing/mprasnakeflow/experiment/downsampling_ratio_df.{project}.log",
    conda:
        getCondaEnv("mpralib.yml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.3--pyhdfd78af_0"
    shell:
        """
        paths=( {output} );
        python {input.script} --reporter-experiment-barcode {input.reporter_experiment_barcode} --output-folder $(dirname "${{paths[0]}}") > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_downsampling_ratio_df:
    input:
        reporter_experiment_barcode=lambda wc: (
            f"results/{wc.project}/preprocessing/mprasnakeflow/activity/downsample_tmp/reporter_barcode_{wc.fraction}.tsv.gz"
            if wc.fraction != "1.0"
            else config.get("mprasnakeflow", {}).get("experiment", {}).get("reporter_experiment_barcode", [])
        ),
        script=getScript("preprocessing/mprasnakeflow/create_downsampling_ratio_df.py"),
    output:
        "results/{project}/preprocessing/mprasnakeflow/activity/downsample_tmp/ratio_df_{fraction}.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/downsampling_ratio_df.{project}.{fraction}.log",
    conda:
        getCondaEnv("mpralib.yml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.3--pyhdfd78af_0"
    shell:
        """
        python {input.script} --reporter-experiment-barcode {input.reporter_experiment_barcode} --output {output} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_downsampling_ratio_copy:
    """
Copy to the final downsample folder. This is necessary to avoid issues with the downstream rules that expect the files to be in a specific location. The downsample folder is also added to the activity_files dataframe for the downstream rules.
"""
    input:
        final=expand(
            "results/{{project}}/preprocessing/mprasnakeflow/activity/downsample_tmp/ratio_df_{fraction:.1f}.csv.gz",
            fraction=[0.1 * i for i in range(1, 11)],
        ),
    output:
        output_path=directory("results/{project}/preprocessing/mprasnakeflow/activity/downsampling_ratio/"),
    log:
        "logs/preprocessing/mprasnakeflow/experiment/downsampling_ratio_copy.{project}.log",
    conda:
        getCondaEnv("default.yml")
    shell:
        """
        mkdir -p {output.output_path}
        for i in {input.final}; do
            cp $i {output.output_path}/$(basename $i);
        done 2> {log}
        """


rule preprocessing_mprasnakeflow_experiment_downsampling_bcalm:
    input:
        reporter_experiment_barcode="results/{project}/preprocessing/mprasnakeflow/activity/downsample_tmp/reporter_barcode_{fraction}.tsv.gz",
        labels=config.get("mprasnakeflow", {}).get("experiment", {}).get("labels", []),
        script=getScript("preprocessing/mprasnakeflow/bcalm_elements.R"),
    output:
        bcalm_result="results/{project}/preprocessing/mprasnakeflow/activity/downsample_tmp/elements_{fraction}.tsv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/downsampling_bcalm.{project}.{fraction}.log",
    conda:
        getCondaEnv("bcalm.yml")
    container:
        "docker://visze/bcalm:0.100.0"
    params:
        test_label=config.get("mprasnakeflow", {}).get("experiment", {}).get("test_label", ""),
        control_label=config.get("mprasnakeflow", {}).get("experiment", {}).get("control_label", ""),
        percentile=1.0 - (config.get("mprasnakeflow", {}).get("experiment", {}).get("fdr", 0.1) / 2.0),
        normalize=("TRUE" if config.get("mprasnakeflow", {}).get("experiment", {}).get("normalize", True) else "FALSE"),
    shell:
        """
        Rscript {input.script} --counts {input.reporter_experiment_barcode} \
        --labels {input.labels} --test-label {params.test_label} --control-label {params.control_label} \
        --percentile {params.percentile} --normalize {params.normalize} \
        --output {output.bcalm_result} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_downsample_activity_df:
    input:
        reporter_experiment_barcode=lambda wc: (
            f"results/{wc.project}/preprocessing/mprasnakeflow/activity/downsample_tmp/reporter_barcode_{wc.fraction}.tsv.gz"
            if wc.fraction != "1.0"
            else config.get("mprasnakeflow", {}).get("experiment", {}).get("reporter_experiment_barcode", [])
        ),
        bcalm_result=lambda wc: (
            f"results/{wc.project}/preprocessing/mprasnakeflow/activity/downsample_tmp/elements_{wc.fraction}.tsv.gz"
            if wc.fraction != "1.0"
            else f"results/{wc.project}/preprocessing/mprasnakeflow/activity/bcalm/elements.tsv.gz"
        ),
        script=getScript("preprocessing/mprasnakeflow/create_downsampling_activity_df.py"),
    output:
        activity_df="results/{project}/preprocessing/mprasnakeflow/activity/downsample_tmp/activity_df_{fraction}.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/downsample_activity_df.{project}.{fraction}.log",
    conda:
        getCondaEnv("mpralib.yml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.3--pyhdfd78af_0"
    params:
        fdr=config.get("mprasnakeflow", {}).get("experiment", {}).get("fdr", 0.1),
    shell:
        """
        python {input.script} \
        --reporter-experiment-barcode {input.reporter_experiment_barcode} --bcalm-statistics {input.bcalm_result} \
        --output {output.activity_df} --fdr {params.fdr} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_downsampling_activity_copy:
    """
Copy to the final downsample folder. This is necessary to avoid issues with the downstream rules that expect the files to be in a specific location. The downsample folder is also added to the activity_files dataframe for the downstream rules.
"""
    input:
        final=expand(
            "results/{{project}}/preprocessing/mprasnakeflow/activity/downsample_tmp/activity_df_{fraction:.1f}.csv.gz",
            fraction=[0.1 * i for i in range(1, 11)],
        ),
    output:
        output_path=directory("results/{project}/preprocessing/mprasnakeflow/activity/downsampling_activity/"),
    log:
        "logs/preprocessing/mprasnakeflow/experiment/downsampling_activity_copy.{project}.log",
    conda:
        getCondaEnv("default.yml")
    shell:
        """
        mkdir -p {output.output_path}
        for i in {input.final}; do
            cp $i {output.output_path}/$(basename $i);
        done 2> {log}
        """

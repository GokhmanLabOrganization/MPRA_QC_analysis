rule activity_main:
    input:
        activity=lookup(
            within=activity_files,
            query="file == 'activity_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_main_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/Activity_distribution.png",
            caption=getReport("activity/Activity_distribution.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Distribution",
                "figure": "Activity Distribution",
            },
        ),
        report(
            "results/{project}/activity/P_value_distribution.png",
            caption=getReport("activity/P_value_distribution.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Distribution",
                "figure": "P-Value Distribution",
            },
        ),
        report(
            "results/{project}/activity/Cumulative_RNA_reads.png",
            caption=getReport("activity/Cumulative_RNA_reads.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Distribution",
                "figure": "Cumulative RNA Reads",
            },
        ),
        report(
            "results/{project}/activity/RNA_vs_DNA_w_bar.png",
            caption=getReport("activity/RNA_vs_DNA_w_bar.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Sanity check",
                "figure": "RNA vs DNA",
            },
        ),
        report(
            "results/{project}/activity/Activity_statistic_vs_count_ratio.png",
            caption=getReport("activity/Activity_statistic_vs_count_ratio.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Sanity check",
                "figure": "Activity Statistic vs Count Ratio",
            },
        ),
    log:
        "logs/activity/main.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
        python {input.script} activity main \
        --activity {input.activity} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_control_boxplots:
    input:
        activity=lookup(
            within=activity_files,
            query="file == 'activity_df'",
            cols=["path"],
        ),
        controls=lookup(
            within=activity_files,
            query="file == 'control_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_control_boxplots_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/Activity_of_controls.png",
            caption=getReport("activity/Activity_of_controls.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Controls",
                "figure": "Activity of Controls",
            },
        ),
    log:
        "logs/activity/control_boxplots.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity control-boxplots --activity {input.activity} --controls {input.controls} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_replicability_by_activity:
    input:
        activity=lookup(
            within=activity_files,
            query="file == 'activity_df'",
            cols=["path"],
        ),
        activity_per_rep=lookup(
            within=activity_files,
            query="file == 'activity_per_rep'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_replicability_by_activity_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/Replicability_by_activity.png",
            caption=getReport("activity/Replicability_by_activity.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Reproducibility",
                "figure": "Replicability by Activity",
            },
        ),
    log:
        "logs/activity/replicability_by_activity.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity replicability-by-activity --activity {input.activity} --activity-per-rep {input.activity_per_rep} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_gc_content_bias:
    input:
        design=lookup(
            within=activity_files,
            query="file == 'cCRE_fasta'",
            cols=["path"],
        ),
        activity=lookup(
            within=activity_files,
            query="file == 'activity_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_gc_content_bias_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/DNA_counts_vs_GC_content.png",
            caption=getReport("activity/DNA_counts_vs_GC_content.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Bias",
                "figure": "GC Content",
            },
        ),
    log:
        "logs/activity/gc_content_bias.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity gc-content-bias --activity {input.activity} --design {input.design} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_ratio_correlation_between_replicates:
    input:
        activity_per_rep=lookup(
            within=activity_files,
            query="file == 'activity_per_rep'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_ratio_correlation_between_replicates_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/Retained_cCREs_and_BCs.png",
            caption=getReport("activity/Retained_cCREs_and_BCs.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Retention",
                "figure": "Retained cCREs and BCs",
            },
        ),
        report(
            [
                "results/{project}/activity/Correlation_between_replicates.png",
                "results/{project}/activity/Correlation_between_replicates_w_bar.png",
            ],
            caption=getReport("activity/Correlation_between_replicates.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Reproducibility",
                "figure": "Correlation between Replicates",
            },
        ),
    log:
        "logs/activity/ratio_correlation_between_replicates.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity ratio-correlation-between-replicates --activity-per-rep {input.activity_per_rep} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_ratio_correlation_with_controls:
    input:
        controls=lookup(
            within=activity_files,
            query="file == 'control_df'",
            cols=["path"],
        ),
        activity_per_rep=lookup(
            within=activity_files,
            query="file == 'activity_per_rep'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_ratio_correlation_with_controls_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/Correlation_between_replicates_controls.png",
            caption=getReport("activity/Correlation_between_replicates_controls.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Activity",
                "type": "Reproducibility",
                "figure": "Correlation with Controls",
            },
        ),
    log:
        "logs/activity/ratio_correlation_with_controls.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity ratio-correlation-with-controls --activity-per-rep {input.activity_per_rep} --controls {input.controls} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_downsampling:
    input:
        downsampling_activity_path=lookup(
            within=activity_files,
            query="file == 'downsampling_activity_path'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_downsampling_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/Activity_by_sequencing_depth.png",
            caption=getReport("activity/Activity_by_sequencing_depth.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Downsampling",
                "type": "Activity",
                "figure": "Activity by Sequencing Depth",
            },
        ),
    log:
        "logs/activity/downsampling.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity downsampling --downsampling-activity-path {input.downsampling_activity_path} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_reproducibility_by_sequencing_depth:
    input:
        downsampling_ratio_path=lookup(
            within=activity_files,
            query="file == 'downsampling_ratio_path'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_reproducibility_by_sequencing_depth_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/BC_retention_by_DNA_RNA_sequencing_depth.png",
            caption=getReport("activity/BC_retention_by_DNA_RNA_sequencing_depth.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Downsampling",
                "type": "Retention",
                "figure": "BC Retention by Sequencing Depth",
            },
        ),
        report(
            "results/{project}/activity/cCRE_retention_by_DNA_RNA_sequencing_depth.png",
            caption=getReport("activity/cCRE_retention_by_DNA_RNA_sequencing_depth.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Downsampling",
                "type": "Retention",
                "figure": "cCRE Retention by Sequencing Depth",
            },
        ),
    log:
        "logs/activity/reproducibility_by_sequencing_depth.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity reproducibility-by-sequencing-depth --downsampling-ratio-path {input.downsampling_ratio_path} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_mimimise_noise:
    input:
        sdt_thresholds=lookup(
            within=activity_files,
            query="file == 'different_std_threshold_analysis'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_mimimise_noise_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/Minimizing_noise.png",
            caption=getReport("activity/Minimizing_noise.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Outliers",
                "type": "Reproducibility",
                "figure": "Minimizing Noise",
            },
        ),
    log:
        "logs/activity/minimise_noise.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity mimimise-noise --sdt-thresholds {input.sdt_thresholds} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_screen_annotations:
    input:
        screen=lookup(
            within=activity_files,
            query="file == 'screen_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_screen_annotations_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/Genomic_annotations.png",
            caption=getReport("activity/Genomic_annotations.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Annotations",
                "type": "Activity",
                "figure": "Genomic",
            },
        ),
    log:
        "logs/activity/screen_annotations.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity screen-annotations --screen {input.screen} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_tss_proximity:
    input:
        tss_distance=lookup(
            within=activity_files,
            query="file == 'tss_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_tss_proximity_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/Proximity_to_TSS.png",
            caption=getReport("activity/Proximity_to_TSS.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Annotations",
                "type": "Activity",
                "figure": "Proximity to TSS",
            },
        ),
    log:
        "logs/activity/tss_proximity.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity tss-proximity --tss-distance {input.tss_distance} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_prediction_vs_activity:
    input:
        activity_prediction=lookup(
            within=activity_files,
            query="file == 'AI_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_prediction_vs_activity_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            [
                "results/{project}/activity/AI_predictions_vs_activity.png",
                "results/{project}/activity/AI_predictions_vs_activity_w_bar.png",
            ],
            caption=getReport("activity/AI_predictions_vs_activity.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Prediction",
                "type": "Activity",
                "figure": "Predictions vs Activity",
            },
        ),
    log:
        "logs/activity/prediction_vs_activity.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity prediction-vs-activity --activity-prediction {input.activity_prediction} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_prediction_vs_differential_activity:
    input:
        differential_activity_prediction=lookup(
            within=activity_files,
            query="file == 'AI_comparative_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_prediction_vs_differential_activity_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            [
                "results/{project}/activity/AI_predictions_vs_differential_activity.png",
                "results/{project}/activity/AI_predictions_vs_differential_activity_w_bar.png",
            ],
            caption=getReport("activity/AI_predictions_vs_differential_activity.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Prediction",
                "type": "Comparative",
                "figure": "Predictions vs Differential Activity",
            },
        ),
    log:
        "logs/activity/prediction_vs_differential_activity.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
        python {input.script} activity prediction-vs-differential-activity \
        --differential-activity-prediction {input.differential_activity_prediction} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_comparative:
    input:
        differential_activity=lookup(
            within=activity_files,
            query="file == 'comparative_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_comparative_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            [
                "results/{project}/activity/Volcano_plot_FC_vs_Pval_zoom.png",
                "results/{project}/activity/Volcano_plot_FC_vs_Pval.png",
            ],
            caption=getReport("activity/Volcano_plot_FC_vs_Pval.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Comparative",
                "type": "P-Value",
                "figure": "Volcano",
            },
        ),
        report(
            "results/{project}/activity/Differential_activity_distribution.png",
            caption=getReport("activity/Differential_activity_distribution.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Comparative",
                "type": "Activity",
                "figure": "Differential activity distribution",
            },
        ),
    log:
        "logs/activity/comparative.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
        python {input.script} activity comparative \
        --differential-activity {input.differential_activity} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_allelic_pairs:
    input:
        allelic_pairs=lookup(
            within=activity_files,
            query="file == 'allelic_pairs_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_allelic_pairs_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            [
                "results/{project}/activity/Cross_validation_allelic_pairs.png",
                "results/{project}/activity/Cross_validation_allelic_pairs_w_bar.png",
            ],
            caption=getReport("activity/Cross_validation_allelic_pairs.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Comparative",
                "type": "Reproducibility",
                "figure": "Cross validation: allelic pairs",
            },
        ),
    log:
        "logs/activity/allelic_pairs.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
        python {input.script} activity allelic-pairs \
        --allelic-pairs {input.allelic_pairs} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_cell_types:
    input:
        cell_types=lookup(
            within=activity_files,
            query="file == 'cell_types_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_cell_types_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            [
                "results/{project}/activity/Cross_validation_cell_types.png",
                "results/{project}/activity/Cross_validation_cell_types_w_bar.png",
            ],
            caption=getReport("activity/Cross_validation_cell_types.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Cross validation",
                "type": "Reproducibility",
                "figure": "Cross validation: cell types",
            },
        ),
    log:
        "logs/activity/cell_types.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
        python {input.script} activity cell-types --cell-types {input.cell_types} --output-path {params.outdir} > {log} 2>&1
        """


rule activity_comparative_replicates:
    input:
        differential_activity_replicates=lookup(
            within=activity_files,
            query="file == 'allelic_pairs_replicates_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_comparative_replicates_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            [
                "results/{project}/activity/Correlation_of_differential_activity_between_replicates.png",
                "results/{project}/activity/Correlation_of_differential_activity_between_replicates_w_bar.png",
            ],
            caption=getReport("activity/Correlation_of_differential_activity_between_replicates.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Comparative",
                "type": "Reproducibility",
                "figure": "Correlation of differential activity",
            },
        ),
    log:
        "logs/activity/comparative_replicates.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
        python {input.script} activity comparative-replicates \
        --differential-activity-replicates {input.differential_activity_replicates} \
        --output-path {params.outdir} > {log} 2>&1
        """


rule activity_sample_clusters:
    input:
        reads_by_group=lookup(
            within=activity_files,
            query="file == 'reads_by_group'",
            cols=["path"],
        ),
        sample_metadata=lookup(
            within=activity_files,
            query="file == 'samples_metadata'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_sample_clusters_plots(activity_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/activity/Sample_clustering.png",
            caption=getReport("activity/Sample_clustering.rst"),
            category="{project}",
            subcategory="Activity",
            labels={
                "analysis": "Cross validation",
                "type": "Reproducibility",
                "figure": "Sample clusters",
            },
        ),
    log:
        "logs/activity/sample_clusters.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
        python {input.script} activity sample-clusters \
        --reads-by-group {input.reads_by_group} --sample-metadata {input.sample_metadata} \
        --output-path {params.outdir} > {log} 2>&1
        """

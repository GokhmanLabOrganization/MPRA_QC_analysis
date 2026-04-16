rule association_final:
    input:
        design=lookup(within=association_files, query="file == 'cCRE_fasta'", cols=["path"]),
        associations=lookup(
            within=association_files,
            query="file == 'final_associations'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/association/{plot}.{file_type}",
            plot=get_association_final_plots(association_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/association/BCs_per_cCRE.png",
            caption=getReport("association/BCs_per_cCRE.rst"),
            category="{project}",
            subcategory="Association",
            labels={
                "analysis": "Retention",
                "type": "Barcodes",
                "figure": "Barcodes per cCRE",
            },
        ),
        report(
            "results/{project}/association/Retained_cCREs.png",
            caption=getReport("association/Retained_cCREs.rst"),
            category="{project}",
            subcategory="Association",
            labels={
                "analysis": "Retention",
                "type": "cCREs",
                "figure": "Retained cCREs",
            },
        ),
        report(
            "results/{project}/association/PCR_bias_GC.png",
            caption=getReport("association/PCR_bias_GC.rst"),
            category="{project}",
            subcategory="Association",
            labels={
                "analysis": "Bias",
                "type": "PCR",
                "figure": "GC content",
            },
        ),
        report(
            "results/{project}/association/PCR_bias_G_stretches.png",
            caption=getReport("association/PCR_bias_G_stretches.rst"),
            category="{project}",
            subcategory="Association",
            labels={
                "analysis": "Bias",
                "type": "PCR",
                "figure": "G stretches",
            },
        ),
    log:
        "logs/association/final.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/association/"),
    shell:
        """
         python {input.script}  association final --design {input.design} --associations {input.associations} --output-path {params.outdir} > {log} 2>&1
        """


rule association_before_minimum_observations:
    input:
        associations=lookup(
            within=association_files,
            query="file == 'associations_before_minimum_observations'",
            cols=["path"],
        ),
        design=lookup(within=association_files, query="file == 'cCRE_fasta'", cols=["path"]),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/association/{plot}.{file_type}",
            plot=get_association_before_minimum_observations_plots(association_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/association/Reads_per_association.png",
            caption=getReport("association/Reads_per_association.rst"),
            category="{project}",
            subcategory="Association",
            labels={
                "analysis": "Retention",
                "type": "cCREs",
                "figure": "Reads per association",
            },
        ),
    log:
        "logs/association/before_minimum_observations.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/association/"),
        design=lambda wc, input: f"--design {input.design[0]}" if input.design else "",
    shell:
        """
         python {input.script} association before-minimum-observations --associations {input.associations} {params.design} --output-path {params.outdir} > {log} 2>&1
        """


rule association_before_promiscuity:
    input:
        associations=lookup(
            within=association_files,
            query="file == 'associations_before_promiscuity'",
            cols=["path"],
        ),
        design=lookup(within=association_files, query="file == 'cCRE_fasta'", cols=["path"]),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/association/{plot}.{file_type}",
            plot=get_association_before_promiscuity_plots(association_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/association/cCREs_per_BC.png",
            caption=getReport("association/cCREs_per_BC.rst"),
            category="{project}",
            subcategory="Association",
            labels={
                "analysis": "Retention",
                "type": "cCREs",
                "figure": "cCREs per BC",
            },
        ),
    log:
        "logs/association/before_promiscuity.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/association/"),
        design=lambda wc, input: f"--design {input.design[0]}" if input.design else "",
    shell:
        """
         python {input.script}  association before-promiscuity --associations {input.associations} {params.design} --output-path {params.outdir} > {log} 2>&1
        """


rule association_downsampling:
    input:
        design=lookup(within=association_files, query="file == 'cCRE_fasta'", cols=["path"]),
        downsampling_path=lookup(
            within=association_files,
            query="file == 'associations_downsampling_path'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/association/{plot}.{file_type}",
            plot=get_association_downsampling_plots(association_files),
            file_type=["pdf", "eps", "svg"],
        ),
        report(
            "results/{project}/association/Downsampling_Retained_cCREs.png",
            caption=getReport("association/Downsampling_Retained_cCREs.rst"),
            category="{project}",
            subcategory="Association",
            labels={
                "analysis": "Downsampling",
                "type": "cCREs",
                "figure": "Retained cCREs",
            },
        ),
        report(
            "results/{project}/association/Downsampling_Barcodes_per_cCRE.png",
            caption=getReport("association/Downsampling_Barcodes_per_cCRE.rst"),
            category="{project}",
            subcategory="Association",
            labels={
                "analysis": "Downsampling",
                "type": "Barcodes",
                "figure": "Barcodes per cCRE",
            },
        ),
    log:
        "logs/association/downsampling.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/association/"),
    shell:
        """
         python {input.script}  association downsampling --design {input.design} --downsampling-path {input.downsampling_path} --output-path {params.outdir} > {log} 2>&1
        """

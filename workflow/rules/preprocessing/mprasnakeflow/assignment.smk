rule preprocessing_mprasnakeflow_assignment_associations_before_promiscuity:
    input:
        barcodes_incl_other=config.get("mprasnakeflow", {}).get("assignment", {}).get("barcodes_incl_other", []),
    output:
        "results/{project}/preprocessing/mprasnakeflow/association/associations_before_promiscuity.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/assignment/associations_before_promiscuity.{project}.log",
    conda:
        getCondaEnv("default.yml")
    shell:
        """
        (
            echo "barcode,cCRE,match_count"; 
            zcat {input.barcodes_incl_other} | egrep -v "\\sother\\sNA" | \
            cut -f 1,2 | uniq -c | awk -v "OFS=," '{{print $2,$3,$1}}';
        ) | gzip -c > {output} 2> {log}
        """


rule preprocessing_mprasnakeflow_assignment_associations_before_minimum_observations:
    input:
        association_before_promiscuity="results/{project}/preprocessing/mprasnakeflow/association/associations_before_promiscuity.csv.gz",
        script=getScript("preprocessing/mprasnakeflow/create_before_minimum_associations.py"),
    output:
        "results/{project}/preprocessing/mprasnakeflow/association/associations_before_minimum_observations.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/assignment/associations_before_minimum_observations.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        threshold=config.get("mprasnakeflow", {}).get("assignment", {}).get("fraction", 0.75),
    shell:
        """
        python {input.script} --input {input.association_before_promiscuity} --output {output} --threshold {params.threshold} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_assignment_final_associations:
    input:
        assignment=config.get("mprasnakeflow", {}).get("assignment", {}).get("assignment_barcodes_with_ambiguous", []),
    output:
        "results/{project}/preprocessing/mprasnakeflow/association/final_associations.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/assignment/final_associations.{project}.log",
    conda:
        getCondaEnv("default.yml")
    shell:
        """
        (
            echo "barcode,cCRE,match_count";
            zcat {input.assignment} | egrep -v "\\sother\\s" |  egrep -v "\\sambiguous\\s" | cut -f 1,2,4 | tr '\\t' ',' | sed 's|/[^/]*$||';
        ) | gzip -c > {output} 2> {log}
        """


rule preprocessing_mprasnakeflow_assignment_downsample:
    input:
        assignment=config.get("mprasnakeflow", {}).get("assignment", {}).get("barcodes_incl_other", []),
        script=getScript("preprocessing/mprasnakeflow/assignment_downsample.py"),
    output:
        temp(
            expand(
                "results/{{project}}/preprocessing/mprasnakeflow/association/downsample_tmp/barcodes_incl_other.{fraction:.1f}.tsv.gz",
                fraction=[0.1 * i for i in range(1, 11)],
            )
        ),
    log:
        "logs/preprocessing/mprasnakeflow/assignment/downsample.{project}.log",
    conda:
        getCondaEnv("default.yml")
    shell:
        """
        paths=( {output} );
        python {input.script} --input {input.assignment} \
        --output-folder $(dirname "${{paths[0]}}") > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_assignment_downsample_filter:
    """
Filter the barcodes file based on the config given in the config-file.
FIXME: Limitation is that oligos cannot have a name ambiguous or other.
"""
    input:
        assignment="results/{project}/preprocessing/mprasnakeflow/association/downsample_tmp/barcodes_incl_other.{fraction}.tsv.gz",
        script=getScript("preprocessing/mprasnakeflow/filterAssignmentTsv.py"),
    output:
        final=temp(
            "results/{project}/preprocessing/mprasnakeflow/association/downsample_tmp/associations_final_{fraction}.csv.gz"
        ),
    log:
        log="logs/preprocessing/mprasnakeflow/assignment/downsample_filter.{project}.{fraction}.log",
        err="logs/preprocessing/mprasnakeflow/assignment/downsample_filter.{project}.{fraction}.err",
    conda:
        getCondaEnv("default.yml")
    params:
        min_support=config.get("mprasnakeflow", {}).get("assignment", {}).get("min_support", 3),
        fraction=config.get("mprasnakeflow", {}).get("assignment", {}).get("fraction", 0.75),
        bc_length=config.get("mprasnakeflow", {}).get("assignment", {}).get("bc_length", 15),
    shell:
        """
        trap "cat {log.err}" ERR;

        (
            echo "barcode,cCRE,match_count";
            zcat  {input.assignment} | \
            awk -v "OFS=\\t" -F"\\t" '{{if (length($1)=={params.bc_length}){{print $0 }}}}' | \
            python {input.script} \
            -m {params.min_support} -f {params.fraction} | \
            awk -v "OFS=,"  -F"\\t" '{{ print $1,$2 }}';
        ) | gzip -c > {output.final} 2> {log.err};

        gzip -l {output.final} | awk 'NR==2 {{exit($2==0)}}' || {{ echo "Error: Empty barcode file {output.final}. No barcodes detected!" >> {log.err}; exit 1; }}
        """


rule preprocessing_mprasnakeflow_assignment_downsample_copy:
    """
Copy to the final downsample folder. This is necessary to avoid issues with the downstream rules that expect the files to be in a specific location. The downsample folder is also added to the association_files dataframe for the downstream rules.
"""
    input:
        final=expand(
            "results/{{project}}/preprocessing/mprasnakeflow/association/downsample_tmp/associations_final_{fraction:.1f}.csv.gz",
            fraction=[0.1 * i for i in range(1, 11)],
        ),
    output:
        output_path=directory("results/{project}/preprocessing/mprasnakeflow/association/downsample/"),
    log:
        "logs/preprocessing/mprasnakeflow/assignment/downsample_copy.{project}.log",
    conda:
        getCondaEnv("default.yml")
    shell:
        """
        mkdir -p {output.output_path}
        for i in {input.final}; do
            cp $i {output.output_path}/$(basename $i);
        done 2> {log}
        """

"""Rule wrapping the variantcentrifuge CLI."""

rule run_variantcentrifuge:
    input:
        vcf=lambda wc: get_vcf_path(wc.sample, VCF_FOLDER, samples_df),
    output:
        tsv=os.path.join(OUTPUT_DIR, "{sample}", "{sample}.vc_analysis.tsv"),
    params:
        output_dir=lambda wc: os.path.join(OUTPUT_DIR, wc.sample),
        config_file=VC_CONFIG["config_file"],
        genes=VC_CONFIG["genes"],
        log_level=VC_CONFIG["log_level"],
        extra_flags=" ".join(build_vc_flags(VC_CONFIG, IGV_CONFIG)),
    log:
        os.path.join(OUTPUT_DIR, LOG_SUBDIR, "{sample}.variantcentrifuge.log"),
    container:
        CONTAINER_IMAGE if CONTAINER_IMAGE else None
    conda:
        "../envs/variantcentrifuge.yml"
    threads:
        VC_CONFIG["threads"]
    retries: 1
    resources:
        mem_mb=lambda wc, attempt: 16000 * attempt,
        runtime=lambda wc, attempt: 1440 * attempt,
    shell:
        """
        mkdir -p {params.output_dir}
        variantcentrifuge \
            -c {params.config_file} \
            -g {params.genes} \
            -v {input.vcf} \
            --output-file {output.tsv} \
            --output-dir {params.output_dir} \
            --log-level {params.log_level} \
            --threads {threads} \
            {params.extra_flags} \
            > {log} 2>&1
        """

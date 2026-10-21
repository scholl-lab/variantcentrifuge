"""Rule wrapping the variantcentrifuge CLI."""


# Download the snpEff database once: parallel jobs would each auto-download it into the shared
# env and read it half-written (genes2bed then returns, and caches, a truncated gene BED).
rule snpeff_download_db:
    output:
        touch(os.path.join(OUTPUT_DIR, LOG_SUBDIR, f"snpeff_db.{SNPEFF_DB}.done")),
    log:
        os.path.join(OUTPUT_DIR, LOG_SUBDIR, f"snpeff_download_db.{SNPEFF_DB}.log"),
    params:
        db=SNPEFF_DB,
    container:
        CONTAINER_IMAGE if CONTAINER_IMAGE else None
    conda:
        "../envs/variantcentrifuge.yml"
    shell:
        # a read-only, pre-populated cache (docs/source/installation.md) cannot download: check it loads instead
        "(snpEff download -v {params.db:q} || snpEff genes2bed {params.db:q} > /dev/null) > {log:q} 2>&1"


rule run_variantcentrifuge:
    input:
        vcf=lambda wc: get_vcf_path(wc.sample, VCF_FOLDER, samples_df),
        snpeff_db=rules.snpeff_download_db.output if SNPEFF_DB else [],
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
        mkdir -p -- {params.output_dir:q}
        variantcentrifuge \
            -c {params.config_file:q} \
            -g {params.genes:q} \
            -v {input.vcf:q} \
            --output-file {output.tsv:q} \
            --output-dir {params.output_dir:q} \
            --log-level {params.log_level:q} \
            --threads {threads} \
            {params.extra_flags} \
            > {log:q} 2>&1
        """

'''
Unify format of raw FASTQ files,
All FASTQ will be gzipped and concated before sending to fastp.
'''


rule concat_fastq_paired:
    input:
        lambda w: config["read_mapping_jobs"]["paired"][w.sample_id][w.end],
    output:
        temp("resources/concat_fastq/paired/{sample_id}.{end}.fastq.gz"),
    log:
        "resources/concat_fastq/paired/{sample_id}.{end}.log"
    resources:
        io = 50
    run:
        # If only one gzipped file, soft link it to save IO.
        if len(input) == 1 and input[0].lower().endswith(".gz"):
            shell("ln -s $(realpath {input[0]}) {output} 2> {log}")
            return()

        for f in input:
            if f.endswith(".gz"):
                shell("cat {f} >> {output} 2>> {log}")
            elif f.endswith(".bz2"):
                shell("bzip2 -dc {f} 2>> {log} | gzip -c >> {output} 2>> {log}")
            elif f.endswith(".fasta") or f.endswith(".fa"):
                shell("gzip -c {f} >> {output} 2>> {log}")


rule concat_fastq_single:
    input:
        lambda w: config["read_mapping_jobs"]["single"][w.sample_id],
    output:
        temp("resources/concat_fastq/single/{sample_id}.single.fastq.gz"),
    log:
        "resources/concat_fastq/single/{sample_id}.single.log"
    resources:
        io = 50
    run:
        # If only one gzipped file, soft link it to save IO.
        if len(input) == 1 and input[0].lower().endswith(".gz"):
            shell("ln -s $(realpath {input[0]}) {output} 2> {log}")
            return()

        for f in input:
            if f.endswith(".gz"):
                shell("cat {f} >> {output} 2>> {log}")
            elif f.endswith(".bz2"):
                shell("bzip2 -dc {f} 2>> {log} | gzip -c >> {output} 2>> {log}")
            elif f.endswith(".fasta") or f.endswith(".fa"):
                shell("gzip -c {f} >> {output} 2>> {log}")


'''
Use fastp to filter FASTQ, and generate report html.
'''


rule fastp_paired:
    input:
        r1 = "resources/concat_fastq/paired/{sample_id}.r1.fastq.gz",
        r2 = "resources/concat_fastq/paired/{sample_id}.r2.fastq.gz",
    output:
        r1 = temp("resources/fastp/paired/{sample_id}.fastp.r1.fastq.gz"),
        r2 = temp("resources/fastp/paired/{sample_id}.fastp.r2.fastq.gz"),
        json = "results/fastp/paired/{sample_id}.fastp.json",
        html = "results/fastp/paired/{sample_id}.fastp.html"
    log:
        "results/fastp/paired/{sample_id}.fastp.log"
    threads:
        4
    priority:
        60
    resources:
        io = 100
    shell:
        '''
        fastp \
            --thread {threads} \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample_id} \
            --verbose 2>> {log}
        '''


rule fastp_single:
    input:
        se = "resources/concat_fastq/single/{sample_id}.single.fastq.gz",
    output:
        se = temp("resources/fastp/single/{sample_id}.fastp.single.fastq.gz"),
        json = "results/fastp/single/{sample_id}.fastp.single.json",
        html = "results/fastp/single/{sample_id}.fastp.single.html"
    log:
        "results/fastp/single/{sample_id}.fastp.log"
    threads:
        4
    priority:
        60
    resources:
        io = 100
    shell:
        '''
        fastp \
            --thread {threads} \
            --in1 {input.se} \
            --out1 {output.se} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample_id} \
            --verbose 2>> {log}
        '''


# STAR's built-in sorting is ineffcient and error prone,
# so I pipe output to samtools to do sorting.


rule star_paired:
    input:
        r1 = branch(
            condition = config["fastp"]["skip"],
            then = rules.fastp_paired.input.r1,
            otherwise = rules.fastp_paired.output.r1
        ),
        r2 = branch(
            condition = config["fastp"]["skip"],
            then = rules.fastp_paired.input.r2,
            otherwise = rules.fastp_paired.output.r2
        ),
        star_index = config["star"]["index"],
    output:
        bam = "results/star_mapping/paired/{sample_id}.bam",
        reads_per_gene = "results/star_mapping/paired/{sample_id}.ReadsPerGene.out.tab"
    log:
        "results/star_mapping/paired/{sample_id}.star_mapping.log",
    threads:
        config["star"]["mapping_threads"]
    params:
        sort_threads = config["samtools"]["sort"]["threads"],
        sort_mem_per_thread = config["samtools"]["sort"]["mem_per_thread"],
        out_prefix = "results/star_mapping/paired/{sample_id}.",
    shell:
        '''
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.star_index} \
            --readFilesCommand zcat \
            --readFilesIn {input.r1} {input.r2} \
            --outFileNamePrefix {params.out_prefix} \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --quantMode GeneCounts \
            --outStd BAM_Unsorted 2>> {log} \
        | samtools fixmate -u -m - - 2>> {log} \
        | samtools sort -u -@ {params.sort_threads} -m {params.sort_mem_per_thread} 2>> {log} \
        | samtools markdup -u - - 2>> {log} \
        | samtools view -b -h -@{threads} -o {output.bam} 2>> {log}
        '''


rule star_single:
    input:
        se = branch(
            condition = config["fastp"]["skip"],
            then = rules.fastp_single.input.se,
            otherwise = rules.fastp_single.output.se
        ),
        star_index = config["star"]["index"],
    output:
        bam = "results/star_mapping/single/{sample_id}.bam",
        reads_per_gene = "results/star_mapping/single/{sample_id}.ReadsPerGene.out.tab"
    log:
        "results/star_mapping/single/{sample_id}.star_mapping.log",
    threads:
        config["star"]["mapping_threads"]
    params:
        sort_threads = config["samtools"]["sort"]["threads"],
        sort_mem_per_thread = config["samtools"]["sort"]["mem_per_thread"],
        out_prefix = "results/star_mapping/single/{sample_id}.",
    shell:
        '''
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.star_index} \
            --readFilesCommand zcat \
            --readFilesIn {input.se} \
            --outFileNamePrefix {params.out_prefix} \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --quantMode GeneCounts \
            --outStd BAM_Unsorted 2>> {log} \
        | samtools fixmate -u -m - - 2>> {log} \
        | samtools sort -u -@ {params.sort_threads} -m {params.sort_mem_per_thread} 2>> {log} \
        | samtools markdup -u - - 2>> {log} \
        | samtools view -b -h -@{threads} -o {output.bam} 2>> {log}
        '''


rule collect_star_counts:
    input:
        [f"results/star_mapping/{end}/{sample_id}.ReadsPerGene.out.tab" for end in config["read_mapping_jobs"] for sample_id in config["read_mapping_jobs"][end]]
    output:
        "results/star_mapping/collect_star_counts.csv"
    params:
        sample_id_list = [sample_id for end in config["read_mapping_jobs"] for sample_id in config["read_mapping_jobs"][end]]
    script:
        "../scripts/collect_star_counts.py"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    shell:
        "samtools index {input}"


rule flagstats:
    input:
        bam = "{prefix}.bam",
        bai = "{prefix}.bam.bai"
    output:
        "{prefix}.flagstats.txt"
    threads:
        4
    shell:
        "samtools flagstats -@ {threads} {input.bam} > {output}"


rule collect_star_flagstats:
    input:
        [f"results/star_mapping/{end}/{name}.flagstats.txt" for end, names in config["read_mapping_jobs"].items() for name in names]
    output:
        "results/star_mapping/flagstats.tsv"
    script:
        "../scripts/collect_flagstats.R"


rule salmon_quant_paired:
    input:
        r1 = branch(
            condition = config["fastp"]["skip"],
            then = rules.fastp_paired.input.r1,
            otherwise = rules.fastp_paired.output.r1
        ),
        r2 = branch(
            condition = config["fastp"]["skip"],
            then = rules.fastp_paired.input.r2,
            otherwise = rules.fastp_paired.output.r2
        ),
        salmon_index = config["salmon"]["index"],
    output:
        "results/salmon_quant/paired/{sample_id}/quant.sf"
    threads:
        config["salmon"]["quant_threads"]
    log:
        "results/salmon_quant/paired/{sample_id}.log"
    shell:
        '''
        salmon quant \
            -i {input.salmon_index} \
            -l A -p {threads} \
            -1 {input.r1} \
            -2 {input.r2} \
            --seqBias --gcBias \
            --validateMappings \
            -o {output} 2>> {log}
        '''


rule salmon_quant_single:
    input:
        se = branch(
            condition = config["fastp"]["skip"],
            then = rules.fastp_single.input.se,
            otherwise = rules.fastp_single.output.se
        ),
        salmon_index = config["salmon"]["index"],
    output:
        "results/salmon_quant/single/{sample_id}/quant.sf"
    threads:
        config["salmon"]["quant_threads"]
    log:
        "results/salmon_quant/single/{sample_id}.log"
    shell:
        '''
        salmon quant \
            -i {input.salmon_index} \
            -l A -p {threads} \
            -r {input.se} \
            --seqBias --gcBias \
            --validateMappings \
            -o {output} 2>> {log}
        '''

rule collect_salmon_counts:
    input:
        quant_files = [f"results/salmon_quant/{end}/{sample_id}/quant.sf" for end in config["read_mapping_jobs"] for sample_id in config["read_mapping_jobs"][end]],
        gtf = config["gtf"]
    output:
        "results/salmon_quant/collect_salmon_counts.csv"
    params:
        sample_id_list = [sample_id for end in config["read_mapping_jobs"] for sample_id in config["read_mapping_jobs"][end]]
    script:
        "../scripts/collect_salmon_counts.R"


def get_counts(w):
    if config["mapper"] == "star":
        return rules.collect_star_counts.output
    elif config["mapper"] == "salmon":
        return rules.collect_salmon_counts.output


rule link_counts:
    input:
        get_counts
    output:
        "results/counts.csv"
    shell:
        "ln -s $(realpath {input}) {output}"

configfile: "config/reads2star-config.yaml"

# paired end is prefered over single end.
ruleorder: star_mapping_pe > star_mapping_se > fastp_pe > fastp_se


rule fastp_pe:
    input:
        r1 = lambda w: config["samples"][w.sample_id]["R1"],
        r2 = lambda w: config["samples"][w.sample_id]["R2"]
    output:
        r1_fastp = temp("resources/fastp/{sample_id}.R1.fastp.fastq.gz"),
        r2_fastp = temp("resources/fastp/{sample_id}.R2.fastp.fastq.gz"),
        json = "results/fastp/{sample_id}.json",
        html = "results/fastp/{sample_id}.html"
    log:
        "logs/fastp/{sample_id}.log"
    threads:
        config["threads"]["fastp_pe"]
    priority:
        60
    shell:
        '''
        fastp \
            --thread {threads} \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1_fastp} \
            --out2 {output.r2_fastp} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample_id} \
            --verbose 2>> {log}
        '''

rule fastp_se:
    input:
        r1 = lambda w: config["samples"][w.sample_id]["R1"]
    output:
        r1_fastp = temp("resources/fastp/{sample_id}.R1.fastp.fastq.gz"),
        json = "results/fastp/{sample_id}.json",
        html = "results/fastp/{sample_id}.html"
    log:
        "logs/fastp/{sample_id}.log"
    threads:
        config["threads"]["fastp_se"]
    priority:
        60
    shell:
        '''
        fastp \
            --thread {threads} \
            --in1 {input.r1} \
            --out1 {output.r1_fastp} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample_id} \
            --verbose 2>> {log}
        '''

rule star_mapping_pe:
    input:
        r1_fastp = "resources/fastp/{sample_id}.R1.fastp.fastq.gz",
        r2_fastp = "resources/fastp/{sample_id}.R2.fastp.fastq.gz",
        star_index_dir = config["star_index_dir"]
    output:
        bam = "results/star_mapping/{sample_id}.Aligned.sortedByCoord.out.bam",
        bai = "results/star_mapping/{sample_id}.Aligned.sortedByCoord.out.bam.bai",
        reads_per_gene = "results/star_mapping/{sample_id}.ReadsPerGene.out.tab"
    params:
        out_prefix = "results/star_mapping/{sample_id}.",
        star_params = config["star_params"]
    log:
        "logs/star_mapping/{sample_id}.log"
    threads:
        config["threads"]["star_mapping_pe"]
    priority:
        70
    shell:
        '''
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.star_index_dir} \
            --readFilesCommand zcat \
            --readFilesIn \
                {input.r1_fastp} \
                {input.r2_fastp} \
            --outFileNamePrefix {params.out_prefix} \
            --outBAMsortingThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes Standard \
            --outSAMunmapped Within \
            {params.star_params} \
            --quantMode GeneCounts 2>> {log} >> {log}
        samtools index {output.bam} 2>> {log}
        '''

rule star_mapping_se:
    input:
        r1_fastp = "resources/fastp/{sample_id}.R1.fastp.fastq.gz",
        star_index_dir = config["star_index_dir"]
    output:
        bam = "results/star_mapping/{sample_id}.Aligned.sortedByCoord.out.bam",
        bai = "results/star_mapping/{sample_id}.Aligned.sortedByCoord.out.bam.bai",
        reads_per_gene = "results/star_mapping/{sample_id}.ReadsPerGene.out.tab"
    params:
        out_prefix = "results/star_mapping/{sample_id}.",
        star_params = config["star_params"]
    log:
        "logs/star_mapping/{sample_id}.log"
    threads:
        config["threads"]["star_mapping_se"]
    priority:
        70
    shell:
        '''
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.star_index_dir} \
            --readFilesCommand zcat \
            --readFilesIn \
                {input.r1_fastp} \
            --outFileNamePrefix {params.out_prefix} \
            --outBAMsortingThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes Standard \
            --outSAMunmapped Within \
            {params.star_params} \
            --quantMode GeneCounts 2>> {log} >> {log}
        samtools index {output.bam} 2>> {log}
        '''


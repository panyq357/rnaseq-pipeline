
'''
Use STAR to map reads to reference.
'''

ruleorder: star_mapping_pe > star_mapping_se

rule star_mapping_pe:
    input:
        r1 = "resources/cat/{sample_id}.R1.fastq.gz" if config["skip_fastp"] else "resources/fastp/{sample_id}.fastp.R1.fastq.gz",
        r2 = "resources/cat/{sample_id}.R2.fastq.gz" if config["skip_fastp"] else "resources/fastp/{sample_id}.fastp.R2.fastq.gz",
        star_index_dir = config["star_index_dir"]
    output:
        bam = "results/star_mapping/{sample_id}.Aligned.sortedByCoord.out.bam",
        reads_per_gene = "results/star_mapping/{sample_id}.ReadsPerGene.out.tab"
    params:
        out_prefix = "results/star_mapping/{sample_id}.",
        star_params = config["star_params"]
    log:
        "logs/star_mapping/{sample_id}.log"
    threads:
        config["threads"]["star_mapping"]
    priority:
        70
    shell:
        '''
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.star_index_dir} \
            --readFilesCommand zcat \
            --readFilesIn {input.r1} {input.r2} \
            --outFileNamePrefix {params.out_prefix} \
            --outBAMsortingThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes Standard \
            --outSAMunmapped Within \
            {params.star_params} \
            --quantMode GeneCounts 2>> {log} >> {log}
        '''

rule star_mapping_se:
    input:
        se = "resources/cat/{sample_id}.R1.fastq.gz" if config["skip_fastp"] else "resources/fastp/{sample_id}.fastp.fastq.gz",
        star_index_dir = config["star_index_dir"]
    output:
        bam = "results/star_mapping/{sample_id}.Aligned.sortedByCoord.out.bam",
        reads_per_gene = "results/star_mapping/{sample_id}.ReadsPerGene.out.tab"
    params:
        out_prefix = "results/star_mapping/{sample_id}.",
        star_params = config["star_params"]
    log:
        "logs/star_mapping/{sample_id}.log"
    threads:
        config["threads"]["star_mapping"]
    priority:
        70
    shell:
        '''
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.star_index_dir} \
            --readFilesCommand zcat \
            --readFilesIn {input.se} \
            --outFileNamePrefix {params.out_prefix} \
            --outBAMsortingThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes Standard \
            --outSAMunmapped Within \
            {params.star_params} \
            --quantMode GeneCounts 2>> {log} >> {log}
        '''


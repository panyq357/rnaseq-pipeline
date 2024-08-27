
'''
Use STAR to map reads to reference.
'''

ruleorder: star_mapping_pe > star_mapping_se

rule star_mapping_pe:
    input:
        r1_fastp = "resources/fastp/{sample_id}.fastp.R1.fastq.gz",
        r2_fastp = "resources/fastp/{sample_id}.fastp.R2.fastq.gz",
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
        config["threads"]["star_mapping"]
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
        r1_fastp = "resources/fastp/{sample_id}.fastp.R1.fastq.gz",
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
        config["threads"]["star_mapping"]
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


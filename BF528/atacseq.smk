import pandas

sample_csv = pandas.read_csv('sample_sheet.csv', index_col='name')

# REPS = ['ATACrep3', 'ATACrep4']
# READS = ['R1', 'R2']

REPS = set(sample_csv['replicate'].tolist())
READS = set(sample_csv['read'].tolist())
PAIRS=['1P', '1U', '2P', '2U']
BAMS = ['sorted', 'rmchrM']

rule all:
	input:
		'results/qc/multiqc_report.html',
		'results/ATAC_annotated.bed',
		'results/ATAC_motifs',
		expand('results/signal_coverage_{reps}_0_100.png', reps=REPS),
		expand('results/signal_coverage_{reps}_180_247.png', reps=REPS)

rule wget_files:
	output: 'samples/{REPS}_{READS}.fastq.gz'
	params:
		link = lambda wildcards: sample_csv.loc['{}_{}'.format(wildcards.replicate, wildcards.read), 'ftp_link'],
		renamed = lambda wildcards: 'samples/{}_{}.fastq.gz'.format(wildcards.replicate, wildcards.read)
	shell: '''
		wget -O {params.renamed} {params.link} 
		'''

rule unzip_gtf:
	input: 'results/gencode.v45.primary_assembly.annotation.gtf.gz'
	output: 'results/gencode.v45.primary_assembly.annotation.gtf'
	shell: '''
		gzip -d {input}
		'''

rule unzip_genome:
	input: 'results/GRCh38.primary_assembly.genome.fa.gz'
	output: 'results/GRCh38.primary_assembly.genome.fa'
	shell: '''
		gzip -d {input}
		'''

rule bowtie2_build_gencode:
	input: 'results/GRCh38.primary_assembly.genome.fa'
	output: expand('results/bowtie/primary_assembly.{middle}.bt2', middle = ['1', '2', '3', '4', 'rev.1', 'rev.2'])
	params: 'results/bowtie/primary_assembly'
	threads: 32
	conda: 'envs/bowtie2_env.yml'
	shell: '''
		bowtie2-build -f {input} {params}
		'''

rule fastqc:
	input: 'samples/{REPS}_{READS}.fastq.gz'
	output: 'results/qc/{REPS}_{READS}_fastqc.html'
	params: 'results/qc/'
	threads: 4
	conda: 'envs/fastqc_env.yml'
	shell: '''
		fastqc {input} -o {params}
		'''

# rule trimmomatic:
# 	input: 'samples/{REPS}_{READS}.fastq.gz'
# 	output: 'results/{REPS}_{READS}.trimmed.fastq.gz'
# 	params: 'results/NexteraPe-PE.fa'
# 	threads: 8
# 	conda: 'envs/trimmomatic_env.yml'
# 	shell: '''
# 		trimmomatic SE -threads {threads} \
# 		{input} {output} \
# 		ILLUMINACLIP:{params}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
# 		'''

rule trimmomatic:
	input:
		r1 = 'samples/{REPS}_R1.fastq.gz',
		r2 = 'samples/{REPS}_R2.fastq.gz'
	output:
		p1 = 'results/{REPS}_1P.fastq.gz',
		p2 = 'results/{REPS}_2P.fastq.gz',
		u1 = 'results/{REPS}_1U.fastq.gz',
		u2 = 'results/{REPS}_2U.fastq.gz'
	params:
		adapter = 'results/NexteraPE-PE.fa',
		outfile = 'results/{REPS}.fastq.gz'
	threads: 8
	conda: 'envs/trimmomatic_env.yml'
	shell: '''
		trimmomatic PE -threads {threads} \
		{input.r1} {input.r2} -baseout {params.outfile} \
		ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
		'''

rule fastqc_trim:
	input: 'results/{REPS}_{PAIRS}.fastq.gz'
	output: 'results/qc/{REPS}_{PAIRS}_fastqc.html'
	params: 'results/qc/'
	threads: 4
	conda: 'envs/fastqc_env.yml'
	shell: '''
		fastqc {input} -o {params}
		'''

rule bowtie2_align:
	input: 
		r1 = 'results/{REPS}_1P.fastq.gz',
		r2 = 'results/{REPS}_2P.fastq.gz',
		build = expand('results/bowtie/primary_assembly.{middle}.bt2', middle = ['1', '2', '3', '4', 'rev.1', 'rev.2'])
	output: 'results/{REPS}.bam'
	threads: 32
	params: ref_index = 'results/bowtie/primary_assembly'
	conda: 'envs/bowtie2_env.yml'
	shell: '''
		bowtie2 -x {params.ref_index} -1 {input.r1} -2 {input.r2} -X 2000 | \
		samtools view -bS - > {output}
		'''

rule samtools_sort:
	input: 'results/{REPS}.bam'
	output: 'results/{REPS}.sorted.bam'
	conda: 'envs/samtools_env.yml'
	threads: 8
	shell: '''
		samtools sort {input} -o {output}
		'''

rule samtools_idx:
	input: 'results/{REPS}.sorted.bam'
	output: 'results/{REPS}.sorted.bam.bai'
	conda: 'envs/samtools_env.yml'
	threads: 8
	shell: '''
		samtools index -b {input} -o {output}
		'''

rule samtools_remove_chrM:
	input: bam='results/{REPS}.sorted.bam',
		idx='results/{REPS}.sorted.bam.bai'
	output: 'results/{REPS}.rmchrM.bam'
	conda: 'envs/samtools_env.yml'
	threads: 8
	shell: '''
		samtools view -o {output} -e 'rname != "chrM"' {input.bam}
		'''

rule samtools_flagstat:
	input: 'results/{REPS}.{BAMS}.bam'
	output: 'results/qc/{REPS}.{BAMS}.txt'
	conda: 'envs/samtools_env.yml'
	threads: 8
	shell: '''
		samtools flagstat {input} > {output}
		'''

rule multiqc:
	input: expand('results/qc/{reps}.{bams}.txt', reps = REPS, bams = BAMS),
		expand('results/qc/{reps}_{reads}_fastqc.html', reps = REPS, reads = READS),
		expand('results/qc/{reps}_{pairs}_fastqc.html', reps = REPS, pairs = PAIRS)
	output: 'results/qc/multiqc_report.html'
	params: 
			indir = 'results/qc',
			outdir = 'results/qc'
	conda: 'envs/multiqc_env.yml'
	shell: '''
		multiqc {params.indir} -o {params.outdir}
		'''

rule samtools_idx_rmchrM:
	input: 'results/{REPS}.rmchrM.bam'
	output: 'results/{REPS}.rmchrM.bam.bai'
	conda: 'envs/samtools_env.yml'
	threads: 8
	shell: '''
		samtools index -b {input} -o {output}
		'''

rule deeptools_shift:
	input: bam = 'results/{REPS}.rmchrM.bam',
		index = 'results/{REPS}.rmchrM.bam.bai'
	output: 'results/{REPS}.shifted.bam'
	conda: 'envs/deeptools_env.yml'
	threads: 8
	shell: '''
		alignmentSieve -b {input.bam} --ATACshift -o {output}
		'''

rule samtools_idx_sort:
	input: 'results/{REPS}.shifted.bam'
	output: 'results/{REPS}.shiftsort.bam'
	conda: 'envs/samtools_env.yml'
	threads: 8
	shell: '''
		samtools sort {input} -o {output}
		'''

rule samtools_idx_shift:
	input: 'results/{REPS}.shiftsort.bam'
	output: 'results/{REPS}.shiftsort.bam.bai'
	conda: 'envs/samtools_env.yml'
	threads: 8
	shell: '''
		samtools index -b {input} -o {output}
		'''

rule call_peaks:
	input: 'results/{REPS}.shifted.bam'
	output: 'results/{REPS}_summits.bed'
	params:
		name = '{REPS}',
		outdir = 'results'
	conda: 'envs/macs3_env.yml'
	shell: '''
		macs3 callpeak -t {input} -f BAM \
		-n {params.name} --outdir {params.outdir}
		'''	

rule intersect_peaks:
	input:
		a = 'results/ATACrep3_summits.bed',
		b = 'results/ATACrep4_summits.bed'
	output: 'results/ATAC_intersect.bed'
	conda: 'envs/bedtools_env.yml'
	shell: '''
		bedtools intersect -a {input.a} -b {input.b} -f 0.5  > {output}
		'''

rule filter_blacklist:
	input:
		peaks = 'results/ATAC_intersect.bed',
		blacklist = 'results/hg38-blacklist.v2.bed'
	output: 'results/ATAC_filtered.bed'
	conda: 'envs/bedtools_env.yml'
	shell: '''
		bedtools intersect -a {input.peaks} -b {input.blacklist} -v > {output}
		'''

rule annotate_peaks:
	input:
		filtered = 'results/ATAC_filtered.bed',
		annotation = 'results/gencode.v45.primary_assembly.annotation.gtf'
	output: 'results/ATAC_annotated.bed'
	conda: 'envs/homer_env.yml'
	shell: '''
		annotatePeaks.pl {input.filtered} mm39 -gtf {input.annotation} > {output}
		'''

rule motifs:
	input:
		peaks = 'results/ATAC_filtered.bed',
		genome = 'results/GRCh38.primary_assembly.genome.fa'
	output: directory('results/ATAC_motifs')
	conda: 'envs/homer_env.yml'
	shell: '''
		findMotifsGenome.pl {input.peaks} {input.genome} {output} -size 200
		'''

# Coverage plots for NFR/NBR
rule deeptools_filter:
	input: bam = 'results/{REPS}.shiftsort.bam',
		index = 'results/{REPS}.shiftsort.bam.bai'
	output: 'results/{REPS}_{min}_{max}.sep.bam'
	conda: 'envs/deeptools_env.yml'
	shell: '''
		alignmentSieve -b {input.bam} -o {output} \
		--minFragmentLength {wildcards.min} \
		--maxFragmentLength {wildcards.max}
		'''

rule samtools_idx_filter:
	input: 'results/{REPS}_{min}_{max}.sep.bam'
	output: 'results/{REPS}_{min}_{max}.sep.bam.bai'
	conda: 'envs/samtools_env.yml'
	shell: '''
		samtools index -b {input} -o {output}
		'''

rule bamCoverage:
	input: 
		bam = 'results/{REPS}_{min}_{max}.sep.bam',
		index = 'results/{REPS}_{min}_{max}.sep.bam.bai'
	output: 'results/{REPS}_{min}_{max}.bw'
	threads: 4
	conda: 'envs/deeptools_env.yml'
	shell: '''
		bamCoverage -b {input.bam} -o {output}
		'''

rule computeMatrix:
	input:
		atac = 'results/{REPS}_{min}_{max}.bw',
		regions = 'results/hg38_genes.bed'
	output: 'results/matrix_{REPS}_{min}_{max}.gz'
	params:
		upstream = 2000,
		downstream = 2000
	conda: 'envs/deeptools_env.yml'
	threads: 4
	shell: '''computeMatrix reference-point -R {input.regions} \
		-S {input.atac} \
		-o {output} \
		-b {params.upstream} -a {params.downstream}
		'''

rule plotMatrix:
	input: 'results/matrix_{REPS}_{min}_{max}.gz'
	output: 'results/signal_coverage_{REPS}_{min}_{max}.png'
	conda: 'envs/deeptools_env.yml'
	shell: '''
		plotProfile -m {input} -o {output}
		'''

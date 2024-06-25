#!/usr/bin/env Rscript

library('ATACseqQC')

rep <- '/projectnb/bf528/students/dltamayo/bf528-individual-project-dltamayo/results/ATACrep3.shiftsort.bam'
rep.labels <- 'ATACrep3'

png('/projectnb/bf528/students/dltamayo/bf528-individual-project-dltamayo/results/ATACseqQC/ATACrep3_lib_complexity.png')
print(estimateLibComplexity(readsDupFreq(rep)))
dev.off()

png('/projectnb/bf528/students/dltamayo/bf528-individual-project-dltamayo/results/ATACseqQC/ATACrep3_frag_size.png')
print(fragSizeDist(rep, rep.labels))
dev.off()


rep <- '/projectnb/bf528/students/dltamayo/bf528-individual-project-dltamayo/results/ATACrep4.shiftsort.bam'
rep.labels <- 'ATACrep4'

png('/projectnb/bf528/students/dltamayo/bf528-individual-project-dltamayo/results/ATACseqQC/ATACrep4_lib_complexity.png')
print(estimateLibComplexity(readsDupFreq(rep)))
dev.off()

png('/projectnb/bf528/students/dltamayo/bf528-individual-project-dltamayo/results/ATACseqQC/ATACrep4_frag_size.png')
print(fragSizeDist(rep, rep.labels))
dev.off()

###########################################################################
##########   041217_BowtieAlignment_MEJH112315_Bowtie2.sh                 ##########
##########   shell script to redo bowtie alignment - with Bowtie2         ##########
############################################################################
#! /bin/bash
bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 041217Unmapped_ME_JH_112315_1105-ATACSlice2-01-20a_S10_L002_R1_001.fastq -x ~/Indexed_Genomes/dmel-r5.57-index -1 ME_JH_112315_1105-ATACSlice2-01-20a_S10_L002_R1_001.fastq -2 ME_JH_112315_1105-ATACSlice2-01-20a_S10_L002_R2_001.fastq -S 041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-01-20a.sam

bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 041217Unmapped_ME_JH_112315_1105-ATACSlice2-02-20p_S11_L002_R1_001.fastq -x ~/Indexed_Genomes/dmel-r5.57-index -1 ME_JH_112315_1105-ATACSlice2-02-20p_S11_L002_R1_001.fastq -2 ME_JH_112315_1105-ATACSlice2-02-20p_S11_L002_R2_001.fastq -S 041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-02-20p.sam

bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 041217Unmapped_ME_JH_112315_1105-ATACSlice2-5-1A_S12_L002_R1_001.fastq -x ~/Indexed_Genomes/dmel-r5.57-index -1 ME_JH_112315_1105-ATACSlice2-5-1A_S12_L002_R1_001.fastq -2 ME_JH_112315_1105-ATACSlice2-5-1A_S12_L002_R2_001.fastq -S 041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-5-1A.sam

bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 041217Unmapped_ME_JH_112315_1105-ATACSlice2-6-1p_S13_L002_R1_001.fastq -x ~/Indexed_Genomes/dmel-r5.57-index -1 ME_JH_112315_1105-ATACSlice2-6-1p_S13_L002_R1_001.fastq -2 ME_JH_112315_1105-ATACSlice2-6-1p_S13_L002_R2_001.fastq -S 041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-6-1p.sam

bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 041217Unmapped_ME_JH_112315_1105-ATACSlice2-7-10whole_S14_L002_R1_001.fastq -x ~/Indexed_Genomes/dmel-r5.57-index -1 ME_JH_112315_1105-ATACSlice2-7-10whole_S14_L002_R1_001.fastq -2 ME_JH_112315_1105-ATACSlice2-7-10whole_S14_L002_R2_001.fastq -S 041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-7-10whole.sam

bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 041217Unmapped_ME_JH_112315_1109-ATACSlice03-01-20Ant_S6_L002_R1_001.fastq -x ~/Indexed_Genomes/dmel-r5.57-index -1 ME_JH_112315_1109-ATACSlice03-01-20Ant_S6_L002_R1_001.fastq -2 ME_JH_112315_1109-ATACSlice03-01-20Ant_S6_L002_R2_001.fastq -S 041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-01-20Ant.sam

bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 041217Unmapped_ME_JH_112315_1109-ATACSlice03-02-20Post_S7_L002_R1_001.fastq -x ~/Indexed_Genomes/dmel-r5.57-index -1 ME_JH_112315_1109-ATACSlice03-02-20Post_S7_L002_R1_001.fastq -2 ME_JH_112315_1109-ATACSlice03-02-20Post_S7_L002_R2_001.fastq -S 041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-02-20Post.sam

bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 041217Unmapped_ME_JH_112315_1109-ATACSlice03-03-1plus2_S8_L002_R1_001.fastq -x ~/Indexed_Genomes/dmel-r5.57-index -1 ME_JH_112315_1109-ATACSlice03-03-1plus2_S8_L002_R1_001.fastq -2 ME_JH_112315_1109-ATACSlice03-03-1plus2_S8_L002_R2_001.fastq -S 041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-03-1plus2.sam

bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 041217Unmapped_ME_JH_112315_1109-ATACSlice03-04-10whole_S9_L002_R1_001.fastq -x ~/Indexed_Genomes/dmel-r5.57-index -1 ME_JH_112315_1109-ATACSlice03-04-10whole_S9_L002_R1_001.fastq -2 ME_JH_112315_1109-ATACSlice03-04-10whole_S9_L002_R2_001.fastq -S 041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-04-10whole.sam

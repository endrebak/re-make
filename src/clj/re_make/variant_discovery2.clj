


(def samples ["A" "B"])

(def config {:prefix "/User/endrebakkenstovner/snakemake-example"})

(def sample-sheet
  [{:sample "A" :fastq "data/samples/A.fastq"}
   {:sample "B" :fastq "data/samples/B.fastq"}
   {:sample "C" :fastq "data/samples/C.fastq"}])



(def input
  {:genome "data/genome.fa"
   :fastq "data/samples/{sample}.fastq"})

(def output
  {:sorted-bam {:protected true}})


;; automatically looks up in sample-sheet using wildcards
(defrule bwa-map
  "Map DNA sequences against a reference genome with BWA."
  {:wildcards [:sample]
   :input [:genome :fastq]
   :output "bwa-map.bam"
   :threads 8
   :params {:rg "@RG\tID:{sample}\tSM:{sample}"}
   :shell "bwa mem -R '{params.rg}' {threads} {genome} {fastq} | samtools view -Sb - > {bwa-map.bam}"})


(defrule samtools-sort
  "Sort the bams."
  {:wildcards [:sample]
   :input "bwa-map.bam"
   :output "sorted.bam"
   :shell "samtools sort -T {sample} -O bam {bwa-map.bam} > {sorted.bam}"})


(defrule samtools-index
  "Index read alignments for random access."
  {:wildcards [:sample]
   :input "sorted.bam"
   :output "sorted.bam.bai"
   :shell "samtools index {sorted.bam}"})


;; it is so common to collect all of a wildcard that it happens by default
;; yes, but what if required custom logic?
;; think about later
(defrule bcftools-call
  "Aggregate mapped reads from all samples and jointly call genomic variants on
  them."
  {:input [:genome "sorted.bam" "sorted.bam.bai"]
   :output "all.vcf"
   :shell "samtools mpileup -g -f {genome} {sorted.bam} | bcftools call -mv - > {all.vcf}"})


; rulenames are unique
; therefore can use as a key to the script
(defrule plot-quals
  {:input "all.vcf"
   :output "quals.svg"})


;; (def x {:sample ["A" "B"]})

(def xs {:sample ["A" "B" "C"]
         :genome ["hg38" "hg19"]})



;; (def in
;;   {:genome "data/genome.fa"
;;    :fastq "data/samples/{sample}.fastq"})


(def out
  {"bam/sorted.bam" {:protected true}})


;; automatically looks up in sample-sheet using x
(defrule bwa-map
  "Map DNA sequences against a reference genome with BWA."
  {:x [:sample :genome]
   :file [:genome :fastq]
   :out "bwa-map.bam"
   :threads 8
   :params {:rg "@RG\tID:{sample}\tSM:{sample}"}
   :shell "bwa mem -R '{params.rg}' {threads} {ext.genome} {ext.fastq} | samtools view -Sb - > {out}"})


(defrule samtools-sort
  "Sort the bams."
  {:x [:sample :genome]
   :in "bwa-map.bam"
   :out "bam/sorted.bam"
   :shell "samtools sort -T {sample} -O bam {in} > {out}"})


(defrule samtools-index
  "Index read alignments for random access."
  {:x [:sample :genome]
   :in "bam/sorted.bam"
   :out "bam/sorted.bam.bai"
   :shell "samtools index {in}"})


;; it is so common to collect all of a wildcard that it happens by default
;; yes, but what if required custom logic?
;; think about later
(defrule bcftools-call
  "Aggregate mapped reads from all samples and jointly call genomic variants."
  {:in ["bam/sorted.bam" "bam/sorted.bam.bai"]
   :out "all.vcf"
   :x [:genome]
   :xf #(collect :sample %)
   :file :genome
   :sh "samtools mpileup -g -f {file.genome} {in.0} | bcftools call -mv - > {out}"})


; rulenames are unique
; therefore can external as a key to the script
(defrule plot-quals
  {:in "all.vcf"
   :x [:genome]
   :out ["quals.svg" "quals.tsv"]
   :shell "plot {in} -o {out.0}"})

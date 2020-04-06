; no target rules

;; when this compiles, it creates a list of dicts of jobs
;; sort order is execution order

;; by default, output-paths are prefix + "/" + rulename + "/" + "/" + wildcards
;; a log is also always created by default at same paths
;; same with benchmark; on by default


(def samples ["A" "B"])

(def config {:prefix "/User/endrebakkenstovner/snakemake-example"})

(def sample-sheet
  {[{:sample "A" :fastq "data/samples/A.fastq"}
    {:sample "B" :fastq "data/samples/B.fastq"}
    {:sample "C" :fastq "data/samples/C.fastq"}]})

;; automatically looks up in sample-sheet using wildcards

(defrule bwa-map
  "Map DNA sequences against a reference genome with BWA."
  {:wildcards [:sample]
   :input {:genome "data/genome.fa"
           :fastq sample-sheet}
   :output ".bam"
   :threads 8
   :params {:rg "@RG\tID:{sample}\tSM:{sample}"}
   :shell "bwa mem -R '{params.rg}' {threads} {input.genome} {input.fastq} | samtools view -Sb - > {output}"})

;; problem: how do you set the same path for outputs?

(def bam-files "")

(defrule samtools-sort
  "Sort the bams."
  {:wildcards [:sample]
   :input bwa-map
   :output {:ext ".bam" :path bam-files}
   :shell "samtools sort -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"})


(defrule samtools-index
  "Index read alignments for random access."
  {:wildcards [:sample]
   :input samtools-sort
   :output {:ext ".bam.bai" :path bam-files :flag :protected}
   :shell "samtools index {input}"})

;; it is so common to collect all of a wildcard that it happens by default
;; yes, but what if required custom logic?
;; think about later
(defrule bcftools-call
  "Aggregate mapped reads from all samples and jointly call genomic variants on
  them."
  {:input {:fa "data/genome.fa"
           :bam samtools-sort
           :bai samtools-index}
   :output ".vcf"
   :shell "samtools mpileup -g -f {input.fa} {input.bam} | bcftools call -mv - > {output}"})

; rulenames are unique
; therefore can use as a key to the script

(defrule plot-quals
  {:input bcftools-call
   :output ".svg"})



;; (defrule 'create-bed
;;   {:run
;;    "#!/usr/bin/env python

;;    import pyranges as pr
;;    gr = pr.random()
;;    gr.to_csv({output[0]}, sep='\t')"})


;; (defrule 'liftover
;;   {:input 'create-bed
;;    :run
;;    "cp {input} {output}"})


;; (defrule 'lengths
;;   {:input 'liftover
;;    :run
;;    "#!/usr/bin/env python

;;    import pyranges as pr
;;    gr = pr.read_bed(input)
;;    gr = gr[gr.lengths() >= 100]
;;    gr.to_csv({output})
;;    "})


;; (defrule 'epic2
;;   {:input {:bed '[lengths create-bed]}
;;    :wildcards '[genome]
;;    :input-switch #(case (:genome wildcards)
;;                     "hg38" 'lengths
;;                     "hg19" 'create-bed)
;;    :run
;;    "cp {input} {output}"})

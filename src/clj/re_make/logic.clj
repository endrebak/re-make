(ns re-make.logic
  )

;; need a collector

;; need endpoints

;; need to hi hi

;; (def getme)

;; by explicitly marking whether you want the input to be a file not created by
;; re-make you can check whether your workflow has DAG errors

(def target :bwa-map)

(def wildcards
  {:samples ["A" "B" "C"]})
  ;; {:chromosomes (map #(str "chr" %) (range 1 23))})

(def bwa-map
  {:input-files {:genome "/Users/endrebakkenstovner/code/re-make/snakemake-tutorial-data-5.4.5/data/genome.fa"
                 :sample "/Users/endrebakkenstovner/code/re-make/snakemake-tutorial-data-5.4.5/data/samples/{sample}.fa"}
   :output "/Users/endrebakkenstovner/deleteme/samples/{sample}.bam"
   :shell "bwa mem {input.genome} {input.sample} | samtools view -Sb - > {output}"})


;; (def )
;; (def fetch-variants
;;   {:wildcards [:chromosomes]
;;    :publish "{chromosome}.vcf.gz"
;;    :shell "axel {url} -q -o {output}"})


;; (def create-index
;;   {:input fetch-variants
;;    :wildcards [:chromosomes]
;;    :publish "{chromosome}.vcf.gz.tbi"
;;    :shell "tabix {input[0]}"})


;; (def split-into-chromosomes
;;   {:input download-genome
;;     :wildcards [:chromosomes]
;;     :shell "blah"})


;; (defn moo-func [args]
;;   (println "hi wildcard1 wildcard2"))


;; (def moo {:input hi
;;           :fn moo-func
;;           :publish "tmp/"})

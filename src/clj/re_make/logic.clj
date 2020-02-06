(ns re-make.logic
  )

;; need a collector

;; need endpoints

;; need to hi hi

;; (def getme)


(def wildcards
  {:chromosomes (map #(str "chr" %) (range 1 23))})

(def fetch-variants
  {:wildcards [:chromosomes]
   :publish "{chromosome}.vcf.gz"
   :shell "axel {url} -q -o {output}"})

(def create-index
  {:input fetch-variants
   :wildcards [:chromosomes]
   :publish "{chromosome}.vcf.gz.tbi"
   :shell "tabix {input[0]}"})

(def split-into-chromosomes
  {:input download-genome
    :wildcards [:chromosomes]
    :shell "blah"})


(def compress-chromosomes
  {})


(defn moo-func [args]
  (println "hi wildcard1 wildcard2"))


(def moo {:input hi
          :fn moo-func
          :publish "tmp/"})

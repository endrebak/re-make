(ns re-make.read-workflows
  (:import [java.io PushbackReader])
  (:require
   [re-make.state :as state]
   [com.stuartsierra.dependency :as dep]
   [clojure.java.io :as io]))




;; (defn defrule [rulename rulebody]
;;   (defrule! rulename rulebody))


;; https://stackoverflow.com/a/24922859/992687
(defn read-all
  [file]
  (let [rdr (-> file io/file io/reader PushbackReader.)]
    (loop [forms []]
      (let [form (try (read rdr) (catch Exception e nil))]
        (if form
          (recur (conj forms form))
          forms)))))


(defn rules [f]
  (nth (read-all f) 0))


(def rules (atom {}))

(defn handle-docs [& body]
  (if (string? (first body))
    (assoc (second body) :doc (first body))
    (first body)))
; (handle-docs "docstring" {:bla "bla"}) -> {:bla "bla", :doc "docstring"}
; (handle-docs {:bla :bla}) -> {:bla :bla}

(defmacro defrule
  [name & body]
  `(do
     (let [body# (handle-docs ~@body)]
       (swap! rules assoc ~(keyword name) (handle-docs body#))
       (def ~(vary-meta name assoc :rule true) (handle-docs body#)))))

(defn name->dependencies
  [rules]
  (for [[name v] rules
        :let [input (:input v)]]
    [name input]))

(defn rulegraph
  [rules]
  (let [m (name->dependencies rules)]
    (for [[n v] m]
      (if (keyword? v)
        [n v]
        [n (vals v)]))))

(def r
  {:bwa-map
   {:input {:genome "data/genome.fa", :fastq "data/samples/{sample}.fastq"}},
   :samtools-sort
   {:input :bwa-map},
   :samtools-index {:input :samtools-sort},
   :bcftools-call {:input {:fa "data/genome.fa", :bam :samtools-sort, :bai :samtools-index}},
   :plot-quals {:input :bcftools-call}})

(rulegraph r)

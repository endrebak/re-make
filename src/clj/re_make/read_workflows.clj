(ns re-make.read-workflows
  (:import [java.io PushbackReader])
  (:require
   [re-make.state :as state]
   [clojure.set :as set]
   [com.stuartsierra.dependency :as dep]
   [clojure.java.io :as io])
  (:use [selmer.parser]))




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


;; (require '[clojure.set :as set])
;; (require '[com.stuartsierra.dependency :as dep])


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
     (let [kw-name# ~(keyword name)
           body# (assoc (handle-docs ~@body) :name kw-name#)]
       (swap! rules assoc kw-name# (handle-docs body#))
       (def ~(vary-meta name assoc :rule true) (handle-docs body#)))))

(defn rules->named-dependencies
  [rules]
  (for [[name v] rules
        :let [input (:input v)]]
    [name input]))

(defn rules->dependencies
  [rules]
  (let [m (rules->named-dependencies rules)]
    (for [[n v] m]
      (if (keyword? v)
        [n v]
        (for [x v] [n v])))))

(def r
  {:bwa-map
   {:input {:genome "data/genome.fa", :fastq "data/samples/{sample}.fastq"}},
   :samtools-sort
   {:input :bwa-map},
   :samtools-index {:input :samtools-sort},
   :bcftools-call {:input {:fa "data/genome.fa", :bam :samtools-sort, :bai :samtools-index}},
   :plot-quals {:input :bcftools-call}})

(defn rulegraph
  [rules]
  (let [r->d (rules->dependencies rules)]
    (reduce (dep/graph) ())))

;; (defn files->rules
;;   [rules])
(defn files->rules
  [])

(def v "baah")
(def vs [1 2 3 4])
(map #(vec [v %]) vs)
(for [x vs]
  [v x])
; [[v, x] for x in vs]


(defn property->rule
  [rules property]
  (partition 2 (flatten (for [[k v] rules :let [p (property v)] :when p]
                          (if (not (coll? (property v)))
                            [(property v) k]
                            (for [f (property v)]
                              [f k]))))))

(defn entries
  [c]
  (if (map? c)
    (vals c)
    c))

(defn rule->property
  [rules property]
  (partition 2
             (flatten
              (for [[k v] rules :let [p (property v)] :when p]
                (if (not (coll? p))
                  [k p]
                  (for [pp (entries p)]
                    [k pp]))))))



(defn add-dependency [g [a b]]
  (dep/depend g a b))


(defn filegraph
  [rules]
  (let [in (rule->property rules :in)
        out (for [[r f] (rule->property rules :out)] [f r])
        g (reduce add-dependency (dep/graph) in)]
    (reduce add-dependency g out)))



(defn rulegraph
  [rules]
  (let [in (rule->property rules :in)
        out (rule->property rules :out)
        f2r (into {} (for [[k v] out] [v k]))
        r2r (for [[k v] in] [k (get f2r v)])]
    (reduce add-dependency (dep/graph) r2r)))

(def rg (rulegraph @rules))

(defn out-targets
  [rules]
  (let [infiles (map first (property->rule rules :in))
        outfiles (map first (property->rule rules :out))]
    (set/difference (set outfiles) (set infiles))))



(defn rule-targets
  [rules] ;; TODO: opt argument out-targets
  (let [out-targets (out-targets rules)
        filegraph (:dependencies (filegraph rules))]
    (apply hash-map
           (flatten
            (for [f out-targets]
              (for [r (filegraph f)]
                [f r]))))))

(def rt (rule-targets @rules))
;; (defn get-xs
;;   [name xs]
;;   (name xs))

;; either have precomputed xs
;; or use cartesian product of xs

(defn expand
  [xs & ks]
  (let [ks (or (first ks) (keys xs))]
    (map #(zipmap ks %)
         (apply combo/cartesian-product (map #(xs %) ks)))))


(defn collect
  [to-collect x xs]
  (let [collected (xs to-collect)]
    (for [c collected]
      (assoc x to-collect c))))


(defn job-targets
  [rules & [precomputed]]
  (let [rule-targets (rule-targets rules)
        missing-targets (set/difference (set (vals rule-targets)) (set (keys precomputed)))
        postcomputed (map vec (combo/cartesian-product missing-targets (expand xs)))]
    (concat precomputed postcomputed)))


(defn jobgraph
  [rules precomputed]
  (let [job-targets (job-targets rules precomputed)
        rulegraph (rulegraph rules)]
    (for [[rulename xs] job-targets
          :let [deps ((:dependencies rulegraph) rulename)]]
      (do
        (println "----")
        (println rulename)
        (println deps)
        (println xs)))))

;; need to iterate over

;; either a rule has precomputed
;; xs or not
;; (defn jobgraph
;;   [rules targets]
;;   (let [rule-targets (rule-targets rules)
;;         missing (difference (set (vals rule-targets)) (set (keys targets)))
;;         ]))

;; (defn job-targets
;;   [rules xs]
;;   (let [rt (rule-targets rules)]
;;     (for [[file rulename] rt]
;;       (let [rule (rules rulename)
;;             ]))))


;; (defn jobgraph [rules xs]
;;   (let [fg (filegraph rules)
;;         ft (out-targets rules)]
;;     ))



;; (defn job->job
;;   [rules wildcards]
;;   (let [rule-targets (rule-targets rules)
;;         filegraph (filegraph rules)]))


;; (defn jobs
;;   [rule]
;;   (let []))
;; (defn all-combos
;;   [rule-wildcards wildcards]
;;   ())

;; (defn fill-rule
;;   [{:keys [in out external shell wildcards params]} all-wildcards]
;;   (let [infiles (fill-file)]))
  ;; (let [code ()])
  ;;       params (:params rule)
  ;;       external (:params external)]
  ;;   ()))

;; (defn format-string
;;   [])


;; (defn jobs
;;   [rule wildcards]
;;   (let [rule-wildcards (rule :wildcards)
;;         job-wildcards (select-keys wildcards rule-wildcards)]
;;     ())



;; what can happen:

;; waiting for job to finish, cannot go on

;; start independent jobs


; (map #(xs %) ks)
; (["A" "B" "C"] ["hg38" "hg19"])

;; (for [[f r] rt]
;;   (for [x (expand )])
;;   {:file f :rule r })


(def deps
  {:samtools-sort #{:bwa-map},
   :samtools-index #{:samtools-sort},
   :bcftools-call #{:samtools-index :samtools-sort},
   :plot-quals #{:bcftools-call}})

;; (defn rulerecur
;;   [rules target m]
;;   (for [r (rules target)]
;;     (recur rules r (assoc m r "whatevz"))))

;; (rulerecur deps :plot-quals {})



I'd like to iterate over a map representing a DAG and build a new map at the
same time. I have tried the below function, but it does not work. The second time it is called

(def deps
  {:samtools-sort #{:bwa-map},
   :samtools-index #{:samtools-sort},
   :bcftools-call #{:samtools-index :samtools-sort},
   :plot-quals #{:bcftools-call}})


(defn rulerecur
  [rules target v i]
  (let [targets (first (vals (select-keys rules target)))
        v (concat v (vec (for [t targets] [t i])))]
    (println (str "vec: " v))
    (println (str "new targets: " targets))
    (if-not (seq targets)
      v
      (recur rules targets v (inc i)))))

(rulerecur deps [:plot-quals] [] 5)

;; map: clojure.lang.LazySeq@a8e77100
;; new targets: (#{:bcftools-call})
;; map: clojure.lang.LazySeq@1
;; new targets:
;; ()

;;                                         ; #object[user$rulerecur$fn__32547 0x942e007 "user$rulerecur$fn__32547@942e007"]


;; (def m {})

;; (reduce m #(assoc {} % "ooo") [:a :b :c])

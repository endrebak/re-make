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


(defn rule->property
  [rules property]
  (partition 2
             (flatten
              (for [[k v] rules :let [p (property v)] :when p]
                (if (not (coll? p))
                  [k p]
                  (for [pp p]
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
        f2r (into {} (for [[k v] out] [v k]))
        r2r (for [[k v] in] [k (get f2r v)])]
    (reduce add-dependency (dep/graph) r2r)))


(defn targets
  [rules]
  (let [infiles (map first (property->rule rules :in))
        outfiles (map first (property->rule rules :out))]
    (set/difference (set outfiles) (set infiles))))

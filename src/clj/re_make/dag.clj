(ns re-make.dag
  (:require
   [com.stuartsierra.dependency :as dep]
   [re-make.state :as state]))



(defn eval-code! [code]
  (do
    (reset! state/rules {})
    (eval code)))


(defn my-flatten [x]
  (if (coll? x)
    (mapcat my-flatten x)
    [x]))


(defn depends-on
  [dependencies]
  (if (not (symbol? dependencies))
    (->> dependencies
         my-flatten
         (filter symbol?))
    [dependencies]))


(defn dag-pairs [rules]
   (for [[rule v] rules
         :let [dependencies (:input v)]
         :when (not (or (nil? rule) (nil? dependencies)))]
     (for [dependency (depends-on dependencies)]
       [rule dependency])))


(defn add-dependency [g [a b]]
  (dep/depend g a b))


(defn rulegraph [code]
  (-> code
      dag-pairs
      (fn [dag-pairs]
        (let [g (dep/graph)
              pairs (->> dag-pairs flatten (partition 2))]
          (reduce add-dependency g pairs)))))

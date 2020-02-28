(ns re-make.dag
  (require
   '[re-make/read-workflows :as rw]))

;; (fn)
;; when

;; need to find the dag
"create dag
check wildcards
"

(def dag (atom {}))

(defn find-dag [rules]

  )

;; (add-watch rw/rules :update-dag
;;            (fn [key atom old-state new-state]
;;              ))


;; (require '[com.stuartsierra.dependency :as dep])
;; Create a new dependency graph:

;; (def g1 (-> (dep/graph)
;;             (dep/depend :b :a)   ; "B depends on A"
;;             (dep/depend :c :b)   ; "C depends on B"
;;             (dep/depend :c :a)   ; "C depends on A"
;;             (dep/depend :d :c))) ; "D depends on C"
;; This creates a structure like the following:

;; need to fetch input and output from
;; need to create DAG to know which functions to run first
;; first create dag

;; :a
;; / |
;; :b  |
;; \ |
;; :c
;; |
;; :d
;; Ask questions of the graph:

;; (dep/transitive-dependencies g1 :d)
;; ;;=> #{:a :c :b}

;; (dep/depends? g1 :d :b)
;; ;;=> true
;; Get a topological sort:

;; (dep/topo-sort g1)
;; ;;=> (:a :b :c :d)

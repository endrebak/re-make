
(defrecord ChIP-seq [input output wildcards])

(defrecord Rule [name input output wildcards])

(defrecord download-bed [input output wildcards])

(defrule parent
  [Rule (= ?input input) (= ?child name) (= ?child-wildcards wildcards)]
  [Rule (= ?input output) (= ?parent name) (= ?parent-wildcards wildcards)]
  [:test (=child-)]
  =>
  (println (str "Output of " ?parent " matching input of " ?child " with " ?child-wildcards)))


(-> (mk-session 'clara.examples.validation)
    (insert (->Rule "ChIP-Seq" "f1" "f2" {:chromosome [1 2]})
            (->Rule "Functional analysis" "f2" "f3" {:chromosome [1 2 3 4]}))
    (fire-rules))


;; Have a lot of unordered rules
;; The relationships are given by DAG
;; can use dependency library to find bottom rules

;; insert rules and wildcards
;; when you have a DAG, can start reasoning about i

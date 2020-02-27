(ns re-make.read-workflows
  (import '[java.io PushbackReader])
  (require '[clojure.java.io :as io]))

;; first, read and eval the main workflow
;;     if more workflows are included, read those in order
;;
;; then find target, try to do that

;; (def read-workflow-file )
;; (eval )
;; (defn read-eval)
;; (eval (read-string "))

(def rules (atom {}))

;; (defmacro defrule [rulename rulebody]
;;   (list swap! assoc 'rules rulename rulebody))

(defn defrule [rulename rulebody] (swap! rules assoc rulename rulebody))

;; so need to pick out defrule
;; from the code

;; https://stackoverflow.com/a/24922859/992687
(defn read-all
  [file]
  (let [rdr (-> file io/file io/reader PushbackReader.)]
    (loop [forms []]
      (let [form (try (read rdr) (catch Exception e nil))]
        (if form
          (recur (conj forms form))
          forms)))))

;; will read the workflow
;; must parse out rules and put in rules-map
;; how should that be done

(defn read-workflow [f]
  (let [code (read-all f)
        to-include nil]
    code))

;; must parse file to select re-make specific stuff and other
;; (defn split-code)


(def wf "src/clj/re_make/workflow.clj")

(def code (read-workflow wf))

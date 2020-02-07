
(ns re-make.watch
  (:require
   [hawk.core :as hawk]))


(def main-workflow-file (atom []))


;; need to watch: the workflow file
;;   and the workflow subfiles they point to

;;

;; when either is changed, need to reeval whole workflow file and subworkflow-files
;; (include rules)
;; these then tell you what script files to watch for changes
;;    and what other produced files to watch for changes

(def workflow-files-to-watch (atom []))


;; (defn reeval-workflow-files-to-watch! []
;;   "read workflow definition and check for more files to watch
;; "
;;   )


;; (defn code-files-to-watch []
;;   "Always have main-workflow as one to watch"
;;   ())

;; (defn watch-main-file
;;   (hawk/watch! [{:paths ["."]
;;                  :handler (fn [ctx e]
;;                             (println "event: " e)
;;                             (println "context: " ctx)
;;                             ctx)}]))

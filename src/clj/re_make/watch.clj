(ns re-make.watch
  (:require
    [hawk.core :as hawk]
    [re-make.read-workflows :as rw]))


(def main-workflow-file (atom ""))

(def main-workflow-code (atom ()))

"find dag
find all rule describing files you want
find what files need to be created"

(def workflow-files-to-watch (atom []))

(defn watcher-read-main-file [ctx e]
  (when (= (:file e) @main-workflow-file)
    (rw/read-workflow @main-workflow-file)))

(def watcher
  (hawk/watch! [{:paths ["."]
                 :handler (fn [ctx e]
                            (println "event: " e)
                            (println "context: " ctx)
                            ctx)}]))

(defn stop [watcher]
  (hawk/stop! watcher))

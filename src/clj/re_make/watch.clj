
(ns re-make.watch
  (:require
   [hawk.core :as hawk]))


(def main-workflow-file (atom []))

"find dag
find all rule describing files you want
find what files need to be created"

(def workflow-files-to-watch (atom []))


(def watcher
  (hawk/watch! [{:paths ["."]
                 :handler (fn [ctx e]
                            (println "event: " e)
                            (println "context: " ctx)
                            ctx)}]))

(defn stop [watcher]
  (hawk/stop! watcher))

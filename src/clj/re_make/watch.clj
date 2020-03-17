(ns re-make.watch
  (:require
    [hawk.core :as hawk]
    [re-make.state :as state]
    ))



;; (defn watcher-read-main-file [ctx e]
;;   (when (= (:file e) @state/workflow-file)
;;     (rw/read-workflow @state/workflow-file)))

(def watcher
  (hawk/watch! [{:paths ["."]
                 :handler (fn [ctx e]
                            (println "event: " e)
                            (println "context: " ctx)
                            ctx)}]))

(defn stop [watcher]
  (hawk/stop! watcher))

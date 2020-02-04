(ns user
  "Userspace functions you can run by default in your local REPL."
  (:require
    [re-make.config :refer [env]]
    [clojure.spec.alpha :as s]
    [expound.alpha :as expound]
    [mount.core :as mount]
    ))

(alter-var-root #'s/*explain-out* (constantly expound/printer))

(add-tap (bound-fn* clojure.pprint/pprint))

;; (defn start
;;   "Starts application.
;;   You'll usually want to run this on startup."
;;   []
;;   (mount/start-without #'re-make.core/repl-server))

;; (defn stop
;;   "Stops application."
;;   []
;;   (mount/stop-except #'re-make.core/repl-server))

;; (defn restart
;;   "Restarts application."
;;   []
;;   (stop)
;;   (start))

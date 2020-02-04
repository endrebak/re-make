(ns re-make.env
  (:require
    [selmer.parser :as parser]
    [clojure.tools.logging :as log]
    [re-make.dev-middleware :refer [wrap-dev]]))

(def defaults
  {:init
   (fn []
     (parser/cache-off!)
     (log/info "\n-=[re-make started successfully using the development profile]=-"))
   :stop
   (fn []
     (log/info "\n-=[re-make has shut down successfully]=-"))
   :middleware wrap-dev})

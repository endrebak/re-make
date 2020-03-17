(ns re-make.read-workflows
  (:import [java.io PushbackReader])
  (:require
   [re-make.state :as state]
   [com.stuartsierra.dependency :as dep]
   [clojure.java.io :as io]))


(defn defrule! [rulename rulebody]
  (let [rulebody (assoc rulebody :name rulename)]
    (swap! state/rules assoc rulename rulebody)))


(defn defrule [rulename rulebody]
  (defrule! rulename rulebody))


;; https://stackoverflow.com/a/24922859/992687
(defn read-all
  [file]
  (let [rdr (-> file io/file io/reader PushbackReader.)]
    (loop [forms []]
      (let [form (try (read rdr) (catch Exception e nil))]
        (if form
          (recur (conj forms form))
          forms)))))


(defn read-workflow [f]
  (let [code (read-all f)
        to-include nil]
    code))

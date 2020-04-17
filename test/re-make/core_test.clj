(ns re-make.core-test
  (:use midje.sweet)               ;; <<==
  (:require [re-make.read-workflows :as rw]))

(facts "about `split`"
       (str/split "a/b/c" #"/") => ["a" "b" "c"]
       (str/split "" #"irrelvant") => [""]
       (str/split "no regexp matches" #"a+\s+[ab]") => ["no regexp matches"])

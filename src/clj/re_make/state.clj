(ns re-make.state)

(def workflow-file (atom "src/clj/re_make/chip_seq.clj"))

(def rules (atom {}))

(def code (atom []))

(def dag (atom {}))

(def rulegraph (atom {}))


;; btw this can be done without atoms if you have an event flow design: you make
;; a function that takes an immutable map of data, does some updates, then
;; passes it to your consumers (maybe a list of functions defined by third
;; parties) and lets them attempt updates (which you then validate and accept or
;; reject)

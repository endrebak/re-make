(ns re-make.state)

(def workflow-file (atom "src/clj/re_make/chip_seq.clj"))

(def rules (atom {}))

(def code (atom []))

(def dag (atom {}))

(def rulegraph (atom {}))

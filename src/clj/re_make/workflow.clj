(ns re-make.workflow)

(def all
  [{:from create-files
    :wildcards {:f1 '(1 2 3 4) :f2 '(1)}}])


(defrule create-files
  {:output "/mnt/work/endrebak/re-make/{file}.tsv"
   ;; :vars {:file }
   :run
   "
#!/usr/bin/env python

import hashlib

f = wildcards.file
seed = int(hashlib.sha1(s.encode()).hexdigest(), 16) % (10 ** 8)

import numpy as np
np.random.seed(seed))

import pyranges as pr
gr = pr.random()
gr.to_csv({output[0]}, sep='\t')
"})


(defrule nearest
  {:input {create-files }
   :output "/mnt/work/endrebak/re-make/nearest_{f1}_{f2}.tsv"
   :run
   "
#!/usr/bin/env python
import pandas as pd
import pyranges as pr
gr = pr.PyRanges(input[0])
gr2 = pr.PyRanges(input[1])
gr = gr.nearest(gr2)
gr.to_csv({output[0]}, sep='\t')
"})

(defrule stats
  {:input nearest
   :output "/mnt/work/endrebak/re-make/stats.tsv"
   })
;; collect desired results
;; collect rules
;; find what output they produce
;; start watching for changes
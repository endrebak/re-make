(def wildcards ())

(defrule 'create-bed
  {:run
   "#!/usr/bin/env python

   import pyranges as pr
   gr = pr.random()
   gr.to_csv({output[0]}, sep='\t')"})

(defrule 'liftover
  {:input 'create-bed
   :run
   "cp {input} {output}"})

(defrule 'lengths
  {:input 'liftover
   :run
   "#!/usr/bin/env python

   import pyranges as pr
   gr = pr.read_bed(input)
   gr = gr[gr.lengths() >= 100]
   gr.to_csv({output})
   "})


(defrule 'epic2
  {:input {:bed '[lengths create-bed]}
   :wildcards '[genome]
   :input-switch #(case (:genome wildcards)
                    "hg38" 'lengths
                    "hg19" 'create-bed)
   :run
   "cp {input} {output}"})


;; how postpone the call?

;; - first eval all def and defn

;; - then collect defrules and add them to rules
;;

;; (defrule :all
;;   [{:from 'create-files
;;     :wildcards {:f1 '(1 2 3 4) :f2 '(1)}}])


;; (defrule :create-files
;;   {:output "/mnt/work/endrebak/re-make/{file}.tsv"
;;    ;; :vars {:file }
;;    :run
;;    "
;; #!/usr/bin/env python

;; import hashlib

;; f = wildcards.file
;; seed = int(hashlib.sha1(s.encode()).hexdigest(), 16) % (10 ** 8)

;; import numpy as np
;; np.random.seed(seed))

;; import pyranges as pr
;; gr = pr.random()
;; gr.to_csv({output[0]}, sep='\t')
;; "})


;; (defrule :nearest
;;   {:input 'create-files
;;    :output "/mnt/work/endrebak/re-make/nearest_{f1}_{f2}.tsv"
;;    :run
;;    "
;; #!/usr/bin/env python
;; import pandas as pd
;; import pyranges as pr
;; gr = pr.PyRanges(input[0])
;; gr2 = pr.PyRanges(input[1])
;; gr = gr.nearest(gr2)
;; gr.to_csv({output[0]}, sep='\t')
;; "})

;; (defrule stats
;;   {:input nearest
;;    :output "/mnt/work/endrebak/re-make/stats.tsv"
;;    })
;; collect desired results
;; collect rules
;; find what output they produce
;; start watching for changes

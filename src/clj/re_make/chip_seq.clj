;; (ns chip-seq)

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

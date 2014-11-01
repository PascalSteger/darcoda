(TeX-add-style-hook "Steger+_2014_Gravlite"
 (lambda ()
    (LaTeX-add-bibliographies
     "thesis")
    (TeX-add-symbols
     '("TODO" 1)
     '("changefont" 3)
     '("bcite" 1)
     "widebar"
     "ltsima"
     "simlt"
     "gtsima"
     "simgt"
     "siglos"
     "siglosi"
     "vztwo"
     "vztwoi"
     "rhodisc"
     "rhodmext"
     "rhodm"
     "rhoeff"
     "nuobs"
     "tot"
     "pc"
     "kpc")
    (TeX-run-style-hooks
     "amsmath"
     "color"
     "ulem"
     "graphicx"
     "myaasmacros"
     ""
     "latex2e"
     "mn2e10"
     "mn2e"
     "usenatbib"
     "useAMS"
     "intro"
     "method"
     "results"
     "conclusion"
     "acks"
     "appendix")))


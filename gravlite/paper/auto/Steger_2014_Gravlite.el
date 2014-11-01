(TeX-add-style-hook "Steger_2014_Gravlite"
 (lambda ()
    (LaTeX-add-bibliographies)
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
     "kpc"
     "GravLite"
     "MultiNest")
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
     "letter"
     "usenatbib"
     "useAMS"
     "intro"
     "method"
     "results"
     "conclusion"
     "acks"
     "appendix")))


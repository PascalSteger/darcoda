(TeX-add-style-hook "main"
 (lambda ()
    (LaTeX-add-bibliographies)
    (TeX-add-symbols
     '("TODO" 1)
     '("changefont" 3)
     '("bcite" 1)
     "ltsima"
     "simlt"
     "gtsima"
     "simgt"
     "siglos"
     "siglosi"
     "gprior"
     "cprior"
     "bprior"
     "lbprior"
     "vztwo"
     "vztwoi"
     "rhodisc"
     "rhodmext"
     "rhodm"
     "rhoeff"
     "nuobs"
     "tot"
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
     "result"
     "conclusion"
     "acks")))


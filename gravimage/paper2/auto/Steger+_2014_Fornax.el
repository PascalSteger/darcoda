(TeX-add-style-hook
 "Steger+_2014_Fornax"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("mn2e" "useAMS" "usenatbib")))
   (TeX-run-style-hooks
    "latex2e"
    "intro"
    "method"
    "results"
    "conclusion"
    "acks"
    "appendix"
    "mn2e"
    "mn2e10"
    "myaasmacros"
    "graphicx"
    "ulem"
    "color"
    "amsmath")
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
    "gprior"
    "cprior"
    "bprior"
    "lbprior"
    "GravImage")
   (LaTeX-add-bibliographies
    "thesis")))


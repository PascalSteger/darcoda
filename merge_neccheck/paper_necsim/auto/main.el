(TeX-add-style-hook "main"
 (lambda ()
    (LaTeX-add-bibliographies)
    (LaTeX-add-labels
     "firstpage"
     "sec:intro"
     "sec:meth"
     "fig:vis"
     "tab:simprop"
     "tab:haloprop"
     "sec:res"
     "fig:dm_prof"
     "fig:prof_dmonly"
     "fig:star_prof"
     "fig:vis_proj_dm_sod_ahf"
     "fig:ov_1"
     "fig:ov_4"
     "fig:hist_shrsph"
     "fig:star_proj_1"
     "sec:sum"
     "sec:ack"
     "sec:app"
     "fig:corr_r_core"
     "lastpage")
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
     "kpc"
     "Msun")
    (TeX-run-style-hooks
     "amsmath"
     "color"
     "graphicx"
     "myaasmacros"
     ""
     "latex2e"
     "mn2e10"
     "mn2e"
     "usenatbib"
     "useAMS")))


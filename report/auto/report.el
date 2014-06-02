(TeX-add-style-hook "report"
 (lambda ()
    (LaTeX-add-labels
     "sec-1"
     "sec-2"
     "sec-3"
     "sec-3-1"
     "sec-3-2"
     "sec-3-3"
     "sec-3-4"
     "sec-3-5"
     "sec-4"
     "sec-5"
     "sec-6")
    (TeX-run-style-hooks
     "latexsym"
     "soul"
     "fancyheadings"
     "moreverb"
     "babel"
     "english"
     "swedish"
     "mathtools"
     "parskip"
     "parfill"
     "hyperref"
     "amssymb"
     "wasysym"
     "marvosym"
     "textcomp"
     "amsmath"
     "ulem"
     "normalem"
     "rotating"
     "wrapfig"
     "float"
     "longtable"
     "graphicx"
     "fixltx2e"
     ""
     "fontenc"
     "T1"
     "inputenc"
     "utf8"
     "latex2e"
     "art11"
     "article"
     "11pt")))


# Package

version       = "0.5.1"
author        = "G. Bareigts"
description   = "Renormalization of colloidal charges of polydipserse dispersions using the Poisson-Boltzmann equation"
license       = "MIT"

# Dependencies

requires "nim >= 1.6.12"
requires "https://github.com/mratsim/constantine#f57d071f1192a4039979a3baf6c835b89841bcfa"

# Build

srcDir = "src"
bin = @["polypbren"]
binDir = "bin"

# Doc

task gendoc, " Generate html documentation.":
  mkDir("docs/html")
  exec "nim rst2html --index:on --out=docs/html/index.html docs/index.rst"
  exec "nim rst2html --index:on --out:docs/html/usage.html docs/usage.rst"
  exec "nim rst2html --index:on --out:docs/html/example.html docs/example.rst"
  exec "nim --index=on doc2 --out=docs/html/pbsolv.html src/polypbrenpkg/pbsolv.nim"
  #exec "nim --index=on doc2 --out=docs/html/pbren.html src/polypbrenpkg/pbren.nim"
  #exec "nim --index=on doc2 --out=docs/html/distribs.html src/polypbrenpkg/distrib.nim"
  #exec "nim --index=on doc2 --out=docs/html/params.html src/polypbrenpkg/params.nim"
  #exec "nim buildIndex --out=docs/html/theindex.html docs/html"

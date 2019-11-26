#********************** Licensing information ****************************
# Copyright (c) 2018 Guillaume Bareigts
# MIT License, see LICENSE
#*************************************************************************
# Main program

import os, math, parseopt, strutils, cpuinfo
import polypbrenpkg / [pbren, params, distrib]

proc gitRevision(): string {.compileTime.} =
  let ret = gorgeEx("git rev-parse --short HEAD")
  if ret.exitCode == 0:
    return "-git" & ret.output
  else:
    return ""

proc getHostName(): string {.compileTime, used.} =
  let ret = when defined(windows): gorgeEx("hostname") else: gorgeEx("uname -n")
  if ret.exitCode == 0:
    return ret.output
  else:
    return ""


const
  version = "0.4.2"
  revision = gitRevision()
  buildInfo = when defined(printBuildInfo) :
                "Compiled on " & getHostName() & ", " & CompileDate & ", " & CompileTime & " with Nim " & NimVersion 
              else:
                ""

proc writeVersion() =
  echo paramStr(0),", version ", version, revision
  echo buildInfo

var 
  fileParas* = "polypbren.cfg" 
  fileDistrib* = "distrib.in"
  maxThreads* = 1

proc writeHelp() =
  echo """Usage: polypbren [OPTIONS]

Options:
  -h, --help                 Print help message
  -v, --verbosity:N          Set the degree of verbosity. Levels are 0,1,2 (default: 0).
  -d, --distribution:FILE    Change the distribution file to FILE (default: distrib.in)
  -p, --parameters:FILE      Change the distribution file to FILE (default: polypbren.cfg)
  --mth, --max-threads:N     Set the maximum of threads to use. 0 means it is auto-detected. (default: 1)
  --version                  Print version

"""
proc parseCmdLine() =
  for kind, key,val in getopt():
    case kind
    of cmdLongOption, cmdShortOption:
      case key
      of "h","help":
        writeHelp()
        quit()
      of "d","distribution":
        fileDistrib = val
      of "p","parameters":
        fileParas = val
      of "v", "verbosity":
        verb = parseInt(val)
      of "version":
        writeVersion()
        quit()
      of "mth", "max-threads":
        maxThreads = parseInt(val)
      else:
        echo "Bad command line argument"
        writeHelp()
        quit(2)
    else: discard

proc main() =
  parseCmdLine()

  if verb >= 1:
    writeVersion()
    echo()

  maxThreads = if maxThreads == 0: countProcessors() else: maxThreads

  loadParams(fileParas)

  if verb >= 1:
    echo "Ion types:"
    for i in 0..ionDensv.high:
      echo ionChv[i],"  ",ionDensv[i]
    echo()

  let dist = readDistrib(fileDistrib)
  if verb >= 1:
    echo "Distribution:"
    for np,r,s in dist.items:
      echo np,' ',r,' ',s
    echo()

  let pbres = doCalculations(dist, maxThreads)

  proc formSci(x: float): string {.inline.} = 
    formatFloat(x, precision=7, format=ffScientific)

  var fout = open("distrib.out",fmWrite)
  fout.writeLine(dist.len)
  for f,fam in dist:
    fout.writeLine(fam.np,' ',fam.r,' ',pbres.sigmaEffv[f].formSci)
  fout.close
  fout = open("peffs.dat",fmWrite)
  let (volmoy, _) = calcVolmNptot(dist)
  let press = kBpnm3 * temperature * (pbres.volfrac/(volmoy) + pbres.rhoEdge - ionDensv.sum)
  fout.writeLine(pbres.volfrac.formSci, ' ', pbres.kappaEff.formSci, ' ', press.formSci)
  fout.close

  for f,fam in dist:
    fout = open("phi" & $f & ".dat",fmWrite)
    for i,phi in pbres.finalPhiv[f]:
      let rr = pbres.finalMesh[f][i] #+ fam.r+i.float*Sbin 
      fout.writeLine(rr, ' ',phi, ' ', calcIonDensity(phi, lambdaB, ionChv, ionDensv))
    fout.close

  if solveCondition==PHCond:
    var fout = open("distrib-bare.out",fmWrite)
    fout.writeLine("  ", dist.len)
    for f,fam in dist:
      fout.writeLine(fam.np,' ',fam.r,' ',pbres.finalSig[f].formSci)
    fout.close
    if verb >= 1:
      echo "\nCharges for PH ", pHInput, " :"
      for f,fam in dist:
        echo(fam.r,' ',pbres.finalSig[f].formSci)

try:
  main()
except:
  stderr.writeLine("Error: ", getCurrentExceptionMsg())
  quit(QuitFailure)

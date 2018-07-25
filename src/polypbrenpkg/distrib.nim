#********************** Licensing information ****************************
# Copyright (c) 2018 Guillaume Bareigts
# MIT License, see LICENSE
#*************************************************************************
## This module helps to import *distrib* files.


import strutils

type
  FamComp* = tuple
    ## **np** : number of representants, 
    ## **r** radius, **ch**:  charge
    np: int       
    r, ch: float
  
  Distrib* = seq[FamComp]

proc readDistrib*(distFile: string): Distrib =
  ## Reads a distribution file of name **distFile**
  ## and returns a *Distrib* object.
  var fdist = open(distFile,fmRead)
  let nfam = fdist.readLine.strip.parseInt
  result = newSeqofCap[FamComp](nfam)
  var dump: seq[string]
  for i in 0..<nfam:
    try:
      dump = fdist.readLine.splitWhiteSpace
    except IOError:
      raise newException(IOError, "Problem in " & distFile & ", line " & $(i+2))
    if dump.len == 3:
      result.add(( parseInt(dump[0]), parseFloat(dump[1]),parseFloat(dump[2])))
    else:
      raise newException(IOError, "Problem cat 2 in " & distFile & ", line " & $(i+2))
  fdist.close

proc readDistrib2*(distFile: string): tuple[np: seq[int], r, ch: seq[float]] =
  ## Reads a distribution file of name **distFile**
  ## and returns the number, radius and charge sequences.
  var fdist = open(distFile,fmRead)
  let nfam = fdist.readLine.strip.parseInt
  result.np = newSeqofCap[int](nfam)
  result.r = newSeqofCap[float](nfam)
  result.ch = newSeqofCap[float](nfam)
  var dump: seq[string]
  for i in 0..<nfam:
    try:
      dump = fdist.readLine.split
    except IOError:
      raise newException(IOError, "Problem in " & distFile & ", line " & $(i+2))
    if dump.len == 3:
      #result.add(Distrib(np: parseInt(dump[0]),r: parseFloat(dump[1]),ch: parseFloat(dump[2])))
      result.np.add(parseInt(dump[0]))
      result.r.add(parseFloat(dump[1]))
      result.ch.add(parseFloat(dump[2]))
    else:
      raise newException(IOError, "Problem cat 2 in " & distFile & ", line " & $(i+2))
  fdist.close

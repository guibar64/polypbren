#********************** Licensing information ****************************
# Copyright (c) 2018 Guillaume Bareigts
# MIT License, see LICENSE
#*************************************************************************
# Input parameters and loaders

import sequtils, strutils, parseutils, math
import pbsolv, distrib

type
  ModelType* = enum
    InvalidModel
    Jellium
    Cellmod
    Jellbid
  ConditionType* = enum
    ChargeCond, PHCond

proc parseModel(s: string): ModelType =
  case s
  of "jellium":
    result = Jellium
  of "cellmod":
    result = Cellmod
  of "jellium2":
    result = Jellbid
  else:
    result = InvalidModel

proc parseCondition(s: string): ConditionType =
  case s
  of "charge":
    result = ChargeCond
  of "pH":
    result = PHCond
  else:
    result = ChargeCond
    
const
  mM* = 0.0006022140857 # in nm^-3
  kBpnm3* = 1.0e27 * 1.38064852e-23
var
  verb*: range[0..2] = 0
  model* = Jellium
  volFracIn* = 0.08
  maxSteps* = 20
  maxIters* = 50
  maxInIters* = 100000
  tolPhi0* = 0.00001
  tolZeff* = 0.001
  tolCellm* = 0.0001
  tolVolFrac* = 0.001
  sbin* = 0.05
  phiD* = 0.5
  lambdaB* = 0.7105
  slength* = 50.0
  eqKind* = RadialPBE
  tolInLoop* = 0.00001
  solveCondition* = ChargeCond
  pHInput* = 9.0
  pKa* = 7.7
  siteDens* = 5.5 # nm^-2
  sternL* = 0.1
  temperature* = 300.0 # K
  max_potential_zeff_jellium* = 0.5
  factCalcZeffJellium* = 3.5   # Unused
  factPhi0ChargeCond* = 1.2

var 
  ionChv*: seq[float] = @[-1.0, 1.0]
  ionDensv*: seq[float] = @[5*mM, 5*mM]

const
  ValueChars = IdentChars + {'.','-'}

iterator configKeyVals*(fileName: string): tuple[key, val: string] =
  template raiseSyntaxError(msg: string) =
    raise newException(KeyError, "Syntax error at line " & $nl & ": " & msg)
  var key, val = newStringOfCap(64)
  var nl = 0
  for line in lines(fileName):
    inc(nl)
    var pos = skipWhiteSpace(line, 0)
    let p1 = parseWhile(line, key, IdentChars, pos)
    if p1 == 0 or pos == line.len or line[pos] == '#': continue # comment
    pos += p1
    pos += skipWhiteSpace(line, pos)
    if pos >= line.len or line[pos] == '#':
      yield (key, "")
    elif line[pos] == '=':
      inc(pos) ; pos += skipWhiteSpace(line, pos)
      if pos >= line.len:
        yield (key, "")
      else:
        setLen(val, 0)
        case line[pos]
        of ValueChars:
          discard parseWhile(line, val, ValueChars, start = pos)
          yield (key, val)
        of '"':
          discard parseUntil(line, val, '"', start = pos+1)
          yield (key, val)
        of '[':
          var ok = false
          for i in pos+1..<line.len:
            if line[i] == ']':
              ok = true
              break
            elif line[i] notin WhiteSpace:
              # Filter whitespace for easier spliting + parsefloat
              val.add(line[i])
          if not ok:
            raiseSyntaxError("'[' must be closed by a ']'")
          else:
            yield (key, val)
        else:
          raiseSyntaxError("Unexpected character: " & line[pos])
    else:
      raiseSyntaxError("")
    
  
proc loadParams*(fileParas: string) =
  try:
    for key, val in configKeyVals(fileParas):
      case normalize(key)
      of "model": model = parseModel(val)
      of "condition": solveCondition = parseCondition(val)
      of "bjerrumlength": lambdaB = parseFloat(val)
      of "volumefraction": volFracIn = parseFloat(val)
      of "tolvolumefraction": tolVolFrac = parseFloat(val)
      of "maxofsteps": maxSteps = parseInt(val)
      of "maxofiterations": maxIters = parseInt(val)
      of "maxofinternaliterations": maxInIters = parseInt(val)
      of "tolinitvalue": tolPhi0 = parseFloat(val)
      of "tolzeff": tolZeff = parseFloat(val)
      of "tolcellcondition": tolCellm = parseFloat(val)
      of "spacebin": sbin = parseFloat(val)
      of "celllength": slength = parseFloat(val)
      of "externpotential": phiD = parseFloat(val)
      of "ioncharges": ionChv =  val.split(',').map(parseFloat)
      of "iondensities": ionDensv = val.split(',').mapIt(parseFloat(it)*mM)
      of "verbosity": verb = parseInt(val)
      of "planar": eqKind = PlanePBE
      of "tolinnerloop": tolInLoop = parseFloat(val)
      of "phinput": pHInput = parseFloat(val)
      of "pka": pKa = parseFloat(val)
      of "sitedensity": siteDens = parseFloat(val)
      of "sternlength": sternL = parseFloat(val)
      of "temperature": temperature = parseFloat(val)
      of "maxpotentialzeffjellium": max_potential_zeff_jellium = parseFloat(val)
      of "factcalczeffjellium": factCalcZeffJellium = parseFloat(val)
      of "factphi0chargecond": factPhi0ChargeCond = parseFloat(val)
      else:
        quit "Error: " & fileParas & ": Unkown parameter : " & key, 1
  except KeyError, ValueError:
    stderr.writeLine "Error in ", fileParas, ": ", getCurrentExceptionMsg()
    quit(QuitFailure)
  except OSError, IOError:
    stderr.writeLine "Error in ", fileParas, ": ", getCurrentExceptionMsg()
    quit(QuitFailure)
  if ionDensv.len != ionChv.len:
    quit "Error: Mismatch between ion densities and ion charges",1


proc calcVolmNptot*(dist: Distrib): (float, int) =
    var volm = 0.0
    var nptot = 0
    for np,r,s in dist.items:
      volm += np.float*pow(r,3)
      nptot += np
    volm *= (4*PI/3.0)/nptot.float
    return (volm, nptot)


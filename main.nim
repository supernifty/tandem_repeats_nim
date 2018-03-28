
import logging
import parseopt
import sets
import strutils
import tables
import times

proc log(s: string) =
  var
    startTimeInfo = getLocalTime(fromSeconds(epochTime()))
    startTimeStr = startTimeInfo.format("HH:mm:ss")
  writeLine(stderr, "$1 $2".format(startTimeStr, s))
  flushFile(stderr)

# find kmer at position
# we check each possible kmer and see if it is at this position
# this isn't efficient
proc find_kmer_at_pos(data: string, pos: int, repeat: int): string =
  var 
    kmer = data[pos .. pos + repeat - 1]
    ok = false
  if kmer.contains('N'):
    return ""
  if repeat > 2 and repeat mod 2 == 0 and kmer[0 .. repeat / 2 - 1] == kmer[repeat / 2 .. ^1]: # abab
    return ""
  if repeat > 1:
    for x in 1 .. kmer.len - 1: # aaa
      if kmer[x] != kmer[0]:
        ok = true
        break
    if not ok:
      return ""
  return kmer

# find the longest tandem run starting from pos
proc find_run(data: string, kmer: string, start: int, repeat: int): int =
  var
    runlen = 0
    pos = start
  while data[pos .. pos + repeat - 1] == kmer:
    runlen += repeat
    pos += repeat
    if pos >= data.len:
      return runlen
  return runlen

# find tandem repeats of a kmer in chrom
# chrom: chromosome name
# kmers: each kmer and its location in the chromosome
# repeat: repeat length 
# minlen: minimum repeat length
# max_pos: last position in kmers
proc find_tandems(chrom: string, chrom_data: string, repeat: int, minlen: int) =
  var
    run_len: int
    run_pos = 0
    max_run = 0
  log("checking $1...".format(chrom))
  for pos in 0 .. chrom_data.len - repeat:
    if pos mod 10000000 == 0:
      log("checked $1:$2...".format(chrom, pos))
    if pos < run_pos:
      continue
    var kmer = find_kmer_at_pos(chrom_data, pos, repeat)
    if kmer != "":
      run_len = find_run(chrom_data, kmer, pos, repeat)
      if run_len > max_run:
        max_run = run_len
      if run_len >= minlen:
        writeLine(stdout, "$1\t$2\t$3\trepeat=$4;length=$5".format(chrom, pos, pos + run_len, kmer, run_len))
        run_pos = pos + run_len
  log("checking $1: done. max run: $2".format(chrom, max_run))
    
# main function - given an input fasta file, find tandem repeats
# fh: input file handle
# repeat: repeat length
# minlen: min length repeat to find
proc find(fh: File, repeat: int, minlen: int) =
  var 
    line: string
    pos = 0
    last_log = 0
    chrom: string
    chrom_data: seq[string]

  # build a list of kmers and their positions in each chromosome
  log("looking for repeat $1 of minlen $2...".format(repeat, minlen))
  while readLine(stdin, line):
    if line.startsWith(">"):
      if chrom_data.len > 0:
        find_tandems(chrom, join(chrom_data, ""), repeat, minlen)
      chrom = line[1..^1]
      pos = 0
      chrom_data = @[]
      log("processing $1" % chrom)
    else:
      line = line.toUpperAscii()
      pos += line.len
      chrom_data.add(line)
      if pos - last_log >= 10000000:
        log("read $1:$2".format(chrom, pos))
        last_log = pos

  if chrom_data.len > 0:
    find_tandems(chrom, join(chrom_data, ""), repeat, minlen)
  log("done")

proc main() =
  var
    repeat = 2
    minlen = 6

  for kind, key, val in getopt():
    case kind
    of cmdArgument:
      discard
    of cmdShortOption, cmdLongOption:
      case key
      of "repeat", "r": repeat = parseInt(val)
      of "min", "m": minlen = parseInt(val)
      else: discard
    of cmdEnd: 
      discard
  find(stdin, repeat, minlen)

when isMainModule:
  main()

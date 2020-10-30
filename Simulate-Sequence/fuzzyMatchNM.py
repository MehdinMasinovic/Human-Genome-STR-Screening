import regex

def fuzzyMatch(seq, sequence):
  
  pattern = "((" + seq + "){6,}){s<1}"

  positionsnm = [(m.start(), m.end(), m.group()) for m in regex.finditer("{}".format(pattern), sequence)]

  return positionsnm
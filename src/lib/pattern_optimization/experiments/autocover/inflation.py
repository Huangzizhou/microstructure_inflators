import paths, os, re, subprocess

def commonCommandArray(pat, dim, constraints):
    args = [paths.inflator(dim), paths.pattern(pat, dim), "-I"]
    for c in constraints: args += ["-C", c]
    return args

def getParameterTypes(pat, dim = 3, constraints=[]):
    cmd = commonCommandArray(pat, dim, constraints)
    inflOut = subprocess.check_output(cmd)
    types = []
    for p in inflOut.strip().split('\n'):
        m = re.match("^Param [^:]+: ([^,]+), [-0-9.]+$", p)
        types += [m.group(1)]
    return types

def inflate(pat, params, outPath, args = ['-S1'], dim = 3, constraints=[]):
    cmd = commonCommandArray(pat, dim, constraints)
    cmd += ["-p", " ".join(map(str, params)), '-o', outPath]
    cmd += args
    subprocess.check_output(cmd)

def isPrintable(pat, params, dim = 3, constraints=[]):
    cmd = commonCommandArray(pat, dim, constraints)
    cmd += ["-P", "-p", " ".join(map(str, params))]
    inflOut = subprocess.check_output(cmd)
    printableString = None
    for l in inflOut.strip().split('\n'):
        m = re.match("^Printable:\s+(\S+)", l)
        if (m): printableString = m.group(1)
    if   (printableString == "1"): return True
    elif (printableString == "0"): return False
    raise Exception("Invalid response from inflator to printability query")

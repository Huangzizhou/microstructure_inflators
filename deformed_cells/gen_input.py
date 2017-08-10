import numpy as np, math, sys;

if len(sys.argv) != 2:
    raise Exception("Usage: gen_input.py nsubdiv");
nsubdiv=int(sys.argv[1]);
thetas = np.linspace(0, np.pi, nsubdiv);

spacing = math.exp(math.log(4) / (nsubdiv / 2));
lambdas = spacing ** np.arange(-(nsubdiv / 2),(nsubdiv / 2 + 1));

for t in thetas:
    for l in lambdas:
        print t, l

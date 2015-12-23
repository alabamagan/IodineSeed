#!/usr/bin/python


import numpy as np
import IodineSeed

def main():
    seedexp = IodineSeed.IodineSeed([0,0,0], 10., [1,0,0], 1)
    # seedexp._energy(10, 90)
    seedexp.SetSize(10,10,10)
    seedexp.SetSpacing(0.2,0.2,0.2)
    out = seedexp.Update()
    print out
    np.save("out", out)

if __name__ == '__main__':
    main()
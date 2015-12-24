#!/usr/bin/python


import numpy as np
import IodineSeed

def main():
    seedexp = IodineSeed.IodineSeed([0,0,0], 4., [0,-1,0], 1)
    # seedexp._energy(5, 100)
    # seedexp._energy(5, 80)
    seedexp.SetSize(10,10,10)
    seedexp.SetSpacing(0.2,0.2,0.2)
    out = seedexp.Update2D()
    print out
    np.save("out", out)
    # print seedexp._mg(2,10)

if __name__ == '__main__':
    main()
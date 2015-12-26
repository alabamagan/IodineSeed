#!/usr/bin/python


import numpy as np
import IodineSeed

def main():
    seedexp = IodineSeed.IodineSeed([0,0,0], 4., [0,-1,0], 1)
    # seedexp._energy(5, 100)
    # seedexp._energy(5, 80)
    seedexp.SetSize(2,1,20)
    seedexp.SetSpacing(0.01,0.01,0.1)
    out = seedexp.Update2D()
    np.save("out", out)
    # print seedexp._mg(1.4,10)

if __name__ == '__main__':
    main()
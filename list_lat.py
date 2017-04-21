#!/usr/bin/env python

from find_elite import elite_int_pars

if __name__ == "__main__":

  lats = range(40,70,5)

  for lat in lats:

    print "%d: %s"%(lat, elite_int_pars('RAYMO/TEST%d'%lat,81))

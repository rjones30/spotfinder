#!/usr/bin/env python3
#
# cobrems_master.py - Celery tasks for parallel computing of coherent
#                     bremsstrahlung spectra, used by the wsgi script
#                     spotfinder.py for real-time computing of coherent
#                     bremsstrahlung photon beam spectra.
#
# author: richard.t.jones at uconn.edu
# version: july 14, 2022

import ROOT
import pickle
import base64
import random
import time

from cobrems_worker import fill_histo

def sum_histo(nprocs=100000, nrand=100000, maxbacklog=900):
   for b in range(0, nprocs, maxbacklog):
      bmax = min(b + maxbacklog, nprocs)
      procs = [fill_histo.delay(i, nrand) for i in range(b, bmax)]
      hsum = 0
      while len(procs) > 0:
         pending = []
         for proc in procs:
            if proc.ready():
               msg = proc.get()
               h = pickle.loads(base64.b64decode(msg))
               if hsum:
                  hsum.Add(h)
               else:
                  hsum = h
            else:
               pending.append(proc)
         print(f"{len(pending)} still pending in batch {b}")
         procs = pending
         time.sleep(1)
   return hsum

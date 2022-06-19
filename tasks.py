#!/usr/bin/env python3

from celery import Celery

app = Celery("tasks", backend="rpc://", broker="amqp://guest@localhost//")

@app.task
def add(i, j):
   return i + j

def parallel_add(nreps, batchsize=900):
   s = 0
   for b in range(0, nreps, batchsize):
      bmax = min(nreps, b + batchsize)
      procs = [add.delay(i,i) for i in range(b, bmax)]
      while len(procs) > 0:
         pending = []
         for proc in procs:
            if proc.ready():
               s += proc.get()
            else:
               pending.append(proc)
         print(f"batch {b} completed with {len(pending)} still pending")
         procs = pending
   print(f"got answer {s}")

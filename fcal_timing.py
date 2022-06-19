#!/usr/bin/env python3

import sys
import hddm_s
import ROOT
import multiprocessing
import tqdm

rf_period = 4.008 # ns
clight = 30.0 # cm/ns
ztarget = 65 # cm
fcal_middle_row = 29
fcal_middle_col = 29

def usage():
   print("usage: fcal_timing.py <infile.hddm> [...]")
   sys.exit(1)

def hist1D(name, title, nbins, x0, x1, xtitle, ytitle, hist):
   hist[name] = ROOT.TH1D(name, title, nbins, x0, x1)
   hist[name].GetXaxis().SetTitle(xtitle)
   hist[name].GetYaxis().SetTitle(ytitle)
   hist[name].SetDirectory(0)

def fresh_histos(hist={}):
   hist1D('t0mc', 't0mc', 501, -25, 25,
          "primary vertex thrown time (s)", "events", hist)
   hist1D('trf', 'tRF', 501, -25, 25,
          "primary vertex RF time (ns)", "events", hist)
   hist1D('tmicro', 'tmicro-tRF', 501, -25, 25,
          "tagm time minus RF time (ns)", "events", hist)
   hist1D('thodo', 'thodo-tRF', 501, -25, 25,
          "tagh time minus RF time (ns)", "events", hist)
   hist1D('tfcal1', 'tfcal-tRF ring 1', 501, 0, 50,
          "fcal time minus RF time (ns)", "events", hist)
   hist1D('tfcal2', 'tfcal-tRF ring 2', 501, 0, 50,
          "fcal time minus RF time (ns)", "events", hist)
   hist1D('tfcal3', 'tfcal-tRF ring 3', 501, 0, 50,
          "fcal time minus RF time (ns)", "events", hist)
   hist1D('tfcal4', 'tfcal-tRF ring 4', 501, 0, 50,
          "fcal time minus RF time (ns)", "events", hist)
   hist1D('tfcal3129', 'tfcal-tRF block 31,29', 501, 0, 50,
          "fcal time minus RF time (ns)", "events", hist)
   hist1D('tfcal2729', 'tfcal-tRF block 27,29', 501, 0, 50,
          "fcal time minus RF time (ns)", "events", hist)
   hist1D('tfcal2931', 'tfcal-tRF block 29,31', 501, 0, 50,
          "fcal time minus RF time (ns)", "events", hist)
   hist1D('tfcal2927', 'tfcal-tRF block 29,27', 501, 0, 50,
          "fcal time minus RF time (ns)", "events", hist)
   hist1D('tpsc', 'tpsc-trf', 1000, -50, 50,
          "psc time minus RF time (ns)", "events", hist)
   hist1D('tstc', 'tstc-trf', 1000, -50, 50,
          "start counter time minus RF time (ns)", "events", hist)
   hist1D('tftof', 'tftof-trf', 1000, -50, 50,
          "time of flight counter time minus RF time (ns)", "events", hist)
   return hist

def event_scan(pid, hddmfile, results,
               events=[0,999999999999,1], goodevents=0):
   hist = fresh_histos()
   hddmstream = hddm_s.istream(hddmfile)
   hddmstream.skip(events[0])
   eventlist = range(events[0], events[1], events[2])
   for count in eventlist:
      try:
         rec = next(hddmstream)
      except:
         break
      eventNo = rec.getPhysicsEvent().eventNo
      if goodevents and not eventNo in goodevents:
         continue
      vtx = rec.getVertices()
      zvtx = vtx[0].getOrigin().vz
      tvtx = vtx[0].getOrigin().t
      t0mc = tvtx + (ztarget - zvtx) / clight
      hist['t0mc'].Fill(t0mc)
      for trf in rec.getRFtimes():
         hist['trf'].Fill(trf.tsync)
      for micro in rec.getMicroChannels():
         for tag in micro.getTaggerHits():
            hist['tmicro'].Fill(tag.t - t0mc)
      for hodo in rec.getHodoChannels():
         for tag in hodo.getTaggerHits():
            hist['thodo'].Fill(tag.t - t0mc)
      for fblock in rec.getFcalBlocks():
         r = ((fblock.column - fcal_middle_col)**2 +
              (fblock.row - fcal_middle_row)**2)**0.5
         for fhit in fblock.getFcalHits():
            if r < 3:
               hist['tfcal1'].Fill(fhit.t - t0mc)
            elif r < 5:
               hist['tfcal2'].Fill(fhit.t - t0mc)
            elif r < 7:
               hist['tfcal3'].Fill(fhit.t - t0mc)
            elif r < 9:
               hist['tfcal4'].Fill(fhit.t - t0mc)
            if fblock.row == 29:
               if fblock.column == 27:
                  hist['tfcal2927'].Fill(fhit.t - t0mc)
               elif fblock.column == 31:
                  hist['tfcal2931'].Fill(fhit.t - t0mc)
            elif fblock.column == 29:
               if fblock.row == 27:
                  hist['tfcal2729'].Fill(fhit.t - t0mc)
               elif fblock.row == 31:
                  hist['tfcal3129'].Fill(fhit.t - t0mc)
      for phit in rec.getPscHits():
         hist['tpsc'].Fill(phit.t - t0mc)
      for stchit in rec.getStcHits():
         hist['tstc'].Fill(stchit.t - t0mc)
      for tofhit in rec.getFtofHits():
         hist['tftof'].Fill(tofhit.t - t0mc)
      try:
         hddmstream.skip(events[2] - 1)
      except:
         break
   results[pid] = hist
   return pid

def parallel_scans(nprocs, goodevents=0):
   nfiles = len(sys.argv) - 1
   manager = multiprocessing.Manager()
   pool = multiprocessing.Pool(nprocs)
   ROOT.EnableImplicitMT(nprocs)
   try:
      ftmp = open(sys.argv[1])
      ftmp = 0
   except:
      usage()
   procs = []
   results = manager.dict()
   for ifile in range(nfiles):
      procs_per_file = max(int((nprocs - len(procs)) / (nfiles - ifile)), 1)
      for i in range(procs_per_file):
         pid = len(procs)
         hddmfile = sys.argv[ifile + 1]
         args={'events': [i, 10000, procs_per_file],
               'goodevents': goodevents,
              }
         #print(f"spawning subprocess {pid}")
         proc = pool.apply_async(event_scan, (pid, hddmfile, results), args)
         procs.append(proc)

   pbar = tqdm.tqdm(total=len(procs))

   hist = {}
   for proc in procs:
      pid = proc.get()
      pbar.update()
      #print(f"reaping subprocess {i}")
      if not pid in results:
         print(f"process {pid} returned no results")
         continue
      for h in results[pid]:
         if h in hist:
            hist[h].Add(results[pid][h])
         else:
            hist[h] = results[pid][h]
   return hist
       

if __name__ == "__main__":
   goodevents = {}
   for line in open("tags"):
      evno = int(line.split()[0])
      if evno in goodevents:
         goodevents[evno] += 1
      else:
         goodevents[evno] = 1
   nprocs = 8
   hist = parallel_scans(nprocs, goodevents)
   fout = ROOT.TFile("fcal_timing.root", "recreate")
   for h in hist:
      hist[h].Write()

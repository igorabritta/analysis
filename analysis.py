#!/usr/bin/env python

import os,math,sys
import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gROOT.SetBatch(True)
from channelsTools import cameraChannel

from snakes import SnakesFactory
from output import OutputTree
from treeVars import AutoFillTreeProducer

class analysis:

    def __init__(self,rfile,options):
        self.xmax = 2048
        self.rebin = options.rebin        
        self.rfile = rfile
        self.options = options
        self.pedfile_name = '{base}_ped_rebin{rb}.root'.format(base=os.path.splitext(self.rfile)[0],rb=self.rebin)
        if not os.path.exists(self.pedfile_name):
            print "WARNING: pedestal file ",self.pedfile_name, " not existing. First calculate them..."
            self.calcPedestal(options.numPedEvents)
        self.pedfile_fullres_name = '{base}_ped_rebin1.root'.format(base=os.path.splitext(self.rfile)[0])
        if not os.path.exists(self.pedfile_fullres_name):
            print "WARNING: pedestal file with full resolution ",self.pedfile_fullres_name, " not existing. First calculate them..."
            self.calcPedestal(options.numPedEvents,1)
        if not options.justPedestal:
           print "Pulling pedestals..."
           # first the one for clustering with rebin
           pedrf = ROOT.TFile.Open(self.pedfile_name)
           self.pedmap = pedrf.Get('pedmap').Clone()
           self.pedmap.SetDirectory(None)
           self.pedmean = pedrf.Get('pedmean').GetMean()
           self.pedrms = pedrf.Get('pedrms').GetMean()
           pedrf.Close()
           # then the full resolution one
           pedrf_fr = ROOT.TFile.Open(self.pedfile_fullres_name)
           self.pedmap_fr = pedrf_fr.Get('pedmap').Clone()
           self.pedmap_fr.SetDirectory(None)
           pedrf_fr.Close()

    # the following is needed for multithreading
    def __call__(self,evrange=(-1,-1)):
        return self.reconstruct(evrange)

    def beginJob(self):
        # prepare output file
        self.outputFile = ROOT.TFile.Open(self.options.outFile, "RECREATE")
        # prepare output tree
        self.outputTree = ROOT.TTree("Events","Tree containing reconstructed quantities")
        self.outTree = OutputTree(self.outputFile,self.outputTree)
        self.autotree = AutoFillTreeProducer(self.outTree)

        self.outTree.branch("run", "I")
        self.outTree.branch("event", "I")
        self.autotree.createPMTVariables()

    def endJob(self):
        self.outTree.write()
        self.outputFile.Close()
        
    def getNEvents(self):
        tf = ROOT.TFile.Open(self.rfile)
        ret = len(tf.GetListOfKeys())
        tf.Close()
        return ret

    def calcPedestal(self,maxImages=-1,alternativeRebin=-1):
        nx=ny=self.xmax
        rebin = self.rebin if alternativeRebin<0 else alternativeRebin
        nx=int(nx/rebin); ny=int(ny/rebin); 
        pedfile = ROOT.TFile.Open(self.pedfile_name,'recreate')
        pedmap = ROOT.TProfile2D('pedmap','pedmap',nx,0,self.xmax,ny,0,self.xmax,'s')
        tf = ROOT.TFile.Open(self.rfile)
        for i,e in enumerate(tf.GetListOfKeys()):
            if maxImages>-1 and i<len(tf.GetListOfKeys())-maxImages: continue
            name=e.GetName()
            obj=e.ReadObj()
            if not obj.InheritsFrom('TH2'): continue
            print "Calc pedestal with event: ",name
            obj.RebinX(rebin); obj.RebinY(rebin); 
            for ix in xrange(nx+1):
                for iy in xrange(ny+1):
                    x = obj.GetXaxis().GetBinCenter(ix+1)
                    y = obj.GetYaxis().GetBinCenter(iy+1)
                    pedmap.Fill(x,y,obj.GetBinContent(ix+1,iy+1)/float(math.pow(self.rebin,2)))

        tf.Close()
        pedfile.cd()
        pedmap.Write()
        pedmean = ROOT.TH1D('pedmean','pedestal mean',500,97,103)
        pedrms = ROOT.TH1D('pedrms','pedestal RMS',500,0,5)
        for ix in xrange(nx):
            for iy in xrange(ny):
               pedmean.Fill(pedmap.GetBinContent(ix,iy)) 
               pedrms.Fill(pedmap.GetBinError(ix,iy)) 
        pedmean.Write()
        pedrms.Write()
        pedfile.Close()
        print "Pedestal calculated and saved into ",self.pedfile_name

    def reconstruct(self,evrange=(-1,-1)):
        
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        tf = ROOT.TFile.Open(self.rfile)
        c1 = ROOT.TCanvas('c1','',600,600)
        cc = cameraChannel()
        print "Reconstructing event range: ",evrange
        # loop over events (pictures)
        for iobj,key in enumerate(tf.GetListOfKeys()) :
            iev = iobj/2
            if options.maxEntries>0 and iev==max(evrange[0],0)+options.maxEntries: break
            if sum(evrange)>-2:
                if iev<evrange[0] or iev>evrange[1]: continue
                
            name=key.GetName()
            obj=key.ReadObj()

            if obj.InheritsFrom('TH2'):
                run,event=(int(name.split('_')[1].split('run')[-1].lstrip("0")),int(name.split('_')[-1].split('ev')[-1]))
                print "Processing run: ",run," event ",event,"..."
                self.outTree.fillBranch("run",run)
                self.outTree.fillBranch("event",event)

                pic_fullres = obj.Clone(obj.GetName()+'_fr')
                # rebin to make the cluster finding fast
                obj.RebinX(self.rebin); obj.RebinY(self.rebin)
                obj.Scale(1./float(math.pow(self.rebin,2)))

                # applying zero-suppression
                h2zs = cc.zs(obj,self.pedmap)

                # Cluster reconstruction on 2D picture
                print "Analyzing picture ",name," contours..."
                snfac = SnakesFactory(h2zs,name,options)
                # this plotting is only the pyplot representation.
                # Doesn't work on MacOS with multithreading for some reason... 
                snakes = snfac.getClusters(plot=False)
                snfac.plotClusterFullResolution(snakes,pic_fullres,self.pedmap_fr)
                snfac.plotProfiles(snakes,pic_fullres,self.pedmap_fr)

                # PMT waveform reconstruction
                from waveform import PeakFinder,PMTSignal
                wform = tf.Get('wfm_'+'_'.join(name.split('_')[1:]))
                #pmtsig = PMTSignal(wform,snakes,self.options)
                #pmtsig.plotNice()

                pf = PeakFinder(wform,6160,6300) # time range ok for FNG
                threshold = 10 # min threshold for a signal
                min_distance_peaks = 10 # number of samples (10 samples = 2ns)
                prominence = 2 # noise seems ~1 mV
                width = 5 # minimal width of the signal
                pf.findPeaks(threshold,min_distance_peaks,prominence)
                pf.plotpy(pdir=options.plotDir)
                self.autotree.fillPMTVariables(pf)
                
                # fill reco tree
                self.outTree.fill()
                             
                # DEPRECATED (GAC)
                #snakes = snfac.getContours(iterations=100)
                #snfac.plotContours(snakes,fill=True)
                #snfac.filledSnakes(snakes)

if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog h5file1,...,h5fileN [opts] ')
    parser.add_option("-o","--out", dest="outFile", type="string", default="reco.root", help="name of the output root file");
    parser.add_option('-j', '--jobs', dest='jobs', default=1, type='int', help='Jobs to be run in parallel')
    parser.add_option('-r', '--rebin', dest='rebin', default=4, type='int', help='Rebin factor (same in x and y)')
    parser.add_option(      '--numPedEvents', dest='numPedEvents', default=-1, type='float', help='Use the last n events to calculate the pedestal. Default is all events')
    parser.add_option(      '--max-entries', dest='maxEntries', default=-1, type='float', help='Process only the first n entries')
    parser.add_option(      '--pdir', dest='plotDir', default='./', type='string', help='Directory where to put the plots')
    parser.add_option('-p', '--pedestal', dest='justPedestal', default=False, action='store_true', help='Just compute the pedestals, do not run the analysis')

    (options, args) = parser.parse_args()

    inputf = args[0]

    if options.justPedestal:
        ana = analysis(inputf,options)
        print "Pedestals with rebin factor = ",options.rebin, "done. Exiting."
        sys.exit(0)
        
    ana = analysis(inputf,options)
    nev = ana.getNEvents()
    print "This run has ",nev," events."
    print "Will save plots to ",options.plotDir
    ana.beginJob()
    
    if options.jobs>1:
        print "WARNING! With multiprocessing the tree filling is not working yet! Only do reconstruction and plotting"
        nj = int(nev/options.jobs)
        chunks = [(i,min(i+nj,nev)) for i in xrange(0,nev,nj)]
        print chunks
        from multiprocessing import Pool
        pool = Pool(options.jobs)
        ret = pool.map(ana, chunks)
    else:
        ana.reconstruct()

    ana.endJob()

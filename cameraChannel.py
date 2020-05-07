#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
import numpy as np
from root_numpy import hist2array
import math
import debug_code.tools_lib as tl

class cameraGeometry:
    def __init__(self):
        # 240 mm = major axis of the ellipse = 2048 pixels
        # ORANGE
        # 1 px - 55*10-6 m
        # 1 ph/ev - 440 V

        # LEMON
        # 1 px - 125*10-6 m
        # 0,12 ph/ev - 460V
        # 0,10 ph/ev - 450V
        # 0,06 ph/ev - 440V
        
        self.pixelwidth = 1#55E-3#125E-3 # mm
        
class cameraTools:
    def __init__(self):
        pass

    def pedsub(self,img,pedarr):
        return img - pedarr
    
    def zsfullres(self,img_sub,noisearr,nsigma=1):
        img_zs = np.where(img_sub > nsigma * noisearr, img_sub, 0)
        return img_zs
        
    def arrrebin(self,img,rebin):
        newshape = int(self.options.pixels/rebin)
        img_rebin = tl.rebin(img,(newshape,newshape))
        return img_rebin
        
    # this returns x,y,z before the zero suppression
    def getImage(self,th2):
        return hist2array(th2)

    def noisearray(self,th2):
        noisearr = np.zeros( (th2.GetNbinsX(),th2.GetNbinsY()) )
        for ix in range(th2.GetNbinsX()):
            for iy in range(th2.GetNbinsY()):
                noisearr[ix][iy] = th2.GetBinError(ix+1,iy+1)
        return noisearr

    def getRestrictedImage(self,th2,xmin,xmax,ymin,ymax):
        nx = th2.GetNbinsX(); ny = th2.GetNbinsY();
        nxp = xmax-xmin
        nyp = ymax-ymin
        th2_rs = ROOT.TH2D(th2.GetName()+'_rs',th2.GetName()+'_rs',nxp,xmin,xmax,nyp,ymin,ymax)
        th2_rs.SetDirectory(None)
        for ix,x in enumerate(range(xmin,xmax)):
            for iy,y in enumerate(range(ymin,ymax)):
                orig_ixb = th2.GetXaxis().FindBin(x)
                orig_iyb = th2.GetYaxis().FindBin(y)
                z = th2.GetBinContent(orig_ixb,orig_iyb)
                th2_rs.SetBinContent(ix+1,iy+1,z)
        return th2_rs


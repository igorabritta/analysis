import os, math, optparse, ROOT
from array import array

ROOT.gStyle.SetOptStat(111111)
ROOT.gROOT.SetBatch(True)

def doLegend(histos,labels,styles,corner="TR",textSize=0.035,legWidth=0.18,legBorder=False,nColumns=1):
    nentries = len(histos)
    (x1,y1,x2,y2) = (.85-legWidth, .7 - textSize*max(nentries-3,0), .90, .91)
    if corner == "TR":
        (x1,y1,x2,y2) = (.85-legWidth, .7 - textSize*max(nentries-3,0), .90, .91)
    elif corner == "TC":
        (x1,y1,x2,y2) = (.5, .75 - textSize*max(nentries-3,0), .5+legWidth, .91)
    elif corner == "TL":
        (x1,y1,x2,y2) = (.2, .75 - textSize*max(nentries-3,0), .2+legWidth, .91)
    elif corner == "BR":
        (x1,y1,x2,y2) = (.85-legWidth, .33 + textSize*max(nentries-3,0), .90, .15)
    elif corner == "BC":
        (x1,y1,x2,y2) = (.5, .33 + textSize*max(nentries-3,0), .5+legWidth, .35)
    elif corner == "BL":
        (x1,y1,x2,y2) = (.2, .33 + textSize*max(nentries-3,0), .33+legWidth, .35)
    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetNColumns(nColumns)
    leg.SetFillColor(0)
    leg.SetFillColorAlpha(0,0.6)  # should make the legend semitransparent (second number is 0 for fully transparent, 1 for full opaque)
    #leg.SetFillStyle(0) # transparent legend, so it will not cover plots (markers of legend entries will cover it unless one changes the histogram FillStyle, but this has other effects on color, so better not touching the FillStyle)
    leg.SetShadowColor(0)
    if not legBorder:
        leg.SetLineColor(0)
        leg.SetBorderSize(0)  # remove border  (otherwise it is drawn with a white line, visible if it overlaps with plots
    leg.SetTextFont(42)
    leg.SetTextSize(textSize)
    for (plot,label,style) in zip(histos,labels,styles): leg.AddEntry(plot,label,style)
    leg.Draw()
    ## assign it to a global variable so it's not deleted
    global legend_
    legend_ = leg
    return leg


def plotDensity():
    tf = ROOT.TFile('reco_run815.root')
    tree = tf.Get('Events')

    histos = []
    colors = [ROOT.kRed,ROOT.kBlue,ROOT.kOrange]
    #histo = ROOT.TH1F('density','',70,0,2500)
    histo = ROOT.TH1F('density','',70,0,20)
    for it in xrange(1,4):
        h = histo.Clone('h_iter{it}'.format(it=it))
        h.Sumw2()
        #tree.Draw('track_integral/track_nhits>>h_iter{it}'.format(it=it),'track_iteration=={it}'.format(it=it))
        #tree.Draw('track_integral/track_length>>h_iter{it}'.format(it=it),'track_iteration=={it}'.format(it=it))
        tree.Draw('track_length>>h_iter{it}'.format(it=it),'track_iteration=={it}'.format(it=it))
        h.Scale(1./h.Integral())
        h.SetFillColor(colors[it-1])
        h.SetLineColor(colors[it-1])
        h.SetFillStyle(3005)
        histos.append(h)


    # legend
    (x1,y1,x2,y2) = (0.7, .70, .9, .87)
    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetFillColor(0)
    leg.SetFillColorAlpha(0,0.6)
    leg.SetShadowColor(0)
    leg.SetLineColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    
    c = ROOT.TCanvas('c','',600,600)
    for ih,h in enumerate(histos):
        h.Draw('hist' if ih==0 else 'hist same')
        h.GetYaxis().SetRangeUser(0,0.25)
        h.GetXaxis().SetTitle('length (mm)')
        leg.AddEntry(h,'iteration {it}'.format(it=ih+1),'f')

    leg.Draw()

    c.SaveAs('density.pdf')


def plotNClusters(iteration=1):

    f = ROOT.TFile.Open('../runs/fng_runs.root')
    nclu_h = ROOT.TH1F('nclu_h','',20,0,80)
    
    for ie,event in enumerate(f.Events):
        #if event.run<46 or event.run>47: continue
        #if event.run<48 or event.run>50: continue
        if event.run<70 or event.run>71: continue
        nclu_it1 = 0
        for icl in range(event.nTrack):
            if event.track_nhits[icl]<100: continue
            if int(event.track_iteration[icl])>2: continue
            if math.hypot(event.track_xmean[icl]-1024,event.track_ymean[icl]-1024)>1000: continue
            #print "cluster x = ", event.track_xmean[icl]
            #print "cluster y = ", event.track_ymean[icl]
            nclu_it1 += 1
        nclu_h.Fill(nclu_it1)

    c = ROOT.TCanvas('c1','',600,400)
    nclu_h.SetLineColor(ROOT.kRed+2)
    nclu_h.SetMarkerColor(ROOT.kBlack)
    nclu_h.SetMarkerSize(1.5)
    nclu_h.SetMarkerStyle(ROOT.kOpenCircle)
    nclu_h.SetLineWidth(2)
    nclu_h.GetXaxis().SetTitle('# clusters / 100 ms')
    nclu_h.GetYaxis().SetTitle('Events')
    nclu_h.Draw('pe')
    nclu_h.Draw('hist same')
    c.SaveAs("nclusters_iter1_run70_71.pdf")


def withinFC(x,y,ax=400,ay=500,shape=2048):
    center = shape/2.
    x1 = x-center
    y1 = (y-center)*1.2
    return math.hypot(x1,y1)<ax

def slimnessCut(l,w,th=0.6):
    
        return (w/l)>th

def integralCut(i,minTh=1800,maxTh=3200):
    
        return (i>minTh) & (i<maxTh)
    
def varChoice(var):
    
    if var == 'integral':
        var1 = 'cl_integral'
        leg = 'Integral [ph]'
        histlimit = 6000
    elif var == 'length':
        var1 = 'cl_length*125E-3'
        leg = 'Length [mm]'
        histlimit = 32
    elif var == 'width':
        var1 = 'cl_width*125E-3'
        leg = 'Width [mm]'
        histlimit = 10
    elif var == 'size':
        var1 = 'cl_size'
        leg = 'Size [px]'
        histlimit = 1600
    elif var == 'slimness':
        var1 = 'cl_width/cl_length'
        leg = 'Slimness [w/l]'
        histlimit = 1.2
    elif var == 'intsize':
        var1 = '(cl_integral/cl_size)'
        leg = 'Integral/Size [ph/px]'
        histlimit = 10
    else:
        exit()
    
    return var1, leg, histlimit
    

def fillSpectra(cluster='cl'):

    ret = {}
    tf_cosmics = ROOT.TFile('../reco_run02165_3D.root')
    tf_fe55 = ROOT.TFile('../reco_run02279_3D.root')
    tfiles = {'fe':tf_fe55,'cosm':tf_cosmics}
    
    ## signal region histograms
    ret[('fe','integral')] = ROOT.TH1F("fe_integral",'',100,-10,8000)
    ret[('fe','length')]   = ROOT.TH1F("fe_length",'',100,0,500)
    ret[('fe','width')]    = ROOT.TH1F("fe_width",'',35,0,70)
    ret[('fe','nhits')]    = ROOT.TH1F("fe_nhits",'',100,0,1500)
    ret[('fe','slimness')] = ROOT.TH1F("fe_slimness",'',100,0,1)

    # x-axis titles
    titles = {'integral': 'photons', 'length':'length (pixels)', 'width':'width (pixels)', 'nhits': 'active pixels', 'slimness': 'width/length'}
    
    ## control region histograms
    ret2 = {}
    for (region,var),h in ret.iteritems():
        ret[(region,var)].Sumw2()
        ret[(region,var)].SetDirectory(None)
        ret[(region,var)].GetXaxis().SetTitle(titles[var])
        ret[(region,var)].GetXaxis().SetTitleSize(0.1)
        ret2[('cosm',var)] = h.Clone('cosm_{name}'.format(name=var))
        ret2[('cosm',var)].SetDirectory(None)
    ret.update(ret2)

    ## now fill the histograms 
    for runtype in ['fe','cosm']:
        for ie,event in enumerate(tfiles[runtype].Events):
            for isc in range(getattr(event,"nSc" if cluster=='sc' else 'nCl')):
                if getattr(event,"{clutype}_iteration".format(clutype=cluster))[isc]!=2:
                    continue
                if not withinFC(getattr(event,"{clutype}_xmean".format(clutype=cluster))[isc],getattr(event,"{clutype}_ymean".format(clutype=cluster))[isc]):
                    continue
                #if not slimnessCut(getattr(event,"{clutype}_length".format(clutype=cluster))[isc],getattr(event,"{clutype}_width".format(clutype=cluster))[isc]):
                 #   continue
                #if not integralCut(getattr(event,"{clutype}_integral".format(clutype=cluster))[isc]):
                #    continue
                for var in ['integral','length','width','nhits']:
                    ret[(runtype,var)].Fill(getattr(event,("{clutype}_{name}".format(clutype=cluster,name=var)))[isc])
                ret[(runtype,'slimness')].Fill(getattr(event,"{clutype}_width".format(clutype=cluster))[isc] / getattr(event,"{clutype}_length".format(clutype=cluster))[isc])

    return ret

def getCanvas():
    c = ROOT.TCanvas('c','',1200,1200)
    lMargin = 0.14
    rMargin = 0.10
    bMargin = 0.15
    tMargin = 0.10
    c.SetLeftMargin(lMargin)
    c.SetRightMargin(rMargin)
    c.SetTopMargin(tMargin)
    c.SetBottomMargin(bMargin)
    c.SetFrameBorderMode(0);
    c.SetBorderMode(0);
    c.SetBorderSize(0);
    return c

def drawOne(histo_sr,histo_cr,plotdir='./'):
    ROOT.gStyle.SetOptStat(0)
    
    c = ROOT.TCanvas('c','',1200,1200)
    lMargin = 0.12
    rMargin = 0.05
    bMargin = 0.30
    tMargin = 0.07
    padTop = ROOT.TPad('padTop','',0.,0.4,1,0.98)
    padTop.SetLeftMargin(lMargin)
    padTop.SetRightMargin(rMargin)
    padTop.SetTopMargin(tMargin)
    padTop.SetBottomMargin(0)
    padTop.SetFrameBorderMode(0);
    padTop.SetBorderMode(0);
    padTop.SetBorderSize(0);
    padTop.Draw()

    padBottom = ROOT.TPad('padBottom','',0.,0.02,1,0.4)
    padBottom.SetLeftMargin(lMargin)
    padBottom.SetRightMargin(rMargin)
    padBottom.SetTopMargin(0)
    padBottom.SetBottomMargin(bMargin)
    padBottom.SetFrameBorderMode(0);
    padBottom.SetBorderMode(0);
    padBottom.SetBorderSize(0);
    padBottom.Draw()

    padTop.cd()
    histo_cr.SetMaximum(1.2*max(histo_cr.GetMaximum(),histo_sr.GetMaximum()))
    histo_cr.GetYaxis().SetTitle('clusters (a.u.)')
    histo_cr.Draw("hist")
    histo_sr.Draw("pe same")

    histos = [histo_sr,histo_cr]
    labels = ['^{55}Fe source 70/30 - Pos 5','^{55}Fe source 60/40 - Pos 5']
    styles = ['p','f']
    
    legend = doLegend(histos,labels,styles,corner="TR")
    legend.Draw()
    
    padBottom.cd()
    ratio = histo_sr.Clone(histo_sr.GetName()+"_ratio")
    ratio.Divide(histo_cr)
    ratio.GetYaxis().SetTitleSize(0.05)
    ratio.GetYaxis().SetTitle("{num} / {den}".format(num=labels[0],den=labels[1]))
    ratio.Draw('pe1')

    line = ROOT.TLine()
    line.DrawLine(ratio.GetXaxis().GetBinLowEdge(1), 1, ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1), 1)
    line.SetLineStyle(3)
    line.SetLineColor(ROOT.kBlack)
    
    for ext in ['png','pdf','root']:
        c.SaveAs("{plotdir}/{var}.{ext}".format(plotdir=plotdir,var=histo_sr.GetName(),ext=ext))
    
    
def drawSpectra(histos,plotdir):
    variables = [var for (reg,var) in histos.keys() if reg=='fe']
    print "variables to plot: ", variables

    for var in variables:
        histos[('cosm',var)].SetFillColor(ROOT.kAzure+6)
        histos[('fe',var)].SetMarkerStyle(ROOT.kFullDotLarge)
        histos[('cosm',var)].Scale(histos[('fe',var)].Integral()/histos[('cosm',var)].Integral())
        drawOne(histos[('fe',var)],histos[('cosm',var)],plotdir)

def plotEnergyVsDistance(plotdir):

    tf_fe55 = ROOT.TFile('runs/reco_run01740_3D.root')
    tree = tf_fe55.Get('Events')

    np = 9 # number of points
    dist = 100 # distance from center in pixels
    x = [dist*i for i in range(np+1)]
    y_mean = []; y_res = []

    cut_base = 'cl_iteration==2'
    integral = ROOT.TH1F('integral','',100,600,3000)
    
    for p in range(np):
        cut = "{base} && TMath::Hypot(cl_xmean-1024,cl_ymean-1024)>{r_min} && TMath::Hypot(cl_xmean-1024,cl_ymean-1024)<{r_max}".format(base=cut_base,r_min=x[p],r_max=x[p+1])
        tree.Draw('cl_integral>>integral',cut)
        mean = integral.GetMean()
        rms  = integral.GetRMS()
        y_mean.append(mean)
        y_res.append(rms/mean)
        integral.Reset()

    gr_mean = ROOT.TGraph(np,array('f',x),array('f',y_mean))
    gr_res = ROOT.TGraph(np,array('f',x),array('f',y_res))
    gr_mean.SetMarkerStyle(ROOT.kOpenCircle)
    gr_res.SetMarkerStyle(ROOT.kOpenCircle)
    gr_mean.SetMarkerSize(2)
    gr_res.SetMarkerSize(2)
    
    gr_mean.SetTitle('')
    gr_res.SetTitle('')
    
    c = ROOT.TCanvas('c','',1200,1200)
    lMargin = 0.12
    rMargin = 0.05
    bMargin = 0.30
    tMargin = 0.07

    gr_mean.Draw('AP')
    gr_mean.GetXaxis().SetTitle('distance from center (pixels)')
    gr_mean.GetYaxis().SetTitle('integral (photons)')

    for ext in ['png','pdf']:
        c.SaveAs("{plotdir}/mean.{ext}".format(plotdir=plotdir,ext=ext))

    gr_res.Draw('AP')
    gr_res.GetXaxis().SetTitle('distance from center (pixels)')
    gr_res.GetYaxis().SetTitle('resolution (rms)')
    gr_res.GetYaxis().SetRangeUser(0.10,0.50)

    for ext in ['png','pdf']:
        c.SaveAs("{plotdir}/rms.{ext}".format(plotdir=plotdir,ext=ext))

    print x
    print y_mean
    print y_res

def plotCameraEnergyVsPosition(plotdir,var='integral',gas = 60):
    
    x = []
    y_mean = []; y_res = []
    
    #gas = 70
    
    if gas == 70:
        gg   = '70-30'
        n    = 7
    else:
        gg   = '60-40'
        n    = 6
    
    for i in range(0,n):
        
        mean, rms, i, leg = plotHistFit(plotdir,var,i,gas)
        x.append(i)
        y_mean.append(mean)
        y_res.append(rms/mean)

    print y_res
    np = len(x)
    
    gr_mean = ROOT.TGraph(np,array('f',x),array('f',y_mean))
    gr_res = ROOT.TGraph(np,array('f',x),array('f',y_res))
    gr_mean.GetXaxis().SetRangeUser(0,21)
    gr_res.GetXaxis().SetRangeUser(0,21)
    gr_mean.GetYaxis().SetRangeUser(0,max(y_mean)*1.2)
    gr_res. GetYaxis().SetRangeUser(0,max(y_res)*1.2)
    gr_mean.SetMarkerStyle(ROOT.kOpenCircle)
    gr_res.SetMarkerStyle(ROOT.kOpenCircle)
    gr_mean.SetMarkerSize(2)
    gr_res.SetMarkerSize(2)
    
    gr_mean.SetTitle('')
    gr_res.SetTitle('')

    c = ROOT.TCanvas('c','',1200,1200)
    lMargin = 0.17
    rMargin = 0.05
    bMargin = 0.15
    tMargin = 0.07
    c.SetLeftMargin(lMargin)
    c.SetRightMargin(rMargin)
    c.SetTopMargin(tMargin)
    c.SetBottomMargin(bMargin)
    c.SetFrameBorderMode(0);
    c.SetBorderMode(0);
    c.SetBorderSize(0);
    c.SetGrid();

    gr_mean.Draw('AP')
    gr_mean.GetXaxis().SetTitle('Source Position [cm]')
    gr_mean.GetYaxis().SetTitle('{leg}'.format(leg=leg))

    for ext in ['pdf']:
        c.SaveAs("{plotdir}/CameraMean_{var}_{gg}.{ext}".format(plotdir=plotdir,var=var,ext=ext,gg=gg))

    gr_res.Draw('AP')
    gr_res.GetXaxis().SetTitle('Source Position [cm]')
    gr_res.GetYaxis().SetTitle('Resolution [rms]')

    for ext in ['pdf']:
        c.SaveAs("{plotdir}/CameraRMS_{var}_{gg}.{ext}".format(plotdir=plotdir,var=var,ext=ext,gg=gg))
    
    return c, gr_mean, gr_res

def calcCameraEnergyVsPosition(plotdir,var='integral',gas = 60):
    
    x = []
    y_mean = []; y_res = []; y_spot = []
    
    if gas == 70:
        n    = 7
    else:
        n    = 6
      
    for i in range(0,n):
        
        mean, rms, i, leg, spotn = plotHistFit(plotdir,var,i,gas)
        x.append(i)
        y_mean.append(mean)
        y_res.append(rms/mean)
        y_spot.append(spotn)
    
    return x, y_mean, y_res, leg, y_spot


def plotCameraEnergyVsPositionGas(plotdir,var='integral'):
    
    x1, y_mean1, y_res1, leg, y_spot1 = calcCameraEnergyVsPosition(plotdir, var, gas = 60)
    x2, y_mean2, y_res2, leg, y_spot2 = calcCameraEnergyVsPosition(plotdir, var, gas = 70)
    
    
    y_mean =y_mean1 + y_mean2
    y_res = y_res1 + y_res2
    y_spot = y_spot1 + y_spot2
    
    ## Creating the CANVAS
    
    c = ROOT.TCanvas('c','',1200,1200)
    lMargin = 0.17
    rMargin = 0.05
    bMargin = 0.15
    tMargin = 0.07
    c.SetLeftMargin(lMargin)
    c.SetRightMargin(rMargin)
    c.SetTopMargin(tMargin)
    c.SetBottomMargin(bMargin)
    c.SetFrameBorderMode(0);
    c.SetBorderMode(0);
    c.SetBorderSize(0);
    c.SetGrid();   
    
    
    ## Creating MultiGraph
    mg1 = ROOT.TMultiGraph("mg","");
    mg2 = ROOT.TMultiGraph("mg2","");
    
    if var=='integral':
        mg3 = ROOT.TMultiGraph("mg3","");
        
        np3 = len(x1)
        ## Creating TGraph SPOT

        gr_spot1 = ROOT.TGraph(np3,array('f',x1),array('f',y_spot1))
        gr_spot1.SetName("He/CF_{4} - 60/40 ");
        gr_spot1.SetDrawOption("AP");
        gr_spot1.GetXaxis().SetRangeUser(0,21)
        gr_spot1.GetYaxis().SetRangeUser(0,max(y_spot)*1.2)
        gr_spot1.SetMarkerStyle(ROOT.kOpenCircle)
        gr_spot1.SetMarkerSize(2)
        gr_spot1.SetTitle('')
        gr_spot1.GetXaxis().SetTitle('Source Position [cm]')
        gr_spot1.GetYaxis().SetTitle('Spot Number [Counts]')
        
        np4 = len(x2)
        ## Creating TGraph SPOT

        gr_spot2 = ROOT.TGraph(np4,array('f',x2),array('f',y_spot2))
        gr_spot2.SetName("He/CF_{4} - 70/30 ");
        gr_spot2.SetDrawOption("AP");
        gr_spot2.GetXaxis().SetRangeUser(0,21)
        gr_spot2.GetYaxis().SetRangeUser(0,max(y_spot)*1.2)
        gr_spot2.SetMarkerStyle(ROOT.kOpenSquare)
        gr_spot2.SetMarkerSize(2)
        gr_spot2.SetTitle('')
        gr_spot2.GetXaxis().SetTitle('Source Position [cm]')
        gr_spot2.GetYaxis().SetTitle('Spot Number [Counts]')
        
        mg3.Add( gr_spot1 );
        gr_spot2.Draw("ALP");
        mg3.Draw("LP");
        c.BuildLegend(0.3,0.3,0.1,0.1);
    
        for ext in ['pdf']:
            c.SaveAs("{plotdir}/CameraGasCompare_SpotNumber.{ext}".format(plotdir=plotdir,var=var,ext=ext))
    
    # 60-40
    
    np1 = len(x1)
    ## Creating TGraph MEAN
    
    gr_mean1 = ROOT.TGraph(np1,array('f',x1),array('f',y_mean1))
    gr_mean1.SetName("He/CF_{4} - 60/40 ");
    gr_mean1.SetDrawOption("AP");
    gr_mean1.GetXaxis().SetRangeUser(0,21)
    gr_mean1.GetYaxis().SetRangeUser(0,max(y_mean)*1.2)
    gr_mean1.SetMarkerStyle(ROOT.kOpenCircle)
    gr_mean1.SetMarkerSize(2)
    gr_mean1.SetTitle('')
    gr_mean1.GetXaxis().SetTitle('Source Position [cm]')
    gr_mean1.GetYaxis().SetTitle('{leg}'.format(leg=leg))
    
    ## Creating TGraph RES
    
    gr_res1 = ROOT.TGraph(np1,array('f',x1),array('f',y_res1))
    gr_res1.SetName("He/CF_{4} - 60/40 ");
    gr_res1.SetDrawOption("AP");
    gr_res1.GetXaxis().SetRangeUser(0,21)
    gr_res1.GetYaxis().SetRangeUser(0,max(y_res)*1.2)
    gr_res1.SetMarkerStyle(ROOT.kOpenCircle)
    gr_res1.SetMarkerSize(2)
    gr_res1.SetTitle('')
    gr_res1.GetXaxis().SetTitle('Source Position [cm]')
    gr_res1.GetYaxis().SetTitle('Resolution [rms]')


    # 70-30
    np2 = len(x2)
    gr_mean2 = ROOT.TGraph(np2,array('f',x2),array('f',y_mean2))
    gr_mean2.SetName("He/CF_{4} - 70/30 ");
    gr_mean2.SetDrawOption("AP");
    gr_mean2.GetXaxis().SetRangeUser(0,21)
    gr_mean2.GetYaxis().SetRangeUser(0,max(y_mean)*1.2)
    gr_mean2.SetMarkerStyle(ROOT.kOpenSquare)
    gr_mean2.SetMarkerSize(2)
    gr_mean2.SetTitle('')
    gr_mean2.GetXaxis().SetTitle('Source Position [cm]')
    gr_mean2.GetYaxis().SetTitle('{leg}'.format(leg=leg))
    
    gr_res2 = ROOT.TGraph(np2,array('f',x2),array('f',y_res2))
    gr_res2.SetName("He/CF_{4} - 70/30 ");
    gr_res2.SetDrawOption("AP");
    gr_res2.GetXaxis().SetRangeUser(0,21)
    gr_res2.GetYaxis().SetRangeUser(0,max(y_res)*1.2)
    gr_res2.SetMarkerStyle(ROOT.kOpenSquare)
    gr_res2.SetMarkerSize(2)
    gr_res2.SetTitle('')
    gr_res2.GetXaxis().SetTitle('Source Position [cm]')
    gr_res2.GetYaxis().SetTitle('Resolution [rms]')

    
    #### RATIO
    if var=='integral':
        import numpy as np
        yratio = ((np.array(y_mean1)/np.array(y_mean2[0:6]))*1100).tolist()
        
        gr_ratio1 = ROOT.TGraph(np1,array('f',x1),array('f',yratio))
        gr_ratio1.SetName("Ratio 60/40 - 70/30 ");
        gr_ratio1.SetDrawOption("AP");
        gr_ratio1.GetXaxis().SetRangeUser(0,21)
        gr_ratio1.GetYaxis().SetRangeUser(0,max(y_mean)*1.2)
        gr_ratio1.SetMarkerStyle(ROOT.kOpenTriangleUp)
        gr_ratio1.SetMarkerSize(2)
        gr_ratio1.SetTitle('')
        gr_ratio1.GetXaxis().SetTitle('Source Position [cm]')
        gr_ratio1.GetYaxis().SetTitle('{leg}'.format(leg=leg))  
        mg1.Add( gr_ratio1 );
        
    mg1.Add( gr_mean1 );
    gr_mean2.Draw("ALP");
    mg1.Draw("LP");
    c.BuildLegend();
    
    for ext in ['pdf']:
        c.SaveAs("{plotdir}/CameraGasCompareMean_{var}.{ext}".format(plotdir=plotdir,var=var,ext=ext))
        
    mg2.Add( gr_res1 );
    gr_res2.Draw("ALP");
    mg2.Draw("LP");
    c.BuildLegend();
    
    for ext in ['pdf']:
        c.SaveAs("{plotdir}/CameraGasCompareRMS_{var}.{ext}".format(plotdir=plotdir,var=var,ext=ext))    
    
def plotPMTEnergyVsPosition(plotdir,gas=60):
    import numpy as np
      
    if gas == 70:
        tf_fe55 = ROOT.TFile('../reco_run02274_to_run02280_pmt.root')
        pos  = [0, 1, 2, 3, 4, 5, 6]
        run  = [2274, 2275, 2276, 2277, 2278, 2279, 2280]
        dist = (23-np.array([19.5, 16.5, 14.3, 12.5, 10.5, 8.2, 6.2])).tolist()
        gg   = '70-30'
    else:
        tf_fe55 = ROOT.TFile('../reco_run02160_to_run02165_pmt.root')
        pos  = [0, 1, 2, 3, 4, 5]
        run  = [2160, 2161, 2162, 2163, 2164, 2165]
        dist = (23-np.array([19.5, 16.5, 14.3, 12.5, 10.5, 8.2])).tolist()
        gg   = '60-40'

    tree = tf_fe55.Get('Events')
    integral = ROOT.TH1F('integral','',100,2000,20000)
        
    np = len(run)+1
    x = []
    y_mean = []; y_res = []
    
    for i,r in enumerate(run):
        cut = 'run=={r} && pmt_tot<100'.format(r=r)
        tree.Draw('pmt_integral>>integral',cut)
        mean = integral.GetMean()
        rms  = integral.GetRMS()
        x.append(i)
        y_mean.append(mean)
        y_res.append(rms/mean)
        integral.Reset()

    return y_mean,y_res,dist

def plotPMT2gas(plotdir):    

    y_mean1, y_res1, x1 = plotPMTEnergyVsPosition(plotdir,gas=60)
    y_mean2, y_res2, x2 = plotPMTEnergyVsPosition(plotdir,gas=70)
    
    y_mean =y_mean1 + y_mean2
    y_res = y_res1 + y_res2
    
    c = ROOT.TCanvas('c','',1200,1200)
    lMargin = 0.17
    rMargin = 0.05
    bMargin = 0.15
    tMargin = 0.07
    c.SetLeftMargin(lMargin)
    c.SetRightMargin(rMargin)
    c.SetTopMargin(tMargin)
    c.SetBottomMargin(bMargin)
    c.SetFrameBorderMode(0);
    c.SetBorderMode(0);
    c.SetBorderSize(0);
    c.SetGrid();
    
    ## Creating MultiGraph
    mg1 = ROOT.TMultiGraph("mg","");
    mg2 = ROOT.TMultiGraph("mg2","");
    
    
    np1 = len(x1)
    ## Creating TGraph MEAN
    
    gr_mean1 = ROOT.TGraph(np1,array('f',x1),array('f',y_mean1))
    gr_mean1.SetName("He/CF_{4} - 60/40 ");
    gr_mean1.SetDrawOption("AP");
    gr_mean1.GetXaxis().SetRangeUser(0,21)
    gr_mean1.GetYaxis().SetRangeUser(0,max(y_mean)*1.2)
    gr_mean1.SetMarkerStyle(ROOT.kOpenCircle)
    gr_mean1.SetMarkerSize(2)
    gr_mean1.SetTitle('')
    gr_mean1.GetXaxis().SetTitle('Source Position [cm]')
    gr_mean1.GetYaxis().SetTitle('PMT integral [mV]')
    
    ## Creating TGraph RES
    
    gr_res1 = ROOT.TGraph(np1,array('f',x1),array('f',y_res1))
    gr_res1.SetName("He/CF_{4} - 60/40 ");
    gr_res1.SetDrawOption("AP");
    gr_res1.GetXaxis().SetRangeUser(0,21)
    gr_res1.GetYaxis().SetRangeUser(0,max(y_res)*1.2)
    gr_res1.SetMarkerStyle(ROOT.kOpenCircle)
    gr_res1.SetMarkerSize(2)
    gr_res1.SetTitle('')
    gr_res1.GetXaxis().SetTitle('Source Position [cm]')
    gr_res1.GetYaxis().SetTitle('Resolution [rms]')
    
    
    # 70-30
    np2 = len(x2)
    
    gr_mean2 = ROOT.TGraph(np2,array('f',x2),array('f',y_mean2))
    gr_mean2.SetName("He/CF_{4} - 70/30 ");
    gr_mean2.SetDrawOption("AP");
    gr_mean2.GetXaxis().SetRangeUser(0,21)
    gr_mean2.GetYaxis().SetRangeUser(0,max(y_mean)*1.2)
    gr_mean2.SetMarkerStyle(ROOT.kOpenSquare)
    gr_mean2.SetMarkerSize(2)
    gr_mean2.SetTitle('')
    gr_mean2.GetXaxis().SetTitle('Source Position [cm]')
    gr_mean2.GetYaxis().SetTitle('PMT integral [mV]')
    
    gr_res2 = ROOT.TGraph(np2,array('f',x2),array('f',y_res2))
    gr_res2.SetName("He/CF_{4} - 70/30 ");
    gr_res2.SetDrawOption("AP");
    gr_res2.GetXaxis().SetRangeUser(0,21)
    gr_res2.GetYaxis().SetRangeUser(0,max(y_res)*1.2)
    gr_res2.SetMarkerStyle(ROOT.kOpenSquare)
    gr_res2.SetMarkerSize(2)
    gr_res2.SetTitle('')
    gr_res2.GetXaxis().SetTitle('Source Position [cm]')
    gr_res2.GetYaxis().SetTitle('Resolution [rms]')

    mg1.Add( gr_mean1 );
    gr_mean2.Draw("ALP");
    mg1.Draw("LP");
    c.BuildLegend();

    for ext in ['pdf']:
        c.SaveAs("{plotdir}/PMTZScanGasComp_mean.{ext}".format(plotdir=plotdir,ext=ext))

    mg2.Add( gr_res1 );
    gr_res2.Draw("ALP");
    mg2.Draw("LP");
    c.BuildLegend();

    for ext in ['pdf']:
        c.SaveAs("{plotdir}/PMTZScanGasComp_rms.{ext}".format(plotdir=plotdir,ext=ext))

    #print x
    #print y_mean
    #print y_res


def plotCameraPMTCorr(outdir):
    tf_fe55 = ROOT.TFile('runs/reco_run01754_to_run01759.root')
    tree = tf_fe55.Get('Events')

    tot_vs_nhits = ROOT.TH2F('nhits_vs_tot','',45,20,100,30,100,400)
    tot_vs_nhits.GetXaxis().SetTitle("T.o.T. (ms)")
    tot_vs_nhits.GetYaxis().SetTitle("supercluster pixels")
    tot_vs_nhits.SetContour(100)
    
    ## fill the 2D histogram
    for event in tree:
        if event.pmt_tot > 100: continue
        for isc in range(event.nSc):
            if event.sc_iteration[isc]!=2: continue
            if event.sc_width[isc]/event.sc_length[isc]<0.7: continue
            if not withinFC(event.sc_xmean[isc],event.sc_ymean[isc],700,700): continue
            tot_vs_nhits.Fill(event.pmt_tot,event.sc_nhits[isc])

    ## profile for better understanding
    profX = tot_vs_nhits.ProfileX()
    profX.SetMarkerStyle(ROOT.kFullCircle)
    profX.SetLineColor(ROOT.kBlack)
    profX.SetMarkerColor(ROOT.kBlack)
    profX.GetYaxis().SetRangeUser(180,310)
    profX.GetYaxis().SetTitle("average pixels in SC")

    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    ROOT.gStyle.SetOptStat(0)
    c = getCanvas()
    tot_vs_nhits.Draw("colz")
    for ext in ['png','pdf','root']:
        c.SaveAs('{plotdir}/{name}.{ext}'.format(plotdir=outdir,name=tot_vs_nhits.GetName(),ext=ext))

    profX.Draw('pe1')
    for ext in ['png','pdf','root']:
        c.SaveAs('{plotdir}/{name}_profX.{ext}'.format(plotdir=outdir,name=tot_vs_nhits.GetName(),ext=ext))
        
def plotHistFit(plotdir,var='integral',i=0, gas=60):
    import numpy as np
    
    ROOT.gStyle.SetOptFit(1011)
    ROOT.gStyle.SetErrorX(0)
    
    #gas=60
    if gas == 70:
        pos  = [0, 1, 2, 3, 4, 5, 6]
        run  = [2274, 2275, 2276, 2277, 2278, 2279, 2280]
        dist = (23-np.array([19.5, 16.5, 14.3, 12.5, 10.5, 8.2, 6.2])).tolist()
        gg   = '70-30'
    else:
        pos  = [0, 1, 2, 3, 4, 5]
        run  = [2160, 2161, 2162, 2163, 2164, 2165]
        dist = (23-np.array([19.5, 16.5, 14.3, 12.5, 10.5, 8.2])).tolist()
        gg   = '60-40'
        
    
    var1, leg, histlimit = varChoice(var)
    left, right, peak = posParametersInt(i,gas)
    
    if gas == 70:
        tf_fe55 = ROOT.TFile('../reco_run02274_to_run02280.root')    
    else:
        tf_fe55 = ROOT.TFile('../reco_run02160_to_run02165.root')
        
    tree = tf_fe55.Get('Events')
    
    c = ROOT.TCanvas('','',800,600)

    if gas == 70:
        if i == 6:
            hist = ROOT.TH1F('hist','%.2f cm between Source and GEM' % (dist[i]),35,0,histlimit)
        else:
            hist = ROOT.TH1F('hist','%.2f cm between Source and GEM' % (dist[i]),100,0,histlimit)
        
    else:
        if var == 'intsize':
            hist = ROOT.TH1F('hist','%.2f cm between Source and GEM' % (dist[i]),100,0,histlimit)
        else:
            hist = ROOT.TH1F('hist','%.2f cm between Source and GEM' % (dist[i]),30,0,histlimit)
        
    hist.Sumw2()
    
    if var == 'integral':
        cut_base = 'cl_iteration==2 && run=={r}'.format(r=run[i])
    else:
        le,ri = posParametersCutInt(i,gas)
        cut_base = 'cl_iteration==2 && run=={r} && cl_integral>={left} && cl_integral<={right}'.format(r=run[i],left=le, right=ri)
    cut = "{base} && TMath::Hypot(cl_xmean-1024,(cl_ymean-1024)*1.2)<{r_max}".format(base=cut_base,r_max=400)

    tree.Draw("{var}>>hist".format(var=var1),cut)
    hist.SetFillStyle(3005)
    mean = hist.GetMean()
    rms  = hist.GetRMS()
    
    # add Polya Fuction
    func='gauss'
    
    if func == 'gauss':
        print("Using Gauss fit")
        f = ROOT.TF1('f','gaus')
        f.SetParameter(1,mean);
        f.SetParLimits(1,mean-3*rms,mean+3*rms);
        f.SetParameter(2,rms);
        #f.SetParLimits(2,300,600);
        if var == 'integral':
            fitRe = hist.Fit(f,'R','',left,right)
        elif var == 'length':
            left1, right1 = posParametersLen(i,gas)
            fitRe       = hist.Fit(f,'R','',left1,right1)
        elif var == 'slimness':
            fitRe       = hist.Fit(f,'R','',0.5,1.2)
        elif var == 'intsize':
            fitRe       = hist.Fit(f,'R','',2,7)    
        else:
            fitRe = hist.Fit(f,'S')
        rMean  = f.GetParameter(1)
        rSigma = f.GetParameter(2)
        spotNumber = f.Integral(0,histlimit)/hist.GetBinWidth(1)
        
    else:
        print("Using Polya fit")
             
        #TF1 *myPolyaE = new TF1("myPolyaE","[0]*ROOT::Math::negative_binomial_pdf(x,[1],[2])+[3]*exp(-x/[4])",0,10000);
        #myPolyaE->SetParameters(10000,0.4,1000,1e3,1e2);
        
        myPolyaE = ROOT.TF1("myPolyaE","[0]*ROOT::Math::negative_binomial_pdf(x,[1],[2])+[3]*exp(-x/[4])",0,10000);
        if var == 'integral':
            if gas == 70:
                myPolyaE.SetParameters(10000,0.4,1000,1e3,1e2);
            else:
                myPolyaE.SetParameters(10000,0.4,peak,1e3,1e2);
        elif var == 'intsize':
            myPolyaE.SetParameters(30000,0.1,3,1e5,1e3);
        else:
            print('here')
            myPolyaE.SetParameters(10000,0.4,1,1e3,1e2);
       
        
        if var == 'intsize':
            hist.Fit(myPolyaE ,"R", "", 0.5,10)
        else:
            hist.Fit(myPolyaE ,"R", "", left,right)
            
        
        myPolya = ROOT.TF1("myPolya","[0]*ROOT::Math::negative_binomial_pdf(x,[1],[2])",0,10000);
        myPolya.SetParameters(myPolyaE.GetParameter(0),myPolyaE.GetParameter(1),myPolyaE.GetParameter(2));
        
        myPolya.Draw('hist same')
        myPolya.SetLineColor(ROOT.kGray+1)
        
        rMean = ((1-myPolya.GetParameter(1))*myPolya.GetParameter(2))/myPolya.GetParameter(1)
        rSigma = np.sqrt(((1-myPolya.GetParameter(1))*myPolya.GetParameter(2))/(myPolya.GetParameter(1)*myPolya.GetParameter(1)))
        rRMS = rSigma/rMean
        spotNumber = myPolya.Integral(0,histlimit)/hist.GetBinWidth(1)
        
        print("Mean: %.2f" % rMean)
        print("Sigma: %.2f" % rSigma)
        print("RMS: %.2f" % rRMS)
        print("Spot Number: %.2f" % spotNumber)
        
    hist.Draw('hist same')
    
    ROOT.gPad.Update()
    #h.GetXaxis().SetRangeUser(0,0.25)
    hist.GetYaxis().SetTitle('Counts')
    hist.GetXaxis().SetTitle('{leg}'.format(leg=leg))
    #'Position %d' % (pos[i])
    c.SetGrid()
    c.Draw()
    c.SaveAs("{plotdir}/hist_{var}_pos{pos}_{gg}_{func}.pdf".format(plotdir=plotdir, var=var, pos=pos[i], gg=gg, func=func))
    hist.Reset()
    
    return rMean,rSigma,dist[i],leg,spotNumber

def posParametersCutInt(pos,gas):
    
    if gas == 60:
        if pos == 0:
            left  = 1000
            right = 2400

        elif pos == 1:
            left  = 1600
            right = 3200

        elif pos == 2:
            left = 2000
            right = 3400

        elif pos == 3:
            left = 2400
            right = 3600

        elif pos == 4:
            left = 2400
            right = 3800

        elif pos == 5:
            left = 2400
            right = 4000

        else:
            exit()
    else:
        if pos == 0:
            left = 900
            right = 2300

        elif pos == 1:
            left = 1200
            right = 2700

        elif pos == 2:
            left = 1600
            right = 2900

        elif pos == 3:
            left = 1600
            right = 3100

        elif pos == 4:
            left = 1800
            right = 3200

        elif pos == 5:
            left = 1900
            right = 3300

        elif pos == 6:
            left = 2000
            right = 3300
    
    return left, right


def posParametersInt(pos,gas):
    peak = 0
    
    if gas == 60:
        if pos == 0:
            left  = 800
            right = 4000
            peak  = 1000

        elif pos == 1:
            left  = 1000
            right = 4200
            peak  = 1200

        elif pos == 2:
            left = 1600
            right = 4200
            peak  = 1500

        elif pos == 3:
            left = 1000
            right = 6000
            peak  = 2200

        elif pos == 4:
            left = 500
            right = 6000
            peak  = 2300

        elif pos == 5:
            left = 500
            right = 6000
            peak  = 2400

        else:
            exit()
    else:
        if pos == 0:
            left = 800
            right = 3000

        elif pos == 1:
            left = 1200
            right = 3400

        elif pos == 2:
            left = 1000#1500
            right = 4000#3500

        elif pos == 3:
            left = 1600
            right = 4000

        elif pos == 4:
            left = 1600
            right = 4000

        elif pos == 5:
            left = 1700
            right = 4000

        elif pos == 6:
            left = 1900
            right = 4300
    
    return left, right, peak

def posParametersLen(pos,gas):
    
    if gas == 60:
        if pos == 0:
            left = 2
            right = 4

        elif pos == 1:
            left = 3
            right = 5

        elif pos == 2:
            left = 3
            right = 5

        elif pos == 3:
            left = 3
            right = 5
            
        elif pos == 4:
            left = 3.5
            right = 5.5

        elif pos == 5:
            left = 3.5
            right = 5.5

        else:
            exit()
    else:
        if pos == 0:
            left = 2
            right = 4

        elif pos == 1:
            left = 2.5
            right = 4.5

        elif pos == 2:
            left = 3
            right = 5

        elif pos == 3:
            left = 3
            right = 5
            
        elif pos == 4:
            left = 3.5
            right = 5.5

        elif pos == 5:
            left = 3.5
            right = 6

        elif pos == 6:
            left = 3.5
            right = 6
    
    return left,right
    

def plotHist2D(plotdir,v1='integral',v2='slimness',i=0):
    import numpy as np
    
    ROOT.gStyle.SetOptFit(1011)
    
    gas=70
    #i=5
    
    if gas == 70:
        pos  = [0, 1, 2, 3, 4, 5, 6]
        run  = [2274, 2275, 2276, 2277, 2278, 2279, 2280]
        dist = (23-np.array([19.5, 16.5, 14.3, 12.5, 10.5, 8.2, 6.2])).tolist()
        gg   = '70-30'
    else:
        pos  = [0, 1, 2, 3, 4, 5]
        run  = [2160, 2161, 2162, 2163, 2164, 2165]
        dist = (23-np.array([19.5, 16.5, 14.3, 12.5, 10.5, 8.2])).tolist()
        gg   = '60-40'
        
    if gas == 70:
        tf_fe55 = ROOT.TFile('../reco_run02274_to_run02280.root')    
    else:
        tf_fe55 = ROOT.TFile('../reco_run02160_to_run02165.root')        
    
    vr1,legy,histlimity = varChoice(v1)
    vr2,legx,histlimitx = varChoice(v2)
        
    var1 = vr1+":"+vr2
        
    tree = tf_fe55.Get('Events')
    
    c = ROOT.TCanvas('','',800,800)

    hist = ROOT.TH2F('hist2D','%.2f cm between Source and GEM' % (dist[i]),1000,0,histlimity,1000,0,histlimitx)
    hist.Sumw2()

    cut_base = 'cl_iteration==2 && run=={r}'.format(r=run[i])
    cut = "{base} && TMath::Hypot(cl_xmean-1024,(cl_ymean-1024)*1.2)<{r_max}".format(base=cut_base,r_max=400)

    tree.Draw("{var}>>hist".format(var=var1),cut)
    hist.SetFillStyle(3005)
    
    #hist.Draw('hist sames') 
    ROOT.gPad.Update()
    #hist.GetYaxis().SetRangeUser(0,10000)
    hist.GetYaxis().SetTitle('{leg}'.format(leg=legy))
    hist.GetXaxis().SetTitle('{leg}'.format(leg=legx))
    c.SetGrid()
    c.Draw()
    c.SaveAs("{plotdir}/hist_{var}_vs_{var2}_pos{pos}_{gg}.pdf".format(plotdir=plotdir,var=v1,var2=v2,pos=pos[i],gg=gg))
    hist.Reset()
    
    
if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('', '--make'   , type='string'       , default='tworuns' , help='run simple plots (options = tworuns, evsdist, pmtvsz, cluvspmt, cluvsz, multiplicity, hist1d, hist2d)')
    parser.add_option('', '--outdir' , type='string'       , default='./'      , help='output directory with directory structure and plots')
    parser.add_option('', '--var' , type='string'       , default='integral'      , help='variable to plot the histogram')
    parser.add_option('', '--pos' , type='int'       , default=0      , help='position of the iron source')
    parser.add_option('', '--var2' , type='string'       , default='slimness'      , help='variable2 to plot the histogram 2D')
    parser.add_option('', '--gas' , type='int'       , default='60'      , help='gas mixture 60 or 70')
   
    (options, args) = parser.parse_args()

    ## make the output directory first
    os.system('mkdir -p {od}'.format(od=options.outdir))
    
    if options.make in ['all','multiplicity']:
        plotNClusters()

    if options.make in ['all','tworuns']:
        histograms = fillSpectra()
        odir = options.outdir+'/clusters/'
        os.system('mkdir -p {od}'.format(od=odir))
        drawSpectra(histograms,odir)
        os.system('cp index.php {od}'.format(od=odir))
    
    if options.make in ['all','evsdist']:
        plotEnergyVsDistance(options.outdir)

    if options.make in ['all','pmtvsz']:
        plotPMTEnergyVsPosition(options.outdir,options.gas)
    
    if options.make in ['all','pmt2gas']:
        plotPMT2gas(options.outdir)

    if options.make in ['all','cluvspmt']:
        plotCameraPMTCorr(options.outdir)
        
    if options.make in ['all','cluvsz']:
        plotCameraEnergyVsPosition(options.outdir, options.var, options.gas)
    
    if options.make in ['all','hist1d']:        
        plotHistFit(options.outdir, options.var, options.pos, options.gas)
        
    if options.make in ['all','hist2d']:
        plotHist2D(options.outdir, options.var, options.var2, options.pos)
    
    if options.make in ['all','gascomp']:
        plotCameraEnergyVsPositionGas(options.outdir, options.var)

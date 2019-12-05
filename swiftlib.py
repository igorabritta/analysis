#!/usr/bin/env python

def swift_root_file(sel, run):
    BASE_URL  = "https://swift.cloud.infn.it:8080/v1/AUTH_1e60fe39fba04701aa5ffc0b97871ed8/Cygnus/"
    file_root = (sel+'/histograms_Run%05d.root' % run)
    return BASE_URL+file_root

def reporthook(blocknum, blocksize, totalsize):
    import sys
    readsofar = blocknum * blocksize
    if totalsize > 0:
        percent = readsofar * 1e2 / totalsize
        s = "\r%5.1f%% %*d / %d" % (
            percent, len(str(totalsize)), readsofar, totalsize)
        sys.stderr.write(s)
        if readsofar >= totalsize: # near the end
            sys.stderr.write("\n")
    else: # total size is unknown
        sys.stderr.write("read %d\n" % (readsofar,))

def swift_download_root_file(url,run):
    import ROOT
    import os
    from urllib import urlretrieve
    tmpname = ("/tmp/histograms_Run%05d.root" % run)
    urlretrieve(url, tmpname, reporthook)
    return tmpname 

def swift_read_root_file(tmpname):
    import ROOT
    f  = ROOT.TFile.Open(tmpname);
    return f

def swift_rm_root_file(tmpname):
    import os
    os.remove(tmpname)
    #print("tmp file removed")

def checkfiletmp(run):
    import os.path
    return os.path.isfile("/tmp/histograms_Run%05d.root" % run)


def root_TH2_name(root_file):
    pic = []
    wfm = []
    for i,e in enumerate(root_file.GetListOfKeys()):
        che = e.GetName()
        if ('pic_run' in str(che)):
            pic.append(che)
        elif ('wfm_run' in str(che)):
            wfm.append(che)
    return pic, wfm
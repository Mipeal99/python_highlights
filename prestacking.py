# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits
import shutil
import os
import matplotlib.pyplot as plt

reduced_folder=r'C:\Users\Usuario\Desktop\reduced\NGC2392'
science_filter_list=['SDSS_r','SDSS_g','SDSS_i','Hb','Ha','HeI','SII','OIII'] #list of all the filters where you observed
for filter in science_filter_list: #open the reduced exposures
    print('-',filter,' band')
    
    for filename in os.listdir(os.path.join(reduced_folder,filter)):
        if os.path.isdir(os.path.join(reduced_folder,filter,filename)):
            continue
        
        hdu=fits.open(os.path.join(reduced_folder,filter,filename))
        print('opening',filename)
        bad_data=hdu[0].data #get the data
        bad_data[:,3700:]=0 #just delete the bad pixels (alternative to cropping that doesn't change the resolution)
        places=bad_data>=0 #also delete the bad pixels with negative counts
        good_data=bad_data*(places*1)
        #sanity check plot, can (and maybe should) comment because it opens a lot of figures
        # plt.figure()
        # img=plt.imshow(good_data,origin='lower', vmin=np.nanpercentile(good_data, 1), vmax=np.nanpercentile(good_data, 99), cmap='inferno') 
        # plt.title(filter+','+ filename)
        # plt.colorbar()
        # plt.show()
        if not os.path.isdir(os.path.join(reduced_folder,filter,'fixed')):
            print('creating folder for ', filter)
            os.makedirs(os.path.join(reduced_folder,filter,'fixed'))
        fits.writeto(os.path.join(reduced_folder,filter,'fixed',filename+'FIXED.fits'), good_data, hdu[0].header, overwrite=True) #save the good data to the desired location
        print('saved ',filename)
        
    if filter=='SDSS_r': #adding the 3 SDSS-r exposures, kind of scuffed and unoptimized but it works
        hdu1=fits.open(os.path.join(reduced_folder,filter,'fixed','ABC0007 (2).fitsFIXED.fits'))
        dat1=hdu1[0].data
        hdu2=fits.open(os.path.join(reduced_folder,filter,'fixed','ABC0008 (2).fitsFIXED.fits'))
        dat2=hdu2[0].data
        hdu3=fits.open(os.path.join(reduced_folder,filter,'fixed','ABC0009 (2).fitsFIXED.fits'))
        dat3=hdu3[0].data
        stack=dat1+dat2+dat3
        fits.writeto(os.path.join(r'C:\Users\Usuario\Desktop\reduced\stack','SDSS_r.fits'), stack, hdu1[0].header, overwrite=True)
        print('Stack saved!')
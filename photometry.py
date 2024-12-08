
# Photometry Function

def photometry(section,sources,r_ap,r_in,r_out,exptime,zeropoint,color):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from photutils.aperture import CircularAnnulus as CircAn, CircularAperture as CircAp, ApertureStats as ApStats, aperture_photometry as ApPho

    # Creating the apertures and establishing the star positions
    for col in sources.colnames:
        if col not in ('id', 'npix'):
            sources[col].info.format = '%.2f'

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircAp(positions, r=r_ap)
    annulus_aperture = CircAn(positions, r_in = r_in, r_out = r_out)

    plt.figure()
    plt.imshow(section, origin = 'lower', norm = LogNorm(), cmap = 'Greys' )
    apertures.plot(color = 'green', lw = 1.5, alpha = 0.5)
    annulus_aperture.plot(color = color, lw = 1.5, alpha = 0.5)
    
    # Background
    bkgstats = ApStats(section, annulus_aperture)
    apstats = ApStats(section, apertures)
    bkg_mean = bkgstats.mean
    aperture_area = apertures.area_overlap(section)
    total_bkg = bkg_mean * aperture_area
    star_data = ApPho(section, apertures)
    star_data['total_bkg'] = total_bkg
    for col in star_data.colnames:
        star_data[col].info.format = '%.8g'

    # Converting to magnitudes & calculating uncertainty
    mags = []
    fluxes = []
    for line in star_data:
        flux = (abs(line[3]-line[4])/exptime)
        fluxes.append(flux)
        mags.append((zeropoint-2.5*np.log10(flux)))
    error = apstats.std/(np.sqrt(len(apstats)))
    unc_mag = (((2.5/np.array(fluxes)))*error)
    sul = np.mean(mags) + 3*np.std(mags)
    star_data['magnitudes'] = mags
    star_data.pprint()
    print('')
    return mags, unc_mag, sul
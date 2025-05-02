## 
## Copyright (C) 2024 by
## Chiara Magosso
## 
## This work is licensed under a  
## 
## Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
## ( http://creativecommons.org/licenses/by-nc-sa/4.0/ )
## 
## Please contact chiara.magosso@polito.it for information.
##

"""
====================================================================
Chiara Magosso - Automated analysis of lamellar BCP images form SEM
====================================================================
"""

############################################################
# import
############################################################
from math import sqrt
from scipy.misc import derivative
import tifffile
import numpy as np
import matplotlib.pyplot as plt
import cv2 
import skimage.color
import skimage.io
from scipy import fftpack
import subprocess
from matplotlib.colors import LogNorm
import scipy.interpolate
import pandas as pd
import fingerprint_enhancer	
import os
import json
from PIL import Image
import shutil
import matplotlib
import gc
import codecs
import scipy.optimize as spopt
from scipy.spatial import distance
import xlsxwriter

np.set_printoptions(threshold=np.inf)

############################################################
# Metadata extractor 
############################################################
def metadata(immagine_sem, indirizzo, nome):
    metadata = open(f'{indirizzo}analysis/metadati/{nome}_metadata.txt', 'w')
    with tifffile.TiffFile(immagine_sem) as tif:
        for page in tif.pages:
            for tag in page.tags:
                tag_name, tag_value = tag.name, tag.value
                if  isinstance(tag_value, dict):      #tag_name in ['CZ_SEM', 'FEI_HELIOS']: 
                    metadata.write(f'{tag_name} \n')
                    # different SEMs               
                    if tag_name == 'CZ_SEM':
                        sem = 'CZ_SEM'
                        InvRisoluzione = tag_value['ap_image_pixel_size'][1] # nm/pix  
                        #magnification = np.multiply(float(tag_value['ap_mag'][1].split(' ', 1)[0]), 1000) #Kx
                        brightness = tag_value['ap_brightness'][1]
                        contrast = tag_value['ap_contrast'][1]
                    if tag_name == 'FEI_HELIOS':
                        sem = 'FEI_HELIOS'
                        InvRisoluzione = np.multiply(tag_value['Scan']['PixelWidth'], 10**9) # nm/pix  
                        #magnification = 30000 ## da fare 
                        #brightness = np.nan
                        #contrast = np.nan
                        try: 
                            brightness = tag_value['ETD']['Brightness']
                            contrast = tag_value['ETD']['Contrast']
                        except KeyError:
                            brightness = tag_value['LVSED']['Brightness']
                            contrast = tag_value['LVSED']['Contrast']
  
                    for key, value in tag_value.items():
                        metadata.write(f'{key} \t {value} \n')
                else:
                    metadata.write(f"{tag_name} \t {codecs.decode(tag_value, 'unicode_escape') if isinstance(tag_value, str) else tag_value} \n")
                    if tag_name == 'ImageDescription':
                        process_param = eval(tag_value)
                        for key, value in zip(process_param.keys(), process_param.values()):
                            if isinstance(value, list):
                                process_param[key]=value[0]
                        df_process_param = pd.DataFrame([process_param])
                    
                    if tag_name == 'Artist':
                        y = json.loads(tag_value)
                        sem = 'Hitachi_SU4800'
                        InvRisoluzione = float(y['PixelSize'][0])                      
                        #magnification = tag_value['Magnification']
                        brightness = np.NaN
                        contrast = np.NaN
                        
    metadata.close()
    return sem, InvRisoluzione, brightness, contrast, df_process_param

############################################################
# Read and plot the image
############################################################
def original_image(immagine_sem):

    im = skimage.color.rgb2gray(skimage.io.imread(immagine_sem)) 
    if np.amax(im) > 1: 
        #im = im/255
        im = (im-np.min(im))/(np.max(im)-np.min(im))
    plt.figure()
    plt.imshow(im, plt.cm.gray)
    plt.title('Original image')
    return im

############################################################
# Crop
############################################################
def crop(im):
    r, c = im.shape  #rows and columns of original image
    bianco_95 = c * 0.95 # mostly white data bar
    bianco_5 = c * 0.05 # mostly black data bar
    riga = im.sum(axis=1)
    legenda = None
    if np.amin(riga)<=bianco_5:
        legenda=np.where(riga<bianco_5)
    elif np.amax(riga)>=bianco_95:
        legenda=np.where(riga>bianco_95)

    if legenda:
        crop_im = np.delete(im, np.s_[np.amin(legenda):r], 0) # create new image without data bar
    else:
        crop_im = im

    plt.figure()
    plt.imshow(crop_im, plt.cm.gray)
    plt.title('Crop image')
    return crop_im

############################################################
# Show the results
############################################################
def plot_spectrum(im_fft):
    plt.imshow(np.fft.fftshift(np.abs(im_fft)), norm=LogNorm(vmin=5))
    plt.colorbar()

############################################################
# Compute the 2d FFT of the input image
############################################################
def fft(crop_im):
    im_fft = fftpack.fft2(crop_im)
    plt.figure()
    plot_spectrum(im_fft)
    plt.title('Fourier transform')
    return im_fft

############################################################
# Radial Mean
############################################################
def Radial_Mean(im_fft, InvRisoluzione):
    Fsh = np.fft.fftshift(im_fft)
    Y1, X1 = Fsh.shape
    N_pixel = X1 * Y1
    Z = np.log(1+abs(Fsh))
    L1 = X1 * InvRisoluzione
    L2 = Y1 * InvRisoluzione
    LFFT1 = 1/L1
    LFFT2 = 1/L2
    FREQ1 = LFFT1 * X1
    FREQ2 = LFFT2 * Y1
    x = np.linspace(-FREQ1/2, FREQ1/2, Y1)
    y = np.linspace(-FREQ1/2, FREQ1/2, X1)
    X, Y = np.meshgrid(y,x)
    radii = np.linspace(0, FREQ1/2, int(X1/2-1))
    valor_medio = []

    for r in radii:
        r2 = r + LFFT1
        pixel_radii = np.add(np.multiply(X,X), np.multiply(Y,Y))
        filtro1 = pixel_radii >= r**2 #(r_max - sigma)**2
        filtro2 = pixel_radii <= r2**2 #(r_max + sigma)**2
        filtro = np.multiply(filtro1, filtro2)
        filtered_Z = np.multiply(Z, filtro)
        integrale = np.sum(filtered_Z)
        pixel_nulli = np.nonzero(filtered_Z==0)
        valor_medio.append(integrale/(N_pixel - len(pixel_nulli[0])))

    valor_medio = np.array(valor_medio)
    r1 = np.transpose(radii)
    plt.figure()
    plt.plot(r1,valor_medio,'.')
    plt.title('Radial Mean')
    plt.xlabel("Ferq [nm^(-1)]")
    plt.ylabel("Intensity")
    return valor_medio, r1

############################################################
# sigma peak
############################################################
def sigma(der_seconde, x_picco, y_picco, spline, freq, intensity): 
    #find the mins before and after the peak
    index_picco = der_seconde[der_seconde['x'] == x_picco].index.values 
    min_SX_picco = der_seconde.iloc[index_picco-1]
    x_min_SX_picco = min_SX_picco.iloc[0]['x']
    y_min_SX_picco = min_SX_picco.iloc[0]['y']
    plt.plot(x_min_SX_picco, y_min_SX_picco, color='yellow', marker='o')
    try:
        min_DX_picco = der_seconde.iloc[index_picco+1]
        x_min_DX_picco = min_DX_picco.iloc[0]['x']
        y_min_DX_picco = min_DX_picco.iloc[0]['y']
        plt.plot(x_min_DX_picco, y_min_DX_picco, color='yellow', marker='o')
    except IndexError:
        print('DECRESCE a inf ultimo')
        x_min_DX_picco = freq[-1]
        y_min_DX_picco = intensity[-1]
        plt.plot(x_min_DX_picco, y_min_DX_picco, color='yellow', marker='o')

    # find the straight line that passes through the two mins of the peak
    # Define the known points
    x = [x_min_SX_picco, x_min_DX_picco]
    y = [y_min_SX_picco, y_min_DX_picco]

    # Calculate the coefficients --> line y = m*x + c
    coefficients = np.polyfit(x, y, 1)

    # Let's compute the values of the line
    polynomial = np.poly1d(coefficients)
    y_axis = polynomial(x_picco)

    #plot the points and the line
    plt.plot(x_picco, y_axis, color='yellow', marker='o')

    #find the half height 
    delta = (y_picco - y_axis)/2
    HM = delta + y_axis
    plt.plot(x_picco, HM, color='pink', marker='o')

    #find the two points where the line y = HM intersects the spline
    coefficients[0] = 0
    coefficients[1] = HM
    line = np.poly1d(coefficients)
    f = lambda x: spline(x) - line(x) # difference function, its zero marks the intersection
    
    # find root via bisection
    left_int = spopt.bisect(f, a = x_min_SX_picco, b = x_picco)
    right_int = spopt.bisect(f, a = x_picco, b = x_min_DX_picco)
    plt.plot(left_int, spline(left_int), color='pink', marker='o')
    plt.plot(right_int,  spline(right_int), color='pink', marker='o')

    #FWHM
    a = (right_int,  spline(right_int))
    b = (left_int, spline(left_int))
    #test = distance.euclidean(a, b)
    FWHM = right_int - left_int 

    #sigma peak
    sigma_xc = FWHM/(2*np.sqrt(2*np.log(2)))
    sigma_picco = np.sqrt(((1/(x_picco**2))*sigma_xc)**2)
    print('l0:',1/x_picco, sigma_picco)

    return sigma_picco

############################################################
# l0 e d
############################################################
def interpolate(valor_medio, r1):
    spl = scipy.interpolate.UnivariateSpline(r1, valor_medio, k=4, s=1)
    plt.figure()
    plt.plot(r1,valor_medio,'.')
    plt.plot(r1, spl(r1), 'g', lw=3)
    plt.title('Radial Mean interpolate')
    plt.xlabel("Ferq [nm^(-1)]")
    plt.ylabel("Intensity")
    plt.grid()

    cr_pts = spl.derivative().roots()
    cr_vals = spl(cr_pts)
    sd_spl = spl.derivative().derivative()
    sd_vals = sd_spl(cr_pts)
    pairs = pd.DataFrame({
        'x': cr_pts,
        'y': cr_vals,
        'ysecder': sd_vals
    })

    minimo = pairs.query('ysecder > 0')
    pairs = pairs.query('ysecder < 0')
    derseconde = pd.concat([pairs, minimo]).sort_values(by=['x'])

    while pairs.iloc[0]['x']<0.01 : # l0 largest possible 100nm
            pairs = pairs.drop(pairs.index[0])
    pairs = pairs.sort_values(by=['y'], ascending=False)
    l0 = 1/pairs.iloc[0]['x']
    sigma_l0 = sigma(derseconde, pairs.iloc[0]['x'], pairs.iloc[0]['y'], spl, r1, valor_medio)
    plt.plot(pairs.iloc[0]['x'], pairs.iloc[0]['y'], color='red', marker='o')

    if len(pairs)>1 :
        for i in range(1, len(pairs)):
            if pairs.iloc[i]['x']>1.8*pairs.iloc[0]['x'] and pairs.iloc[i]['x']<2.2*pairs.iloc[0]['x']: #and len(pairs)>1:  
                d = 1/pairs.iloc[i]['x']
                plt.plot(pairs.iloc[i]['x'], pairs.iloc[i]['y'], color='red', marker='o')              
                sigma_d = np.NaN
                #sigma_d = sigma(derseconde, pairs.iloc[i]['x'], pairs.iloc[i]['y'], spl)
            else:
                d=l0/2
                sigma_d = np.NaN
                #sigma_d = np.sqrt((0.5 * sigma_l0)**2)
    else :  
        d=l0/2  
        sigma_d = np.NaN
        #sigma_d = np.sqrt((0.5 * sigma_l0)**2)   

    return l0, sigma_l0, d, sigma_d

############################################################
# remove black edges of Fingerprint
############################################################
def trim(frame):
    #crop top
    if not np.sum(frame[0]):
        return trim(frame[1:])
    #crop bottom
    elif not np.sum(frame[-1]):
        return trim(frame[:-2])
    #crop left
    elif not np.sum(frame[:,0]):
        return trim(frame[:,1:]) 
    #crop right
    elif not np.sum(frame[:,-1]):
        return trim(frame[:,:-2])    
    return frame

############################################################
# resize image
############################################################
def image_resize(img, scale_percent):
                    
    print('Original Dimensions : ',img.shape)
    
    #scale_percent percent of original size
    width = int(img.shape[1] * scale_percent / 100)
    height = int(img.shape[0] * scale_percent / 100)
    dim = (width, height)
    
    # resize image
    resized = cv2.resize(img, dim, interpolation = cv2.INTER_AREA)
    
    print('Resized Dimensions : ',resized.shape)
    return resized

#############################################
#Find edges \ developed by Stefano Carignano - stefano.carignano@bsc.es
#############################################

def find_squares(img):

    #img = cv2.imread('shapes.jpg')
    gray = img # cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    ret,thresh = cv2.threshold(gray,50,255,0)
    contours,hierarchy = cv2.findContours(thresh, 1, 2)
    print("Number of contours detected:", len(contours))

    for cnt in contours:
        x1,y1 = cnt[0][0]
        approx = cv2.approxPolyDP(cnt, 0.01*cv2.arcLength(cnt, True), True)
        if len(approx) == 4:
            x, y, w, h = cv2.boundingRect(cnt)
            ratio = float(w)/h
            if ratio >= 0.9 and ratio <= 1.1:
                img = cv2.drawContours(img, [cnt], -1, (0,255,255), 3)
                cv2.putText(img, 'Square', (x1, y1), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 0), 2)
            else:
                cv2.putText(img, 'Rectangle', (x1, y1), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 255, 0), 2)
                img = cv2.drawContours(img, [cnt], -1, (0,255,0), 3)

    cv2.imshow("Shapes", img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

#############################################
#PUF crop \ developed by Stefano Carignano & Chiara Magosso 
#############################################
 
def PUF_crop(in_fig, fancy=False):
    #print(in_fig.shape[1])
    if fancy:
        squares = find_squares(in_fig)
        cv2.drawContours(in_fig, squares, -1, (0, 255, 0), 3)
        plt.show()
        #raise NotImplementedError
    elif in_fig.shape[1] == 2048:
        return in_fig[500:-500,500:-500]
    elif in_fig.shape[1] == 1024:
        #plt.imshow(in_fig) ; plt.show()
        #plt.imshow(in_fig[250:-250,250:-250]) ; plt.show()
        return in_fig[240:-240,240:-240]
    else:
        print('Add case in function PUF_crop')
        sys.exit(1)
    

############################################################
# MAIN
############################################################
def main(): 
    indirizzo ="C:/Users/machinelearning/Desktop/gui_metadata/image_with_process_metadata/" #modify accordingly
    cartelle = ['ADAblock', 'crop', 'PUF_matching_crop','Fingerprint-Enhancement', 'metadati', 'radial_mean', 'spettri', 'image_with_all_metadata']
    partial_database_process = pd.DataFrame()
    partial_database_python = pd.DataFrame()
    database_process = pd.DataFrame()
    database_python = pd.DataFrame()
    database_ADAblock = pd.DataFrame()
    database = pd.DataFrame()
    for car in cartelle:
        os.makedirs(os.path.join(indirizzo, 'analysis', car), exist_ok=True)
        os.makedirs(os.path.join(indirizzo, 'not_analyzed', 'image_with_some_metadata'), exist_ok=True)
        os.makedirs(os.path.join(indirizzo, 'done'), exist_ok=True)
    for file in os.listdir(indirizzo):
        if file.endswith(".tif"):

            nome = os.path.splitext(file)[0]
            print('\n'f'{nome}')
            immagine_sem = (f'{indirizzo}{file}')
            
            print('Metadata extraction') 
            sem, InvRisoluzione, brightness, contrast, df_process_param = metadata(immagine_sem, indirizzo, nome) #nm/pix

            colonne = list(df_process_param.columns)
            a = 'Molecular_weight_A_g_over_mol'
            b = 'Molecular_weight_B_g_over_mol'
            c = 'Molecular_weight_C_g_over_mol'
            if a and b and c in colonne: 
                df_process_param = df_process_param.rename(columns={'Molecular_weight_A_g_over_mol':'Molecular_weight_BCP_A_g_over_mol','Molecular_weight_B_g_over_mol':'Molecular_weight_BCP_B_g_over_mol','Molecular_weight_C_g_over_mol':'Molecular_weight_BCP_C_g_over_mol'})

            #if InvRisoluzione > 15 or InvRisoluzione <= 1.1:
            #    shutil.move(os.path.join(indirizzo, file), f'{indirizzo}not_analyzed')
            
            im = original_image(immagine_sem)
            
            print('SEM image Crop')
            crop_im = crop(im)
            plt.imsave(f'{indirizzo}analysis/crop/{nome}_crop.tif', crop_im, format='tiff', cmap='gray')
            
            #crop_im = PUF_crop(crop_im, fancy=False) # activate only if you are doing PUF matching
            #plt.imsave(f'{indirizzo}analysis/PUF_matching_crop/{nome}_crop.tif', crop_im, format='tiff', cmap='gray') # activate only if you are doing PUF matching

            print('Radial Mean')
            im_fft = fft(crop_im)
            plt.savefig(f'{indirizzo}analysis/spettri/{nome}_spectrum.png')
            
            valor_medio, r1 = Radial_Mean(im_fft, InvRisoluzione)
            
            workbook = xlsxwriter.Workbook(f'{indirizzo}analysis/radial_mean/{nome}_radial_mean.xlsx')
            worksheet = workbook.add_worksheet()
            
            # Start from the first cell.
            # Rows and columns are zero indexed.
            row = 0
            column = 0
            
            # iterating through content list
            for item in valor_medio :
            
                # write operation perform
                worksheet.write(row, column, item)
            
                # incrementing the value of row by one with each iterations
                row += 1
            row = 0    
            for item in r1 :
            
                # write operation perform
                worksheet.write(row, column+2, item)
                # incrementing the value of row by one with each iterations
                row += 1

            workbook.close()

            plt.savefig(f'{indirizzo}analysis/radial_mean/{nome}_radial_mean.png')
            
            l0, sigma_l0, d, sigma_d = interpolate(valor_medio, r1) #nm
            plt.savefig(f'{indirizzo}analysis/radial_mean/{nome}_radial_mean_interp.png')
            
            data_analysis = open(f'{indirizzo}analysis/metadati/{nome}_analysis.txt', 'w')
            data_analysis.write(f'{nome} \nInvRisoluzione = {InvRisoluzione:.3f} \nl0 = {l0:.3f} sl0 = {sigma_l0:.3f} \nd = {d:.3f} sd = {sigma_d:.3f}\n')
            data_analysis.close()
            
            python_param = {
                "Image": nome, 
                "SEM" : sem,
                #"magnification" : magnification,
                "brightness" : brightness,
                "contrast" : contrast,
                "InvRisoluzione" : InvRisoluzione,
                "l0" : l0,  
                "sl0" : sigma_l0,
                "d" : d, 
                "sd" : sigma_d
            }

            df_python_data=pd.DataFrame([python_param])
            
            print('Fingerprint binarization')
            print('wavelength : ', l0/InvRisoluzione)
            try: 
                out = fingerprint_enhancer.enhance_Fingerprint(crop_im, resize=False, ridge_segment_blksze=16, ridge_segment_thresh=0.1, gradient_sigma=1, block_sigma=7, orient_smooth_sigma=7,
                        ridge_freq_blksze=round(5*l0/InvRisoluzione), ridge_freq_windsze=5, min_wave_length=l0/InvRisoluzione-(0.05*l0/InvRisoluzione), max_wave_length=l0/InvRisoluzione+(0.05*l0/InvRisoluzione), kx=0.50, ky=0.50, angleInc=3.0, ridge_filter_thresh=-3)		# enhance the fingerprint image

            except IndexError:
                print('FP non funziona')

                image = Image.open(immagine_sem)

                py_analysis = df_python_data.to_dict('r')[0]
                py_analysis['brightness']=[py_analysis['brightness'], r'\one']
                py_analysis['contrast']=[py_analysis['contrast'], r'\one']
                py_analysis['InvRisoluzione']=[py_analysis['InvRisoluzione'], r'\nano\metre\pixel\tothe{-1}']
                py_analysis['l0']=[py_analysis['l0'], r'\nano\metre']
                py_analysis['sl0']=[py_analysis['sl0'], r'\nano\metre']
                py_analysis['d']=[py_analysis['d'], r'\nano\metre']
                py_analysis['sd']=[py_analysis['sd'], r'\nano\metre']

                py_analysis = json.dumps(py_analysis)
                image.save(f'{indirizzo}not_analyzed/image_with_some_metadata/{file}', tiffinfo=image.tag, software=py_analysis)
                image.close()
                
                print('Partial Databese creation')  
                              
                partial_database_python = partial_database_python.append(df_python_data, ignore_index = True)
                partial_database_process = partial_database_process.append(df_process_param, ignore_index = True)
                    
                partial_database = pd.concat([partial_database_process,partial_database_python], axis=1)

                partial_database.to_csv(f'{indirizzo}not_analyzed/image_with_some_metadata/partial_analisi_database.csv', sep=';')

                shutil.move(os.path.join(indirizzo, file), f'{indirizzo}not_analyzed')
                continue

            out = trim(out)

            #out = cv2.bitwise_not(out) # activate only if you are doing PUF
                
            nome_FP = f'{indirizzo}analysis/Fingerprint-Enhancement/{nome}_FE_InvRisol={InvRisoluzione:.3f}_d={d:.3f}_l0={l0:.3f}'
            plt.imsave(f'{nome_FP}.tif', out, format='tiff', cmap='gray')

            pixel_FP = open(f'{indirizzo}analysis/metadati/{nome}_pixel_FP.txt', 'w')
            pixel_FP.write(str(out))
            pixel_FP.close()

            print("ADAblock analysis")
            
            # The ImageJ The software can be downloaded at https://imagej.net/ - We recommend version ImageJ 1.51w
            # The medified script ADAblock.ijm can be found at https://github.com/ChiaraMagosso/ADAblock/blob/ijMacro_updates_only/ADAblock.ijm
            # ADAblock was originally published by Murphy, J. N., Harris, K. D. & Buriak, J. M. Automated Defect and Correlation Length Analysis of Block Copolymer Thin Film Nanopatterns. PLOS ONE 10, e0133088 (2015) - https://doi.org/10.1371/journal.pone.0133088. Follow some of the instructions in the "S1 Instructions - Use of ADAblock" of the above cited article. Specifically, the "Set up ImageJ" and "Sorting Through The Data - List of output files" sections are also valid for the modified version. 
            
            subprocess.run(f'C:/Users/machinelearning/Documents/Dottorato/INRiM/programmi/ImageJ/ImageJ.exe --console -macro C:/Users/machinelearning/Documents/Dottorato/INRiM/programmi/database_article/ADAblock_Chiara/ADAblock/ADAblock.ijm "{nome_FP} {InvRisoluzione} {indirizzo}analysis/ADAblock/"').returncode #modify accordingly
            
            ADAb = f'{indirizzo}analysis/ADAblock/{nome}_FE_InvRisol={InvRisoluzione:.3f}_d={d:.3f}_l0={l0:.3f}/outputTD_horizontal.xls'
            df1 = pd.DataFrame()
            if os.path.isfile(ADAb):
                df1 = pd.read_csv(ADAb, sep = '\t')
                column_list = ['Image_Title.String','nm_per_pixel','wfft_Period_nm','Width_final','Height_final', 'opa_Hermans',
                        'correlation_length_nm','correlation_length_linear_nm','LER_sigma_avg_nm',
                        'LWR_sigma_avg_nm','PTE','NTE','PT','NT','PDE',
                        'NDE','PD','ND','PJ3','NJ3','PJ4','NJ4',
                        'PJx','NJx','Ptot','Ntot','Total_Defects', 'Total_Area_nm',
                        'Defect_Density_nm','Defect_Density_um']
                for col in column_list:
                    if col not in df1.columns:
                        df1[col]=np.NaN

                df1 = df1[column_list] 
                df1 = df1.rename(columns={'Image_Title.String':'Image_Name','PTE':'Pos_Terminals-Edge','NTE':'Neg_Terminals-Edge','PT':'Pos_Terminals',
                                    'NT':'Neg_Terminals', 'PDE':'Pos_Dots-Edge','NDE':'Neg_Dots-Edge','PD':'Pos_Dots','ND':'Neg_Dots',
                                    'PJ3':'Pos_Junctions3','NJ3':'Neg_Junctions3','PJ4':'Pos_Junctions4', 'NJ4':'Neg_Junctions4','PJx':'Pos_Junctionsx',
                                    'NJx':'Neg_Junctionsx','Ptot':'Pos_Defects','Ntot':'Neg_Defects'})
                image = Image.open(immagine_sem)
                py_analysis = df_python_data.to_dict('r')[0]
                py_analysis['brightness']=[py_analysis['brightness'], r'\one']
                py_analysis['contrast']=[py_analysis['contrast'], r'\one']
                py_analysis['InvRisoluzione']=[py_analysis['InvRisoluzione'], r'\nano\metre\pixel\tothe{-1}']
                py_analysis['l0']=[py_analysis['l0'], r'\nano\metre']
                py_analysis['sl0']=[py_analysis['sl0'], r'\nano\metre']
                py_analysis['d']=[py_analysis['d'], r'\nano\metre']
                py_analysis['sd']=[py_analysis['sd'], r'\nano\metre']

                ADAblock_analysis = df1.to_dict('r')[0]
                ADAblock_analysis['nm_per_pixel']=[ADAblock_analysis['nm_per_pixel'], r'\nano\metre\pixel\tothe{-1}']
                ADAblock_analysis['wfft_Period_nm']=[ADAblock_analysis['wfft_Period_nm'], r'\nano\metre']
                ADAblock_analysis['Width_final']=[ADAblock_analysis['Width_final'], r'\pixel']
                ADAblock_analysis['Height_final']=[ADAblock_analysis['Height_final'], r'\pixel']
                ADAblock_analysis['opa_Hermans']=[ADAblock_analysis['opa_Hermans'], r'\one']
                ADAblock_analysis['correlation_length_nm']=[ADAblock_analysis['correlation_length_nm'], r'\nano\metre']
                ADAblock_analysis['correlation_length_linear_nm']=[ADAblock_analysis['correlation_length_linear_nm'], r'\nano\metre']
                ADAblock_analysis['LER_sigma_avg_nm']=[ADAblock_analysis['LER_sigma_avg_nm'], r'\nano\metre']
                ADAblock_analysis['LWR_sigma_avg_nm']=[ADAblock_analysis['LWR_sigma_avg_nm'], r'\nano\metre']
                ADAblock_analysis['Pos_Terminals-Edge']=[ADAblock_analysis['Pos_Terminals-Edge'], r'\one']
                ADAblock_analysis['Neg_Terminals-Edge']=[ADAblock_analysis['Neg_Terminals-Edge'], r'\one']
                ADAblock_analysis['Pos_Terminals']=[ADAblock_analysis['Pos_Terminals'], r'\one']
                ADAblock_analysis['Neg_Terminals']=[ADAblock_analysis['Neg_Terminals'], r'\one']
                ADAblock_analysis['Pos_Dots-Edge']=[ADAblock_analysis['Pos_Dots-Edge'], r'\one']
                ADAblock_analysis['Neg_Dots-Edge']=[ADAblock_analysis['Neg_Dots-Edge'], r'\one']
                ADAblock_analysis['Pos_Dots']=[ADAblock_analysis['Pos_Dots'], r'\one']
                ADAblock_analysis['Neg_Dots']=[ADAblock_analysis['Neg_Dots'], r'\one']
                ADAblock_analysis['Pos_Junctions3']=[ADAblock_analysis['Pos_Junctions3'], r'\one']
                ADAblock_analysis['Neg_Junctions3']=[ADAblock_analysis['Neg_Junctions3'], r'\one']
                ADAblock_analysis['Pos_Junctions4']=[ADAblock_analysis['Pos_Junctions4'], r'\one']
                ADAblock_analysis['Neg_Junctions4']=[ADAblock_analysis['Neg_Junctions4'], r'\one']
                ADAblock_analysis['Pos_Junctionsx']=[ADAblock_analysis['Pos_Junctionsx'], r'\one']
                ADAblock_analysis['Neg_Junctionsx']=[ADAblock_analysis['Neg_Junctionsx'], r'\one']
                ADAblock_analysis['Pos_Defects']=[ADAblock_analysis['Pos_Defects'], r'\one']
                ADAblock_analysis['Neg_Defects']=[ADAblock_analysis['Neg_Defects'], r'\one']
                ADAblock_analysis['Total_Defects']=[ADAblock_analysis['Total_Defects'], r'\one']
                ADAblock_analysis['Total_Area_nm']=[ADAblock_analysis['Total_Area_nm'], r'\nano\metre\tothe{2}']
                ADAblock_analysis['Defect_Density_nm']=[ADAblock_analysis['Defect_Density_nm'], r'\nano\metre\tothe{-2}']
                ADAblock_analysis['Defect_Density_um']=[ADAblock_analysis['Defect_Density_um'], r'\micro\metre\tothe{-2}']
                

                dic_complete = {**py_analysis, **ADAblock_analysis}
                dic_complete = json.dumps(dic_complete)
                image.save(f'{indirizzo}analysis/image_with_all_metadata/{file}', tiffinfo=image.tag, software=dic_complete)
                image.close()

            else:
                columns=['Image_Name', 'nm_per_pixel', 'wfft_Period_nm', 'Width_final', 'Height_final', 'opa_Hermans', 
                        'correlation_length_nm', 'correlation_length_linear_nm', 'LER_sigma_avg_nm',
                        'LWR_sigma_avg_nm', 'Pos_Terminals-Edge', 'Neg_Terminals-Edge', 'Pos_Terminals', 'Neg_Terminals', 'Pos_Dots-Edge', 
                        'Neg_Dots-Edge', 'Pos_Dots', 'Neg_Dots', 'Pos_Junctions3', 'Neg_Junctions3', 'Pos_Junctions4', 'Neg_Junctions4', 
                        'Pos_Junctionsx', 'Neg_Junctionsx', 'Pos_Defects', 'Neg_Defects', 'Total_Defects', 'Total_Area_nm', 
                        'Defect_Density_nm', 'Defect_Density_um']
                df1 = pd.DataFrame([[np.NaN for i in range(len(columns))]], columns=columns)
            
            print('Databese creation')  
            database_ADAblock = database_ADAblock.append(df1, ignore_index=True)
            database_python = database_python.append(df_python_data, ignore_index = True)
            database_process = database_process.append(df_process_param, ignore_index = True)
                    
            database = pd.concat([database_process,database_python,database_ADAblock], axis=1)
            database['Molecular_weight_Copolymer_g_over_mol'] = database['Molecular_weight_BCP_A_g_over_mol']+database['Molecular_weight_BCP_B_g_over_mol']+database['Molecular_weight_BCP_C_g_over_mol']
            database['scarto_percentuale'] = np.abs(100*((database['l0']-database['wfft_Period_nm'])/(database['l0']+database['wfft_Period_nm'])))
            #database = database.drop('Image_Name', axis=1)
            database['LER_su_d']= database['LER_sigma_avg_nm']/(database['d']*0.5)
            database['LWR_su_d']= database['LWR_sigma_avg_nm']/(database['d']*0.5)
            database['Terminals-Edge_Density_um']=0
            database['Terminals_Density_um']=(0.5*(database['Pos_Terminals']+database['Neg_Terminals'])/database['Total_Area_nm']*1000000)
            database['Dots-Edge_Density_um']=(0.5*(database['Pos_Dots-Edge']+database['Neg_Dots-Edge'])/database['Total_Area_nm']*1000000)
            database['Dots_Density_um']=((database['Pos_Dots']+database['Neg_Dots'])/database['Total_Area_nm']*1000000)
            database['Junctions3_Density_um']=(0.5*(database['Pos_Junctions3']+database['Neg_Junctions3'])/database['Total_Area_nm']*1000000)
            database['Junctions4_Density_um']=((database['Pos_Junctions4']+database['Neg_Junctions4'])/database['Total_Area_nm']*1000000)
            database['Junctionsx_Density_um']=(0.5*(database['Pos_Junctionsx']+database['Neg_Junctionsx'])/database['Total_Area_nm']*1000000)
            
            database.rename(columns={"nm_per_pixel": "nm_over_pixel", "InvRisoluzione": "InvResolution", "LER_su_d": "LER_over_d","LWR_su_d": "LWR_over_d","scarto_percentuale": "relative_uncertainty"}, inplace=True)

            database.to_csv(f'{indirizzo}analysis/database.csv', sep=';')

            shutil.move(os.path.join(indirizzo, file), f'{indirizzo}done')
            plt.cla()
            plt.clf()
            plt.close('all')
            matplotlib.use('Agg')
            print('Done!')

    print('Done for all images in the folder!')
    print('Database:')
    print(database)
    try: 
        print('Partial_database:')
        print(partial_database)
    except UnboundLocalError: 
        print('No images partially analzed')
    
main()  
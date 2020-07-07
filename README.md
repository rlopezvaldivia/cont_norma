# cont_norma
Python 3 code to interactivaly normalized a spectra. The flux is binned through a sigma-clipping process to fit a polynomial of order n and obtain a continuum fuction that is use to normalized the spectra. 

# Basic useage
It is requiered to create previously two folders: final and regions

          from cont_norma import *
          
          norma(file_list, file_folder, orden = 1, bins_ini = 20, ebins = 0, yof_ini = 0, valoff_ini = 0.005, reg_ini = 1, regs=[[15600,15650], [16720,16800], [22010,22130], [22170,22280],[22530,22730], [22920,23040]])

To see a simple example run: **python norma.py**
          
 # Description
   **file_list**: Coma-separated list with the names of the ascii files of the spectra to normalized and the name of the object.   
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Spectra needs to have the following columns: wavelength (angstroms), flux, snr.  
   **orden**: polynomial order to fit. Default value = 1.  
   **bins_ini** = number of flux bins. Default value = 20.  
   **ebins** = number of bins excluded. This excludes the ebins with the lower flux value. Default value = 0.  
   **yof_in** = vertical offset aplied to the normalized spectra. Default value = 0.  
   **valoff_ini** = increase for vertical offset. Default value = 0.005.   
   **reg_ini** = spectral regions to start with. Default = 1.  
   **regs** = spectral regions to normalized.
 
# Output
  In **regs** folder is saved each normalized spectral region as "obj_name-#.dat"   
  In **final** folder is saved in the file "obj_name_reg.dat".  
  **log.txt** saved the name of the object, polynomial order, number of bins, number of excluded bins, y offset, spectral region number and date.
    
    

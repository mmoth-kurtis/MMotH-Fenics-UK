## How to Generate Meshes for FEniCS from MRI Data
* Download MeVisLab, and appropriate files from the 'mesh_generation' directory.
* Load 'LV_Segmentation_UKY.mlab' in MeVisLab.
* Using the DicomImport module in MeVisLab, convert raw MRI images to DICOM format.
* Follow instructions from FEniCS_LV_Segmentation.pdf. This should yield a '.stl' file.
* Run the 1.createLV.py script on the stl file. This generates ??
* Run LV_Test.py on _something_ to yield the HDF5 file. This is what is needed for FEniCS. 

**Fill more information in as I work through the process**
